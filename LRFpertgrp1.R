library(mvtnorm)
library(MASS)
library(pracma)
library(expm)
OpenBlasThreads::set_num_threads(1)

invfact <- function(lg, n){
  ret <- diag(p)
  for(i in 1:n){
    ret <- ret + lg %^% i
  }
  return(ret)
}

matrixouter <- function(A, B, r, J){
  nr <- nrow(A)
  nc <- ncol(B)
  mat <- matrix(0, nr * nc, r*J)
  
  k   <- 1
  for(i in 1:nc){
    for(j in 1:nr){
      mat[k, ] <- outer(A[j, ], B[, i])
      k        <- k + 1
    }
  }
  return(mat)
}

matrixouter1 <- function(A, B, r, J){
  nr <- nrow(A)
  #nc <- ncol(B)
  mat <- matrix(0, nr, r*J)
  
  k   <- 1
  for(j in 1:nr){
    mat[k, ] <- outer(A[j, ], B)
    k        <- k + 1
  }
  return(mat)
}

QYpr <- function(i, mat = Y){
  temp <- matrix(Qlist[, ceiling(i/50)], p, p)
  return(temp%*%mat[, i])
}

set.seed(1)
# data generation
n <- 500
p  <- 21
r  <- 5
grp <- 10

eta0 <- matrix(rnorm(r*n), r, n) 

lambda0           <- matrix(0, p, r)
lambda0[1:5, 1]   <- rnorm(5, mean = 5, sd = 1)
lambda0[5:9, 2]   <- rnorm(5, mean = 5, sd = 1)
lambda0[9:13, 3]  <- rnorm(5, mean = 5, sd = 1)
lambda0[13:17, 4] <- rnorm(5, mean = 5, sd = 1)
lambda0[17:21, 5] <- rnorm(5, mean = 5, sd = 1)

image(t(lambda0))
Y <- matrix(rnorm(p*n, mean = lambda0 %*% eta0, sd = 1), p, n)#matrix(rnorm(p*n, mean = lambda0 %*% eta0, sd = 1), p, n) #t(rmvnorm(n, sigma = solve(pdmat)))#
#Plambda0 <- diag(p) - lambda0 %*% solve(crossprod(lambda0)) %*% t(lambda0)
Qlist <- matrix(array(diag(p)), p^2, grp)
po <- 1
Qmean <- array(1*diag(p))

set.seed(100)
eta0test <- matrix(rnorm(r*n), r, n)
Ytest <- matrix(rnorm(p*n, mean = lambda0 %*% eta0test, sd = 1), p, n)#matrix(rnorm(p*n, mean = lambda0 %*% eta0, sd = 1), p, n) #t(rmvnorm(n, sigma = solve(pdmat)))#
#Plambda0 <- diag(p) - lambda0 %*% solve(crossprod(lambda0)) %*% t(lambda0)
Qlist <- matrix(array(diag(p)), p^2, grp)
po <- 1
Qmean <- array(1*diag(p))
set.seed(1)

for(i in 2:grp){
  Qlist[, i] <- array(matrix(rnorm(p^2, Qmean, sd = sqrt(0.0001)), p, p))
}


Y <- parallel::mcmapply(1:n, FUN = QYpr, MoreArgs = list(mat=Y)) #matrix(Q %*% array(Y), p, n)

Ytest <- parallel::mcmapply(1:n, FUN = QYpr, MoreArgs = list(mat=Ytest))

pdmat0 <- lambda0 %*% t(lambda0) + diag(p)
seed <- 10
set.seed(seed)
RO <- randortho(p)
# parameters initialization
set.seed(seed)
r = p
nu <- 1
a1 <- 2
eta.var2 <- rep(0, r)
var <- eta.var2
delta <- rnorm(r)
eta.var2[r:1] <- rep(1, r)#cumsum(exp(delta[r:1])) / sum(exp(delta))
set.seed(seed)
phi <- matrix(rgamma(p * r, nu / 2, nu / 2), p, r)

psi0 <- c(rgamma(1,5,1), rgamma(r-1, 5, 1))
tau0 <- exp(cumsum(log(psi0)))

psi <- c(rgamma(1,5,1), rgamma(r-1, 5, 1))
tau <- exp(cumsum(log(psi)))

phi2 <- matrix(rgamma(p * r, nu / 2, nu / 2), p, r)

psi2 <- c(rgamma(1,5,1), rgamma(r-1, 5, 1))
tau2 <- exp(cumsum(log(psi2)))

lambda <- Y %*% ginv(crossprod(Y)) %*% t(Y)
eta <- ginv(crossprod(lambda)) %*% (crossprod(lambda, Y))
lambdaginv <- ginv(crossprod(lambda)) %*% t(lambda)

gamma <- lambdaginv 

sigma1 <- apply(Y - lambda %*% eta, 1, sd) 

sigma2 <- apply(eta - gamma %*% Y, 1, sd)

Yhat  <- lambda %*% eta
for(i in 1:p){
  al       <- 0.1 + n /2
  be       <- 0.1 + sum((Y[i, ]-Yhat[i, ])^2)/2
  sigma1[i] <- sqrt(1/rgamma(1, al, be))
}

etahat  <- gamma %*% Y

for(i in 1:r){
  al       <- 0.1 + n /2
  be       <- 0.1 + sum((eta[i, ]-etahat[i, ])^2)/2
  sigma2[i] <- sqrt(1/rgamma(1, al, be))
}

Total_itr <- 5000
eta_p   <- list()
lambda_p <- list()
sigma1_p <- list()
sigma2_p <- list()
itr <- 0
sds1 <- 1e-3
sds2 <- 1e-3
sdlam <- 1e-1
sdgam <- 1e-3
R <- 100
incre <- 4
secparam <- 10*(1:r)
ard1 <- 0
ard2 <- 0
arl <- 0
arg <- 0
L2 <- 10
epsilon2 <- t(rmvnorm(n, sigma = diag(sigma2)))
#Q <- diag(p*n)
po <- 1
Qmean <- array(1*diag(p))
Q_p <- list()
QV <- 0.01
QV_p <- rep(0, Total_itr)
Qiup <- function(i, QVf){
  Qtemp    <- matrix(Qlist[, i], p, p)
  index <- 1:50 + (i-1)*50
  for(j in 1:p){
    Qmean    <- rep(0, p)
    Qmean[j] <- 1
    mean.lami <- -(rowSums((Qtemp[, -j]%*%Y[-j, index] - lambda %*% eta[, index])*matrix(Y[j, index], p, 50, byrow = T))/sigma1^2) + Qmean/QVf
    var.lami  <- 1/(sum(Y[j, index]^2)/sigma1^2 + rep(1/QVf, p)) #phi[, i]*tau[i]) #rep(1/100, p)))
    mean.lami <- var.lami * mean.lami
    
    Qtemp[, j] <- rnorm(p, mean.lami, sqrt(var.lami))
  }
  
  return(Qtemp)
}
out <- array(0, Total_itr)
Qmeanmat <- matrix(array(1*diag(p)), p^2, grp)
while (itr < Total_itr) {
  itr <- itr + 1
  QY <- parallel::mcmapply(1:n, FUN = QYpr, MoreArgs = list(mat=Y)) #matrix(Q %*% array(Y), p, n)
  Yhatred <- QY - lambda %*% epsilon2
  
  for(i in 1:r){
    al       <- 0.1 + n /2
    be       <- 0.1 + sum((Yhatred[i, ])^2)/2
    sigma1[i] <- sqrt(1/rgamma(1, al, be))
  }
 
  for(i in 1:r){
    al       <- 100 + n /2
    be       <- 0.1 + sum((epsilon2[i, ])^2)/2
    sigma2[i] <- sqrt(1/rgamma(1, al, be))
  }

  var.pm <- ginv(crossprod(lambda / sigma1) + diag(1 / sigma2^2))
  var.pm <- (var.pm + t(var.pm)) / 2
  
  mean.etai <- t(lambda/sigma1^2) %*% (QY)
  mean.etai <- var.pm %*% mean.etai
  
  temp <- apply(mean.etai, 2, FUN=function(x){rmvnorm(1, x, var.pm)})
  epsilon2  <- temp
  
  sigma1_p[[itr]] <- sigma1
  sigma2_p[[itr]] <- sigma2
  
  eta <- epsilon2
  
  for(i in 1:r){
    mean.lami <- rowSums((QY - lambda[, -i] %*% eta[ - i, ])*matrix(eta[i, ], p, n, byrow=T)/sigma1^2)
    var.lami  <- 1/(sum(eta[i, ]^2)/sigma1^2 + phi[, i]*tau[i]) #phi[, i]*tau[i]) #rep(1/100, p)))
    var.lami  <- ( var.lami + t(var.lami) ) / 2
    mean.lami <- var.lami * mean.lami
    
    lambda[, i] <- rnorm(p, mean.lami, sqrt(var.lami))
    
    phi[, i]    <- rgamma(p, (nu + 1) / 2, ((lambda[, i]^2) * tau[i] + nu) / 2)
  }
  
  temp1 <- colSums((lambda ^ 2) * phi)
  
  tauprime1 <- tau / (psi[1])
  
  psi[1] <- rgamma(1, a1 + p*r/2, 1 + sum(temp1 * tauprime1) / 2)
  tau    <- tauprime1 * psi[1]
  
  for(i in 2:r){
    tauprime1 <- tau / (psi[i])
    psi[i] <- rgamma(1, a1 + p*(r - i + 1) / 2, 1 + sum(temp1[i:r] * tauprime1[i:r]) / 2)
    tau    <- tauprime1 * psi[i]
  }
  
  Qlist[, 2:grp] <- parallel::mcmapply(2:grp, FUN = Qiup, MoreArgs = list(QVf=QV))
  
  Qlistmat <- (Qlist)
  
  #update QV
  al       <- 0.1 + (grp-1)*p*p /2
  be       <- 0.1 + sum((Qlistmat-Qmeanmat)^2)/2
  QV       <- (1/rgamma(1, al, be))
  
  QV_p[itr] <- QV
 
  Q_p[[itr]] <- Qlist
  
  eta.var2[r:1] <- sigma2[r:1]^2
  
  lambda_p[[itr]] <- lambda
  
  if(itr %% R==0){
    u <- runif(1)
    if(u < exp(-1 - itr * 5 * 10^(-4) )){
      temp    <- colMeans(abs(lambda))
      c <- (which(temp < 1e-4))
      if(r-length(c)<3){c <- NULL}
      if(r < 3) {c <- NULL}
      if(length(c)> 0){
        r <- r - length(c)
        lambda <- lambda[, -c]
        phi <- phi[, -c]
        tau   <- tau[-c]
        psi   <- psi[-c]
        phi2 <- phi2[, -c]
        tau2   <- tau2[-c]
        psi2   <- psi2[-c]
        tau0   <- tau0[-c]
        psi0   <- psi0[-c]
        eta  <- eta[-c, ]
        epsilon2  <- epsilon2[-c, ]
        eta.var2 <- eta.var2[ -c]
        sigma2 <- eta.var2
        gamma    <- gamma[-c, ]
        delta    <- delta[ - c]
      }
    }
    R     <- R + incre
  }
  
  image(t(lambda))
  
  var <- lambda %*% diag(sigma2^2) %*% t(lambda) + diag(sigma1^2)
  Ytestn <- parallel::mcmapply(1:n, FUN = QYpr, MoreArgs = list(mat=Ytest))
  print(out[itr] <- sum(apply(Ytestn, 2, dmvnorm, sigma = var, log = T)))
}

save(out, lambda_p, sigma1_p, sigma2_p, Q_p, QV_p, file = paste("0001QVranperterb100LRFgrp",seed,".rda", sep = ""))