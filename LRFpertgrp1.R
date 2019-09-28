library(mvtnorm)
library(MASS)
library(pracma)
library(expm)
OpenBlasThreads::set_num_threads(1)
HMC = function (U, grad_U, epsilon, L = 30, current_q, ard)
{
  q = current_q
  p = rnorm(length(q),0,1)  # independent standard normal variates
  current_p = p
  
  # Make a half step for momentum at the beginning
  
  p = p - epsilon * grad_U(q) / 2
  
  # Alternate full steps for position and momentum
  
  for (i in 1:L)
  {
    # Make a full step for the position
    
    q = q + epsilon * p
    
    # Make a full step for the momentum, except at end of trajectory
    
    if (i!=L) p = p - epsilon * grad_U(q)
  }
  
  # Make a half step for momentum at the end.
  
  p = p - epsilon * grad_U(q) / 2
  
  # Negate momentum at end of trajectory to make the proposal symmetric
  
  p = -p
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  
  current_U = U(current_q)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q)
  proposed_K = sum(p^2) / 2
  
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  
  if (log(runif(1)) < (current_U-proposed_U+current_K-proposed_K))
  {
    ard <- ard + 1
    return (list(up=exp(q), ard = ard))  # accept
  }
  else
  {
    return (list(up=exp(current_q), ard = ard))  # reject
  }
}

U1 <- function(delta){
  cov <- diag(1/exp(delta)) %*% lambda + t(gamma) %*% diag(1/sigma2^2)
  invetavar <- (t(lambda)%*%diag(1/exp(delta))%*%lambda + diag(1/sigma2^2))
  invdatavar <- t(gamma)%*%diag(1/sigma2^2)%*%gamma + diag(1/exp(delta)) 
  pres1 <- rbind(cbind(invdatavar, cov), cbind(t(cov), invetavar))
  da <- rbind(Y, eta)
  part2 <- sapply(1:n, FUN = function(j){(t(da[, j]) %*% pres1 %*% da[, j])})
  part2 <- sum(unlist(part2))
  detp <- determinant(pres1)
  part1 <- (1/2)*detp$modulus*detp$sign
  ret <- - part1 + part2 / 2 + sum(delta^2)/(2 * 100)
  return(ret)
}

U1n <- function(delta){
  cov <- diag(1/(delta)) %*% lambda + t(gamma) %*% diag(1/sigma2^2)
  invetavar <- (t(lambda)%*%diag(1/(delta))%*%lambda + diag(1/sigma2^2))
  invdatavar <- t(gamma)%*%diag(1/sigma2^2)%*%gamma + diag(1/(delta)) 
  pres1 <- invdatavar - cov %*% ginv(invetavar) %*% t(cov)
  da <- Y
  part2 <- sapply(1:n, FUN = function(j){(t(da[, j]) %*% pres1 %*% da[, j])})
  part2 <- sum(unlist(part2))
  detp <- determinant(pres1)
  part1 <- (1/2)*detp$modulus*detp$sign
  ret <- - part1 + part2 / 2 - sum(dgamma(1/delta, 0.1, 0.1, log = T))
  return(-ret)
}

U1nn <- function(delta){
  cov <- diag(1/(delta)) %*% lambda + t(gamma) %*% diag(1/sigma2^2)
  invetavar <- (t(lambda)%*%diag(1/(delta))%*%lambda + diag(1/sigma2^2))
  invdatavar <- t(gamma)%*%diag(1/sigma2^2)%*%gamma + diag(1/(delta)) 
  pres1 <- rbind(cbind(invdatavar, cov), cbind(t(cov), invetavar))
  da <- rbind(Y, eta)
  part2 <- sapply(1:n, FUN = function(j){(t(da[, j]) %*% pres1 %*% da[, j])})
  part2 <- sum(unlist(part2))
  detp <- determinant(pres1)
  part1 <- (1/2)*detp$modulus*detp$sign
  ret <- - part1 + part2 / 2 - sum(dgamma(1/delta, 0.1, 0.1, log = T))
  return(-ret)
}



grad_U1 <- function(delta){
  cov <- diag(1/exp(delta)) %*% lambda + t(gamma) %*% diag(1/sigma2^2)
  invetavar <- (t(lambda)%*%diag(1/exp(delta))%*%lambda + diag(1/sigma2^2))
  invdatavar <- t(gamma)%*%diag(1/sigma2^2)%*%gamma + diag(1/exp(delta)) 
  pres1 <- rbind(cbind(invdatavar, cov), cbind(t(cov), invetavar))
  da <- rbind(Y, eta)
  covar <- ginv(pres1)
  l <- length(delta)
  ret <- rep(0, l)
  for(i in 1:l){
    ndelta <- rep(1, length(delta)) #exp(- delta)
    ndelta[ - i] <- 0
    cov <- diag(ndelta) %*% lambda 
    invetavar <- (t(lambda)%*%diag(ndelta)%*%lambda )
    invdatavar <- diag(ndelta) 
    presder <- rbind(cbind(invdatavar, cov), cbind(t(cov), invetavar))
    
    part2 <- sapply(1:n, FUN = function(j){(t(da[, j]) %*% presder %*% da[, j])})
    part2 <- - sum(unlist(part2)) * exp(-delta[i])
    detp <- sum(diag(covar %*% (presder* exp( - delta[i]))))
    part1 <-  - (1/2)*detp
    ret[i] <- - part1 + part2 / 2 + (delta[i])/(100)
  }
  return(ret)
}

U2 <- function(delta){
  cov <- diag(1/sigma1^2) %*% lambda + t(gamma) %*% diag(1/exp(delta))
  invetavar <- (t(lambda)%*%diag(1/sigma1^2)%*%lambda + diag(1/exp(delta)))
  invdatavar <- t(gamma)%*%diag(1/exp(delta))%*%gamma + diag(1/sigma1^2)
  
  pres1 <- rbind(cbind(invdatavar, cov), cbind(t(cov), invetavar))
  da <- rbind(Y, eta)
  part2 <- sapply(1:n, FUN = function(j){(t(da[, j]) %*% pres1 %*% da[, j])})
  part2 <- sum(unlist(part2))
  detp <- determinant(pres1)
  part1 <- (1/2)*detp$modulus*detp$sign
  ret <- - part1 + part2 / 2 + sum(delta^2)/(2 * 100)
  return(ret)
}

U2n <- function(delta){
  cov <- diag(1/sigma1^2) %*% lambda + t(gamma) %*% diag(1/(delta))
  invetavar <- (t(lambda)%*%diag(1/sigma1^2)%*%lambda + diag(1/(delta)))
  invdatavar <- t(gamma)%*%diag(1/(delta))%*%gamma + diag(1/sigma1^2)
  
  pres1 <- invdatavar - cov %*% ginv(invetavar) %*% t(cov)
  da <- Y
  part2 <- sapply(1:n, FUN = function(j){(t(da[, j]) %*% pres1 %*% da[, j])})
  part2 <- sum(unlist(part2))
  detp <- determinant(pres1)
  part1 <- (1/2)*detp$modulus*detp$sign
  ret <- - part1 + part2 / 2 - sum(dgamma(1/delta, 0.1, 0.1, log = T))
  return(-ret)
}

U2nn <- function(delta){
  cov <- diag(1/sigma1^2) %*% lambda + t(gamma) %*% diag(1/(delta))
  invetavar <- (t(lambda)%*%diag(1/sigma1^2)%*%lambda + diag(1/(delta)))
  invdatavar <- t(gamma)%*%diag(1/(delta))%*%gamma + diag(1/sigma1^2)
  
  pres1 <- rbind(cbind(invdatavar, cov), cbind(t(cov), invetavar))
  da <- rbind(Y, eta)
  part2 <- sapply(1:n, FUN = function(j){(t(da[, j]) %*% pres1 %*% da[, j])})
  part2 <- sum(unlist(part2))
  detp <- determinant(pres1)
  part1 <- (1/2)*detp$modulus*detp$sign
  ret <- - part1 + part2 / 2 - sum(dgamma(1/delta, 0.1, 0.1, log = T))
  return(-ret)
}

Ulam <- function(delta, i){
  lambdac <- lambda 
  lambdac[, i] <- delta
  cov <- diag(1/sigma1^2) %*% lambdac + t(gamma) %*% diag(1/sigma2^2)
  invetavar <- (t(lambdac)%*%diag(1/sigma1^2)%*%lambdac + diag(1/sigma2^2))
  invdatavar <- t(gamma)%*%diag(1/sigma2^2)%*%gamma + diag(1/sigma1^2) 
  pres1 <- rbind(cbind(invdatavar, cov), cbind(t(cov), invetavar))
  da <- rbind(Y, eta)
  part2 <- sapply(1:n, FUN = function(j){(t(da[, j]) %*% pres1 %*% da[, j])})
  part2 <- sum(unlist(part2))
  detp <- determinant(pres1)
  part1 <- (1/2)*detp$modulus*detp$sign
  ret <- - part1 + part2 / 2 + sum(delta^2*phi[, i]*tau[i] / 2)
  return(-ret)
}

Ugam <- function(delta, i){
  gammac <- gamma 
  gammac[, i] <- delta
  cov <- diag(1/sigma1^2) %*% lambda + t(gammac) %*% diag(1/sigma2^2)
  invetavar <- (t(lambda)%*%diag(1/sigma1^2)%*%lambda + diag(1/sigma2^2))
  invdatavar <- t(gammac)%*%diag(1/sigma2^2)%*%gammac + diag(1/sigma1^2) 
  pres1 <- rbind(cbind(invdatavar, cov), cbind(t(cov), invetavar))
  da <- rbind(Y, eta)
  part2 <- sapply(1:n, FUN = function(j){(t(da[, j]) %*% pres1 %*% da[, j])})
  part2 <- sum(unlist(part2))
  detp <- determinant(pres1)
  part1 <- (1/2)*detp$modulus*detp$sign
  ret <- - part1 + part2 / 2 + sum(delta^2*phi[i, ]*tau / 2) #*phi[i, ]*tau
  return(-ret)
}

Ulamn <- function(delta, i){
  lambdac <- lambda 
  lambdac[, i] <- delta
  cov <- diag(1/sigma1^2) %*% lambdac + t(gamma) %*% diag(1/sigma2^2)
  invetavar <- (t(lambdac)%*%diag(1/sigma1^2)%*%lambdac + diag(1/sigma2^2))
  invdatavar <- t(gamma)%*%diag(1/sigma2^2)%*%gamma + diag(1/sigma1^2) 
  pres1 <- invdatavar - cov %*% ginv(invetavar) %*% t(cov)
  da <- Y
  part2 <- sapply(1:n, FUN = function(j){(t(da[, j]) %*% pres1 %*% da[, j])})
  part2 <- sum(unlist(part2))
  detp <- determinant(pres1)
  part1 <- (1/2)*detp$modulus*detp$sign
  ret <- - part1 + part2 / 2 + sum(delta^2*phi[, i]*tau[i] / 2)
  return(-ret)
}

Ugamn <- function(delta, i){
  gammac <- gamma 
  gammac[, i] <- delta
  cov <- diag(1/sigma1^2) %*% lambda + t(gammac) %*% diag(1/sigma2^2)
  invetavar <- (t(lambda)%*%diag(1/sigma1^2)%*%lambda + diag(1/sigma2^2))
  invdatavar <- t(gammac)%*%diag(1/sigma2^2)%*%gammac + diag(1/sigma1^2) 
  pres1 <- invdatavar - cov %*% ginv(invetavar) %*% t(cov)
  da <- Y
  part2 <- sapply(1:n, FUN = function(j){(t(da[, j]) %*% pres1 %*% da[, j])})
  part2 <- sum(unlist(part2))
  detp <- determinant(pres1)
  part1 <- (1/2)*detp$modulus*detp$sign
  ret <- - part1 + part2 / 2 + sum(delta^2*phi[i, ]*tau / 2) #*phi[i, ]*tau
  return(-ret)
}

grad_U2 <- function(delta){
  cov <- diag(1/sigma1^2) %*% lambda + t(gamma) %*% diag(1/exp(delta))
  invetavar <- (t(lambda)%*%diag(1/sigma1^2)%*%lambda + diag(1/exp(delta)))
  invdatavar <- t(gamma)%*%diag(1/exp(delta))%*%gamma + diag(1/sigma1^2)
  
  pres1 <- rbind(cbind(invdatavar, cov), cbind(t(cov), invetavar))
  da <- rbind(Y, eta)
  covar <- ginv(pres1)
  l <- length(delta)
  ret <- rep(0, l)
  for(i in 1:l){
    ndelta <- rep(1, length(delta)) #exp(- delta)
    ndelta[ - i] <- 0
    cov <- t(gamma) %*% diag(ndelta)
    invetavar <-  diag(ndelta)
    invdatavar <- t(gamma)%*%diag(ndelta)%*%gamma
    presder <- rbind(cbind(invdatavar, cov), cbind(t(cov), invetavar))
    
    part2 <- sapply(1:n, FUN = function(j){(t(da[, j]) %*% presder %*% da[, j])})
    part2 <- - sum(unlist(part2)) * exp(-delta[i])
    detp <- sum(diag(covar %*% (presder* exp( - delta[i]))))
    part1 <-  - (1/2)*detp
    ret[i] <- - part1 + part2 / 2 + (delta[i])/(100)
  }
  return(ret)
}

Posdef <- function (n, ev = runif(n, 0, 3))
{
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp)
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

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

pdmat <- Posdef(n=p, ev=runif(p, 0, 10))

pdmat[which(abs(pdmat) < 1)] <- 0

eta0 <- matrix(rnorm(r*n), r, n) 
# mat1 <- matrix(rnorm(7*3, mean = 5, sd = 1),7, 3)
# mat2 <- matrix(rnorm(7*3, mean = 10, sd = 2),7, 3)
# mat3 <- matrix(rnorm(7*4, mean = 12, sd = 2),7, 4)
# 
# lambda0 <- as.matrix(bdiag(mat1, mat2, mat3))#matrix(rnorm(p*r, mean = 0, sd = 1),p, r)

lambda0           <- matrix(0, p, r)
lambda0[1:5, 1]   <- rnorm(5, mean = 5, sd = 1)
lambda0[5:9, 2]   <- rnorm(5, mean = 5, sd = 1)
lambda0[9:13, 3]  <- rnorm(5, mean = 5, sd = 1)
lambda0[13:17, 4] <- rnorm(5, mean = 5, sd = 1)
lambda0[17:21, 5] <- rnorm(5, mean = 5, sd = 1)

# lambda0[1:64, 1]                  <- rnorm(64, mean = 5, sd = 1)
# lambda0[c(1+0:29, 65+0:29), 2]       <- rnorm(60, mean = 5, sd = 1)
# temp <- c(0:12 + 1, 0:12+30, 0:12+65, 0:12+96)
# lambda0[temp, 3] <- rnorm(52, mean = 5, sd = 1)
# 
# temp <- c(0:4 + 1, 0:4 + 14, 0:4+30,0:4 + 43)
# temp <- c(temp, 0:4 + 65, 0:4 + 78, 0:4+96,0:4 + 109)
# lambda0[temp, 4] <- rnorm(40, mean = 5, sd = 1)
# 
# temp <- c(0:1 + 1, 0:1 + 6, 0:1+14, 0:1 + 19, 0:1 + 30, 0:1+35, 0:1+43, 0:1+48)
# temp <- c(temp, 0:1 + 65, 0:1 + 70, 0:1+78, 0:1 + 83, 0:1 + 96, 0:1+101, 0:1+109, 0:1+114)
# 
# lambda0[temp, 5] <- rnorm(32, mean = 5, sd = 1)

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

#lambda <- Plambda0%*%matrix(rnorm(p*r, 0, sd = 1), p, r)#* RO #1 / sqrt(phi * matrix(tau, p, r, byrow = T)

lambda <- Y %*% ginv(crossprod(Y)) %*% t(Y)
eta <- ginv(crossprod(lambda)) %*% (crossprod(lambda, Y))
lambdaginv <- ginv(crossprod(lambda)) %*% t(lambda)

gamma <- lambdaginv #matrix(0, r, p) #ginv(tcrossprod(Y)) %*% (tcrossprod(Y, eta))

sigma1 <- apply(Y - lambda %*% eta, 1, sd) #rep(1, p)#apply(Y - lambda %*% eta, 1, sd) #

sigma2 <- apply(eta - gamma %*% Y, 1, sd) #rep(1, p)#apply(eta - gamma %*% Y, 1, sd) #

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
  # temp1 <- matrixouter1(diag(p), Y[, i], p, p)
  # 
  # Yresid <- lambda %*% epsilon2[, i]
  # Phimat <- matrix(1/0.0001, p, p)#t(phi2) * matrix(tau2, r, p)#matrix(1/100, r, p)#t(phi2) * tau2
  # 
  # vargamma1  <- crossprod(temp1/((sigma1))) + diag(array(Phimat))
  # file <- try(vargamma1  <- ginv(vargamma1), silent = T)
  # if(class(file) == "try-error"){vargamma1  <- solve(vargamma1)}
  # vargamma1  <- (vargamma1 + t(vargamma1)) / 2
  # meangamma1 <- array(t(temp1) %*% (array(Yresid) / ((sigma1) ^ 2))) + Qmean/0.0001
  # gen        <- mvtnorm::rmvnorm(1, vargamma1 %*% meangamma1, vargamma1)
  #index    <- p*(i-1) + 1:p
  Qtemp    <- matrix(Qlist[, i], p, p)#Q[index, index]
  index <- 1:50 + (i-1)*50
  for(j in 1:p){
    Qmean    <- rep(0, p)
    Qmean[j] <- 1
    mean.lami <- -(rowSums((Qtemp[, -j]%*%Y[-j, index] - lambda %*% eta[, index])*matrix(Y[j, index], p, 50, byrow = T))/sigma1^2) + Qmean/QVf
    var.lami  <- 1/(sum(Y[j, index]^2)/sigma1^2 + rep(1/QVf, p)) #phi[, i]*tau[i]) #rep(1/100, p)))
    #var.lami  <- ( var.lami + t(var.lami) ) / 2
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
  
  # M <- (Y - solve(diag(p)-lambda %*% gamma) %*% lambda %*% epsilon2)
  # al <- 0.1 + n * p / 2
  # be <- 0.1 + sum(sapply(1:n, FUN = function(j){t(M[, j])%*%((diag(p)-lambda %*% gamma)%^%2)%*%(M[, j])})) / 2
  # if(be > 0){
  #   gen <- sqrt(1/rgamma(1, al, be))
  #   
  #   sigma1 <- rep(gen, p)
  # }
  # 
  
  # temp1 <- rowSums((epsilon2 ^ 2))
  # 
  # tauprime1 <- tau0 / (psi0[1])
  # 
  # psi0[1] <- rgamma(1, 1 + n*r/2, 1 + sum(temp1 * tauprime1) / 2)
  # tau0    <- tauprime1 * psi0[1]
  # 
  # for(i in 2:r){
  #   tauprime1 <- tau0 / (psi0[i])
  #   psi0[i] <- rgamma(1, 2 + n*(r - i + 1) / 2, 1 + sum(temp1[i:r] * tauprime1[i:r]) / 2)
  #   tau0    <- tauprime1 * psi0[i]
  # }
  # 
  # sigma2 <- sqrt(1/tau0)
  
  #print(sigma2)
  
  for(i in 1:r){
    al       <- 100 + n /2
    be       <- 0.1 + sum((epsilon2[i, ])^2)/2
    sigma2[i] <- sqrt(1/rgamma(1, al, be))
  }
  
  #print(max(sigma1))
  #print(max(sigma2))
  
  # cov <- diag(1/sigma1^2) %*% lambda + t(gamma) %*% diag(1/sigma2^2)
  # invetavar <- ginv(t(lambda)%*%diag(1/sigma1^2)%*%lambda + diag(1/sigma2^2))
  # 
  # new <- t(gamma)%*%diag(1/sigma2^2)%*%gamma + diag(1/sigma1^2) - cov %*% invetavar %*% t(cov) 
  
  
  
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
    #var.lami  <- ( var.lami + t(var.lami) ) / 2
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
  
  # for(i in 1:n){
  #   
  #   index      <- 1:p + l
  #   Q[index, index]      <-  matrix(out[, i], p, p)#diag(p) #
  #   l <- l + p
  # }
  
  Q_p[[itr]] <- Qlist
  
  # for(i in 1:p){
  #   mean.gami <- Y %*% array(lambda[i, ] %*% epsilon2) / sigma1[i]^2
  #   var.gami  <- ginv(tcrossprod(Y)/sigma1[i]^2 + diag(rep(1/100, p))) #phi[, i]*tau[i]) #rep(1/100, p))
  #   var.gami  <- ( var.gami + t(var.gami) ) / 2
  #   mean.gami <- var.gami %*% mean.gami
  #   
  #   Q[i, ] <- rmvnorm(1, mean.lami, var.lami)
  #   
  #   #phi[, i]    <- rgamma(p, (nu + 1) / 2, ((lambda[, i]^2) * tau[i] + nu) / 2)
  # }
  
  # temp1 <- colSums((t(gamma)^2) * phi2)
  # 
  # tauprime1 <- tau2 / (psi2[1])
  # 
  # psi2[1] <- rgamma(1, a1 + p*r/2, 1 + sum(temp1 * tauprime1) / 2)
  # tau2    <- tauprime1 * psi2[1]
  # 
  # for(i in 2:r){
  #   tauprime1 <- tau / (psi2[i])
  #   psi2[i] <- rgamma(1, a1 + p*(r - i + 1) / 2, 1 + sum(temp1[i:r] * tauprime1[i:r]) / 2)
  #   tau2    <- tauprime1 * psi2[i]
  # }
  
  # file <- try(temp <- HMC(U, grad_U, epsilon = sdl, L = L2, delta, ard), silent = T)
  # if (class(file) == "try-error"){
  #   file <- try(temp <- HMC(U, grad_U, epsilon = sdl, L = 1, delta, ard), silent = T)
  #   if (class(file) == "try-error"){
  #     temp <- list(up = delta, ard = ard)
  #   }
  # }
  # delta <- temp$up
  # ard  <- temp$ard
  
  eta.var2[r:1] <- sigma2[r:1]^2#cumsum(exp(delta[r:1])) / sum(exp(delta))
  
  #sigma2_p[[itr]] <- sqrt(eta.var2)
  
  #eta_p[[itr]]  <- eta
  lambda_p[[itr]] <- lambda
  #gamma_p[[itr]] <- gamma
  
  # if(itr%%100 == 0){
  #   Plam <- lambda_p[[itr-99]] %*% ginv(crossprod(lambda_p[[itr-99]])) %*% t(lambda_p[[itr-99]])
  #   lamtemp <- Plam %*% lambda_p[[itr]]
  #   M <- crossprod(lambda_p[[itr-99]], lamtemp)
  #   sv <- svd(M)
  #   Dg <- sv$v %*% t(sv$u)
  #   Dg <- Dg[, 1:nrow(Dg)]
  #   print(mean((Dg-diag(nrow(Dg)))^2))
  #   
  #   
  #   new1 <- lambda %*% diag(sigma2^2) %*% t(lambda) + diag(sigma1^2)
  #   #print(mean((ginv(new)-solve(pdmat))^2))
  #   
  #   #lg <- solve(Q)
  #   
  #   #print(mean((((lg %*% (new1) %*% t(lg)))-(pdmat0))^2))
  #   #print(sigma2)
  # }
  
  #print(itr)
  
  if(itr %% R==0){
    u <- runif(1)
    if(u < exp(-1 - itr * 5 * 10^(-4) )){
      temp    <- colMeans(abs(lambda))
      #print(temp)
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
    #incre <- 2*incre
  }
  #eta.var2 <- sigma2^2
  image(t(lambda))
  #print(itr)
  
  var <- lambda %*% diag(sigma2^2) %*% t(lambda) + diag(sigma1^2)
  Ytestn <- parallel::mcmapply(1:n, FUN = QYpr, MoreArgs = list(mat=Ytest))
  print(out[itr] <- sum(apply(Ytestn, 2, dmvnorm, sigma = var, log = T)))
}

#gamma <- ginv(crossprod(lambda)) %*% t(lambda) %*% (diag(p) - Q) 

#lambda %*% gamma
#save(eta_p, lambda_p, sigma1_p, sigma2_p, Q_p, file = paste("perterb10LRF",seed,".rda", sep = ""))
#save(lambda_p, file = paste("perterbLRF",seed,".rda", sep = ""))

save(out, lambda_p, sigma1_p, sigma2_p, Q_p, QV_p, file = paste("0001QVranperterb100LRFgrp",seed,".rda", sep = ""))

# s2 <- NULL
# for(i in 1:5000){
#   s2 <- c(s2, max(sigma2_p[[i]]))
# }
# 
# plot(s2[2500:5000])

# s1 <- NULL
# for(i in 1:5000){
#   s1 <- c(s1, max(sigma1_p[[i]]))
# }
# 
# plot(s1[2500:5000])
# 
# 
# pres <- list()
# ind <- 0
# for(i in 2501:5000){
#   lambda <- lambda_p[[i]]
#   sigma1 <- sigma1_p[[i]]
#   sigma2 <- sigma2_p[[i]]
#   gamma <- gamma_p[[i]]
#   cov <- diag(1/sigma1^2) %*% lambda + t(gamma) %*% diag(1/sigma2^2)
#   invetavar <- solve(t(lambda)%*%diag(1/sigma1^2)%*%lambda + diag(1/sigma2^2))
#   pres[[i-2500]] <- t(gamma)%*%diag(1/sigma2^2)%*%gamma + diag(1/sigma1^2) - cov %*% invetavar %*% t(cov) 
# }
# 
# varesim <- list()
# lg <- list()
# lgi <- list()
# ind <- 0
# for(i in 2501:5000){
#   lambda <- lambda_p[[i]]
#   sigma1 <- sigma1_p[[i]]
#   sigma2 <- sigma2_p[[i]]
#   gamma <- gamma_p[[i]]
#   r <- ncol(lambda)
#   varest <- solve(diag(r)-gamma %*% lambda)%*%(gamma%*%diag(sigma1^2)%*%t(gamma) + diag(sigma2^2))%*%solve(diag(r)-gamma %*% lambda)
#   #pres[[i-2500]] <- (diag(p)-lambda %*% gamma)%*%solve(lambda%*%diag(sigma2^2)%*%t(lambda) + diag(sigma1^2))%*%(diag(p)-lambda %*% gamma)
#   new1 <- lambda %*% diag(sigma2^2) %*% t(lambda) + diag(sigma1^2)
#   #print(mean((ginv(new)-solve(pdmat))^2))
#   
#   lgi[[i-2500]] <- ginv(diag(p) - lambda %*% gamma)
#   lg[[i-2500]] <- (diag(p) - lambda %*% gamma)
#   varesim[[i -2500]] <- new1
# }
# 
# varesimp <- Reduce('+', varesim) / 2500
# 
# lgp <- Reduce('+', lg) / 2500 
# 
# #lambdap <- Reduce('+', lambda_p[2501:5000]) / 2500
# #gammap <- Reduce('+', gamma_p[2501:5000]) / 2500
# 
# #lgp <- ginv(diag(p) - lambdap %*% gammap)
# 
# mean((solve(lgp) %*% varesimp %*% solve(t(lgp)) - solve(pdmat))^2)
# 
# mean((varesimp - lgp %*% solve(pdmat) %*% t(lgp))^2)

varesim <- list()
pres <- list()

for(i in 2501:5000){
  lambda <- lambda_p[[i]]
  sigma1 <- sigma1_p[[i]]
  sigma2 <- sigma2_p[[i]]
  Q <- Q_p[[i]]
  new1 <- lambda %*% diag(sigma2^2) %*% t(lambda) + diag(sigma1^2)
  lg <- solve(Q)
  
  varesim[[i-2500]] <- lg %*% (new1) %*% t(lg)
  
  pres[[i-2500]] <- Q %*% solve(new1) %*% t(Q)
}
varesimp <- Reduce('+', varesim) / 2500

mean((varesimp-solve(pdmat))^2)

# presp <- Reduce('+', pres) / 2500
# 
# truepres <- pdmat #solve(lambda0 %*% t(lambda0) + diag(p))
# 
# index <- which(truepres != 0)
# 
# mean((presp[index]-truepres[index])^2)
# 
# mean((presp-truepres)^2)

robcor <- list()

for(i in 1:4999){
  M <- crossprod(lambda_p[[i]], lambda_p[[i+1]])
  sv <- svd(M)
  Dg <- sv$v %*% t(sv$u)
  Dg <- Dg[, 1:nrow(Dg)]
  robcor[[i]] <- mean((Dg-diag(nrow(Dg)))^2)#min(cancor(roblam_p1[[i]], roblam_p1[[i+1000]])$cor)#
}

robcorm <- unlist(robcor)#unlist(lapply(robcor, FUN=function(x){diff(range(x))}))
plot(robcorm[2500:4999])
mean(robcorm[2500:4999])

l <- 0
Bigmat <- matrix(0, p*(p+1)/2, 2500)


for(i in 1:p){
  for(j in i:p){
    temp <- rep(0, 2500)
    for(k in 1:2500){
      mat <- pres[[k]]
      temp[k] <- mat[i, j]
    }
    l <- l+1
    Bigmat[l, ] <- temp
  }
}

mse <- rep(0, 2500)
for(i in 1:2500){
  mse[i] <- mean((varesim[[i]]-solve(pdmat))^2)
}

mean(mse)

plot(apply(Bigmat, 1, sd))
# 
# #roblam1 <- matsig(lambda_p)
# load("robLVM1.rda")
# roblam_p1 <- lambda_p
# robeta_p1 <- eta_p
# load("robLVM10.rda")
# roblam_p10 <- lambda_p
# robeta_p10 <- eta_p
# 
# mean((roblam10-roblam1)^2)
# 
# #roblam10 <- matsig(lambda_p)
# load("robLVM100.rda")
# #roblam100 <- matsig(lambda_p)
# 
# 
# #cancor(t(lambda_p[[3000]]), t(lambda_p[[4000]]))$cor
# #lam1 <- matsig(lambda_p)
# load("LVM10.rda")
# lam_p10 <- lambda_p
# eta_p10 <- eta_p
# #lam10 <- matsig(lambda_p)
# 
# 
# load("robLVM10.rda")
# roblam_p1 <- lambda_p
# robeta_p1 <- eta_p
# load("LVM10.rda")
# lam_p1 <- lambda_p
# eta_p1 <- eta_p
# robcor <- list()
# 
# cor <- list()
# 
# for(i in 1:2999){
#   M <- crossprod(roblam_p1[[i]], roblam_p1[[i+1000]])
#   sv <- svd(M)
#   Dg <- sv$v %*% t(sv$u)
#   robcor[[i]] <- min(cancor(roblam_p1[[i]], roblam_p1[[i+1000]])$cor)#mean((Dg-diag(p))^2)
#   M <- crossprod(lam_p1[[i]], lam_p1[[i+1000]])
#   sv <- svd(M)
#   Dg <- sv$v %*% t(sv$u)
#   cor[[i]] <- min(cancor(lam_p1[[i]], lam_p1[[i+1000]])$cor)#mean((Dg-diag(p))^2)
# }
# 
# robcorm <- unlist(robcor)#unlist(lapply(robcor, FUN=function(x){diff(range(x))}))
# corm <- unlist(cor)#unlist(lapply(cor, FUN=function(x){diff(range(x))}))
# plot(robcorm)
# plot(corm)