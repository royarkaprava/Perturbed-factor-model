<<<<<<< HEAD
#'set.seed(1)
# data generation
#'

QYprG <- function(i, mat = Y, vec = grpind, Ql = Qlist){
  p <- nrow(Y)
  temp <- matrix(Ql[, vec[i]], p, p)
  return(temp%*%mat[, i])
}
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

Y <- matrix(rnorm(p*n, mean = lambda0 %*% eta0, sd = 1), p, n)#matrix(rnorm(p*n, mean = lambda0 %*% eta0, sd = 1), p, n)
Qlist <- matrix(array(diag(p)), p^2, grp)
Qmean <- array(1*diag(p))

set.seed(100)
Qlist <- matrix(array(diag(p)), p^2, grp)
Qmean <- array(1*diag(p))
set.seed(1)

#Generating the perturbation matrices
for(i in 2:grp){
  Qlist[, i] <- array(matrix(rnorm(p^2, Qmean, sd = sqrt(0.0001)), p, p))
}
grpind = rep(1:10, each = 50)
Y <- parallel::mcmapply(1:n, FUN = QYprG, MoreArgs = list(mat=Y)) #matrix(Q %*% array(Y), p, n)
fit <- PFA(Y, grpind = rep(1:10, each = 50), ini.PCA = T)
sigma2p <- Reduce('+', fit$Latentsigma[1:1000])/length(fit$Latentsigma[1:1000])
lambdap <- (Reduce('+', fit$Loading[1:1000])/length(fit$Latentsigma[1:1000])) %*% diag(sigma2p)
fields::image.plot(t(lambdap))

####Posterior samples of Shared variance
sharevar <- list()
for(i in 1:1000){
  sharevar[[i]] <- fit$Loading[i]%%diag(fit$Latentsigma[i]^2)%%t(fit$Loading[i]) + diag(fit$Errorsigma[i]^2)
}

######Posterior samples of Study specific variances
grpind = rep(1:10, each = 50)
QSpr <- function(i, mat, vec = grpind, Ql = Qlist){
  temp <- solve(matrix(Ql[, vec[i]], p, p))
  return(temp%*%mat%*%t(temp))
}

vars <- list()
for(i in 1:1000){
  Qlist <- fit$Pertmat[[i]]
  
  vars[[i]] <- parallel::mcmapply(1:n, FUN = QSpr, MoreArgs = list(mat=sharevar[[i]], vec = grpind, Ql = Qlist))
  
=======
#'set.seed(1)
# data generation
#'

QYprG <- function(i, mat = Y, vec = grpind, Ql = Qlist){
  p <- nrow(Y)
  temp <- matrix(Ql[, vec[i]], p, p)
  return(temp%*%mat[, i])
}
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

Y <- matrix(rnorm(p*n, mean = lambda0 %*% eta0, sd = 1), p, n)#matrix(rnorm(p*n, mean = lambda0 %*% eta0, sd = 1), p, n)
Qlist <- matrix(array(diag(p)), p^2, grp)
Qmean <- array(1*diag(p))

set.seed(100)
Qlist <- matrix(array(diag(p)), p^2, grp)
Qmean <- array(1*diag(p))
set.seed(1)

#Generating the perturbation matrices
for(i in 2:grp){
  Qlist[, i] <- array(matrix(rnorm(p^2, Qmean, sd = sqrt(0.0001)), p, p))
}
grpind = rep(1:10, each = 50)
Y <- parallel::mcmapply(1:n, FUN = QYprG, MoreArgs = list(mat=Y)) #matrix(Q %*% array(Y), p, n)
fit <- PFA(Y, grpind = rep(1:10, each = 50), ini.PCA = T)
sigma2p <- Reduce('+', fit$Latentsigma[1:1000])/length(fit$Latentsigma[1:1000])
lambdap <- (Reduce('+', fit$Loading[1:1000])/length(fit$Latentsigma[1:1000])) %*% diag(sigma2p)
fields::image.plot(t(lambdap))

####Posterior samples of Shared variance
sharevar <- list()
for(i in 1:1000){
  sharevar[[i]] <- fit$Loading[i]%%diag(fit$Latentsigma[i]^2)%%t(fit$Loading[i]) + diag(fit$Errorsigma[i]^2)
}

######Posterior samples of Study specific variances
grpind = rep(1:10, each = 50)
QSpr <- function(i, mat, vec = grpind, Ql = Qlist){
  temp <- solve(matrix(Ql[, vec[i]], p, p))
  return(temp%*%mat%*%t(temp))
}

vars <- list()
for(i in 1:1000){
  Qlist <- fit$Pertmat[[i]]
  
  vars[[i]] <- parallel::mcmapply(1:n, FUN = QSpr, MoreArgs = list(mat=sharevar[[i]], vec = grpind, Ql = Qlist))
  
>>>>>>> e7070f62bad2d869fcc62f10e783c756816c5f84
}