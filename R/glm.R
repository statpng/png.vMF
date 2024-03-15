#' @export Cq
Cq <- function(theta){
  q <- length(theta)
  NORM <- norm(theta, "2")
  
  # (NORM)^(q/2-1) / ( (2*pi)^(q/2) * besselI(NORM, q/2-1) )
  
  exp(  (q/2-1) * log(NORM) - ( (q/2) * log(2*pi) + log( besselI(NORM, q/2-1, expon.scaled=TRUE) ) + NORM )  )
}


#' @export logCq
logCq <- function(theta){
  q <- length(theta)
  NORM <- norm(theta, "2")
  
  # (NORM)^(q/2-1) / ( (2*pi)^(q/2) * besselI(NORM, q/2-1) )
  
  (q/2-1) * log(NORM) - ( (q/2) * log(2*pi) + log( besselI(NORM, q/2-1, expon.scaled=TRUE) ) + NORM )
}


#' @export Bq
Bq <- function(theta){
  q <- length(theta)
  NORM <- norm(theta, "2")
  
  # besselI(NORM, nu=q/2, expon.scaled=FALSE) == besselI(NORM, nu=q/2, expon.scaled=TRUE)/exp(-NORM)
  # besselI(NORM, nu=q/2-1, expon.scaled=FALSE) == besselI(NORM, nu=q/2-1, expon.scaled=TRUE)/exp(-NORM)
  
  besselI(NORM, nu=q/2, expon.scaled=TRUE) / besselI(NORM, nu=q/2-1, expon.scaled=TRUE)
}


#' @export Hq
Hq <- function(theta){
  
  q <- length(theta)
  NORM <- norm(theta, "2")
  
  I0 <- besselI(NORM, nu=q/2-1, expon.scaled=TRUE)
  I1 <- besselI(NORM, nu=q/2, expon.scaled=TRUE)

  1 - (I1/I0)^2 - (q-1)/NORM * I1/I0
  # numerator <- I0^2 - I1^2 - (q-1)/(NORM) * I0 * I1
  # denumerator <- I0^2
  # 
  # numerator / denumerator
}


#' @export subgrad
subgrad <- function(theta){
  q <- length(theta)
  NORM <- norm(theta, "2")
  
  if(NORM == 0){
    v <- runif(q,-1,1)
    v / (norm(v, "2")*1.1)
  } else {
    theta / NORM
  }
}


#' @export b1.vMF
b1.vMF <- function(theta){
  Bq(theta) * subgrad(theta)
}


#' @export b2.vMF
b2.vMF <- function(theta){
  q <- length(theta)
  NORM <- norm(theta, "2")
  s <- subgrad(theta)
  
  tcrossprod(s) * Hq(theta) + Bq(theta) * 1/NORM * (diag(1,q) - tcrossprod(theta/NORM) )
}


if(FALSE){
  set.seed(1)
  theta <- rnorm(3)
  
  b1(theta)
  
  Cq(theta)
}




glm_vmf <- function(X, Y, beta0=NULL, maxit=1000, eps=1e-8, standardize=TRUE){
  library(magic)
  
  if(FALSE){
    ii <- pvec.list[[d]]
    X=U.joint[,-1]
    Y=X_joint[,ii]
  }
  
  if(FALSE){
    X=simdata$X
    Y=simdata$Y
    maxit=100
    eps=1e-8
    standardize=TRUE
  }
  
  
  
  
  
  n <- nrow(Y)
  q <- ncol(Y)
  p <- ncol(X)
  
  X0 <- scale(X, center=TRUE, scale=FALSE)
  
  if(standardize){
    sdx.inv <- apply(X0,2,sd) %>% {diag(1/., length(.), length(.))}
    X <- X0 %*% sdx.inv
  } else {
    X <- X0
  }
  
  
  MU <- FrechetMean(Y)
  
  mu0 <- colSums(Y) %>% {./norm(.,"2")}
  rbar <- colSums(Y) %>% {norm(.,"2")/n}
  kappa0 <- (rbar*q - rbar^3) / (1-rbar^2)
  
  MU <- mu0 * kappa0
  
  
  
  beta <- lm(Y ~ -1 + X, offset=tcrossprod(rep(1,n),MU)) %>% coef %>% as.vector
  beta <- c(MU, beta)
  # beta[1:3] <- beta[1:3] %>% {./norm(.,"2")}
  X <- cbind(1, X)
  MU <- 0
  
  
  
  l1 <- 1
  beta_old <- beta
  beta_new <- beta + 1
  beta.list <- loglik.list <- crit.list <- NULL
  crit <- 1
  while( crit > eps & l1 <= maxit ){
    
    beta.list[[l1]] <- beta
    
    beta_old <- beta
    
    Xt.list <- apply(X, 1, function(Xi) kronecker(diag(1,q), Xi), simplify=FALSE)
    Xt <- do.call("cbind", Xt.list)
    yt <- as.vector(t(Y))
    
    b1.list <- lapply(1:n, function(i) b1.vMF( MU + crossprod(Xt.list[[i]], beta) ))
    eta <- do.call("rbind", b1.list)
      
    b2.list <- lapply(1:n, function(i) b2.vMF( MU + crossprod(Xt.list[[i]], beta) ))
    W <- do.call("adiag", b2.list) # magic::adiag
    
    Z <- crossprod(Xt, beta) + solve(W) %*% (yt - eta)
    
    beta <- solve(Xt %*% W %*% t(Xt)) %*% Xt %*% W %*% Z
    
    # beta[1:3] <- beta[1:3] %>% {./norm(.,"2")}
    
    beta_new <- beta
    
    loglik <- Reduce("+", lapply(1:n, function(i){
      theta <- MU + crossprod(Xt.list[[i]], beta)
      crossprod(theta, Y[i,]) + log( Cq(theta) )
    }))
    
    crit <- norm(beta_old-beta_new, "2")
    
    loglik.list[l1] <- loglik
    crit.list[l1] <- crit
    
    l1 <- l1 + 1
  }
  
  # print( matrix(beta,p,q) )
  # print( l1 )
  
  
  # plot(loglik.list, type="l")
  
  beta.list <- lapply(beta.list, function(beta) matrix(beta,p+1,q,byrow=FALSE) )
  beta_new <- matrix(beta_new,p+1,q,byrow=TRUE)
  
  if(standardize){
    beta.list <- lapply(beta.list, function(b) diag(c(1,sdx.inv)) %*% b)
    beta_new <- diag(c(1,sdx.inv)) %*% beta_new
  }
  
  
  
  list(mu=MU, 
       beta=beta_new,
       beta.list=beta.list,
       # beta=matrix(beta,p,q,byrow=TRUE),
       # beta.list=beta.list %>% lapply(function(x) matrix(x,p,q,byrow=TRUE)),
       # beta=matrix(beta,p+1,q,byrow=TRUE),
       # beta.list=beta.list %>% lapply(function(x) matrix(x,p+1,q,byrow=TRUE)),
       loglik.list=loglik.list[-1], crit.list=crit.list[-1])
  
}










glm_vmf_FixedMean <- function(X, Y, lambda=0, beta0=NULL, maxit=1000, eps=1e-8, standardize=TRUE){
  
  if(FALSE){
    ii <- pvec.list[[d]]
    X=U.joint[,-1]
    Y=X_joint[,ii]
  }
  
  
  
  library(magic)
  
  n <- nrow(Y)
  q <- ncol(Y)
  p <- ncol(X)
  
  X0 <- scale(X, center=TRUE, scale=FALSE)
  
  if(standardize){
    sdx.inv <- apply(X0,2,sd) %>% {diag(1/., length(.), length(.))}
    X <- X0 %*% sdx.inv
  } else {
    X <- X0
  }
  
  
  MU <- FrechetMean(Y)
  
  mu0 <- colSums(Y) %>% {./norm(.,"2")}
  rbar <- colSums(Y) %>% {norm(.,"2")/n}
  kappa0 <- (rbar*q - rbar^3) / (1-rbar^2)
  
  MU <- mu0 * kappa0
  
  
  # 1
  # {
  #   beta <- lm(Y ~ -1 + X, offset=tcrossprod(rep(1,n),MU)) %>% coef %>% as.vector
  #   beta <- c(MU, beta)
  #   # beta[1:3] <- beta[1:3] %>% {./norm(.,"2")}
  #   X <- cbind(1, X)
  #   MU <- 0
  # }
  
  
  # 2
  {
    beta <- lm(Y ~ -1 + X, offset=tcrossprod(rep(1,n),MU)) %>% coef %>% as.vector
    # beta <- lm(Y ~ 1 + X) %>% coef %>% {.[-1,]} %>% as.vector
  }
  
  
  
  
  
  l1 <- 1
  beta_old <- beta
  beta_new <- beta + 1
  beta.list <- loglik.list <- crit.list <- NULL
  crit <- 1
  while( crit > eps & l1 <= maxit ){
    
    beta.list[[l1]] <- beta
    
    beta_old <- beta
    
    Xt.list <- apply(X, 1, function(Xi) kronecker(diag(1,q), Xi), simplify=FALSE)
    Xt <- do.call("cbind", Xt.list)
    yt <- as.vector(t(Y))
    
    b1.list <- lapply(1:n, function(i) b1.vMF( MU + crossprod(Xt.list[[i]], beta) ))
    eta <- do.call("rbind", b1.list)
    
    b2.list <- lapply(1:n, function(i) b2.vMF( MU + crossprod(Xt.list[[i]], beta) ))
    W <- do.call("adiag", b2.list) # magic::adiag
    
    Z <- crossprod(Xt, beta) + solve(W) %*% (yt - eta)
    
    beta <- solve(Xt %*% W %*% t(Xt) + diag(lambda,p,p)) %*% Xt %*% W %*% Z
    
    # beta[1:3] <- beta[1:3] %>% {./norm(.,"2")}
    
    beta_new <- beta
    
    loglik <- Reduce("+", lapply(1:n, function(i){
      theta <- MU + crossprod(Xt.list[[i]], beta)
      crossprod(theta, Y[i,]) + log( Cq(theta) )
    }))
    
    crit <- norm(beta_old-beta_new, "2")
    
    loglik.list[l1] <- loglik
    crit.list[l1] <- crit
    
    l1 <- l1 + 1
  }
  
  # print( matrix(beta,p,q) )
  # print( l1 )
  
  
  # plot(loglik.list, type="l")
  
  beta.list <- lapply(beta.list, function(b) matrix(beta,p,q,byrow=FALSE) )
  beta_new <- matrix(beta_new,p,q,byrow=FALSE)
  
  if(standardize){
    beta.list <- lapply(beta.list, function(b) sdx.inv %*% b)
    beta_new <- sdx.inv %*% beta_new
  }
  
  
  
  list(mu=MU, 
       beta=beta_new,
       beta.list=beta.list,
       # beta=matrix(beta,p,q,byrow=TRUE),
       # beta.list=beta.list %>% lapply(function(x) matrix(x,p,q,byrow=TRUE)),
       # beta=matrix(beta,p+1,q,byrow=TRUE),
       # beta.list=beta.list %>% lapply(function(x) matrix(x,p+1,q,byrow=TRUE)),
       loglik.list=loglik.list[-1], crit.list=crit.list[-1])
  
}










glm_vmf_FixedMean_Offset <- function(X, Y, MU=NULL, Offset=NULL, lambda=1e-6, beta0=NULL, maxit=1000, eps=1e-8, standardize=TRUE){
  
  if(FALSE){
    X=simdata$X; Y=simdata$Y; tau=tau; Offset=NULL; lambda=0; maxit=100; eps=1e-8; standardize=TRUE
  }
  
  
  if(FALSE){
    ii <- pvec.list[[d]]
    X=U.joint[,-1]
    Y=X_joint[,ii]
    Offset=Theta.inds[,ii]
    
    beta0=NULL; maxit=1000; eps=1e-8; standardize=TRUE
  }
  
  
  
  library(magic)
  
  n <- nrow(Y)
  q <- ncol(Y)
  p <- ncol(X)
  
  if(is.null(Offset)){
    Offset <- matrix(0,n,q)
  }
  
  X0 <- scale(X, center=TRUE, scale=FALSE)
  
  if(standardize){
    sdx.inv <- apply(X0,2,sd) %>% {diag(1/., length(.), length(.))}
    X <- X0 %*% sdx.inv
  } else {
    X <- X0
  }
  
  
  if(is.null(MU)){
    MU <- FrechetMean(Y)
    
    mu0 <- colSums(Y) %>% {./norm(.,"2")}
    rbar <- colSums(Y) %>% {norm(.,"2")/n}
    kappa0 <- (rbar*q - rbar^3) / (1-rbar^2)
    
    MU <- mu0 * kappa0
  }
  
  
  
  # 1
  # {
  #   beta <- lm(Y ~ -1 + X, offset=tcrossprod(rep(1,n),MU)) %>% coef %>% as.vector
  #   beta <- c(MU, beta)
  #   # beta[1:3] <- beta[1:3] %>% {./norm(.,"2")}
  #   X <- cbind(1, X)
  #   MU <- 0
  # }
  
  
  # 2
  {
    beta <- lm(Y ~ -1 + X, offset=tcrossprod(rep(1,n),MU)+Offset) %>% coef %>% as.vector
    # beta <- lm(Y ~ 1 + X) %>% coef %>% {.[-1,]} %>% as.vector
  }
  
  
  
  l1 <- 1
  beta_old <- beta
  beta_new <- beta + 1
  beta.list <- loglik.list <- crit.list <- NULL
  crit <- 1
  while( crit > eps & l1 <= maxit ){
    
    beta.list[[l1]] <- beta
    
    beta_old <- beta
    
    Xt.list <- apply(X, 1, function(Xi) kronecker(diag(1,q), Xi), simplify=FALSE)
    Xt <- do.call("cbind", Xt.list)
    yt <- as.vector(t(Y))
    
    b1.list <- lapply(1:n, function(i) b1.vMF( MU + Offset[i,] + crossprod(Xt.list[[i]], beta) ))
    eta <- do.call("rbind", b1.list)
    
    b2.list <- lapply(1:n, function(i) b2.vMF( MU + Offset[i,] + crossprod(Xt.list[[i]], beta) ))
    W <- do.call("adiag", b2.list) # magic::adiag
    
    Z <- crossprod(Xt, beta) + solve(W) %*% (yt - eta)
    
    beta <- solve(Xt %*% W %*% t(Xt) + diag(lambda,p*q,p*q)) %*% Xt %*% W %*% Z
    
    # beta[1:3] <- beta[1:3] %>% {./norm(.,"2")}
    
    beta_new <- beta
    
    loglik <- Reduce("+", lapply(1:n, function(i){
      theta <- MU + Offset[i,] + crossprod(Xt.list[[i]], beta)
      crossprod(theta, Y[i,]) + log( Cq(theta) )
    }))
    
    crit <- norm(beta_old-beta_new, "2")
    
    loglik.list[l1] <- loglik
    crit.list[l1] <- crit
    
    l1 <- l1 + 1
  }
  
  # print( matrix(beta,p,q) )
  # print( l1 )
  
  
  # plot(loglik.list, type="l")
  
  beta.list <- lapply(beta.list, function(bb) matrix(bb, p, q, byrow=FALSE) )
  beta_new <- matrix(beta_new, p, q, byrow=FALSE)
  
  if(standardize){
    beta.list <- lapply(beta.list, function(b) sdx.inv %*% b)
    beta_new <- sdx.inv %*% beta_new
  }
  
  
  
  list(mu=MU, 
       beta=beta_new,
       beta.list=beta.list,
       offset=Offset,
       lambda=lambda,
       # beta=matrix(beta,p,q,byrow=TRUE),
       # beta.list=beta.list %>% lapply(function(x) matrix(x,p,q,byrow=TRUE)),
       # beta=matrix(beta,p+1,q,byrow=TRUE),
       # beta.list=beta.list %>% lapply(function(x) matrix(x,p+1,q,byrow=TRUE)),
       loglik.list=loglik.list[-1], crit.list=crit.list[-1])
  
}












#' @export glm_vmf_FixedMean_Offset2
glm_vmf_FixedMean_Offset2 <- function(X, Y, MU=NULL, Offset=NULL, lambda=1e-12, beta0=NULL, maxit=1000, eps=1e-8, standardize=TRUE){
  
  if(FALSE){
    X=simdata$X; Y=simdata$Y; MU=NULL; Offset=NULL; lambda=1e-2; maxit=100; eps=1e-8; standardize=TRUE
  }
  
  
  if(FALSE){
    ii <- pvec.list[[d]]
    X=U.joint[,-1]
    Y=X_joint[,ii]
    Offset=Theta.inds[,ii]
    
    beta0=NULL; maxit=1000; eps=1e-8; standardize=TRUE
  }
  
  
  
  library(magic)
  
  n <- nrow(Y)
  q <- ncol(Y)
  p <- ncol(X)
  
  if(is.null(Offset)){
    Offset <- matrix(0,n,q)
  }
  
  X0 <- scale(X, center=TRUE, scale=FALSE)
  
  if(standardize){
    sdx.inv <- apply(X0,2,sd) %>% {diag(1/., length(.), length(.))}
    X <- X0 %*% sdx.inv
  } else {
    X <- X0
  }
  
  
  if(is.null(MU)){
    MU <- FrechetMean(Y)
    
    muhat <- colSums(Y) %>% {./norm(.,"2")} # mle of mean direction in vMF
    rbar <- colSums(Y) %>% {norm(.,"2")/n}
    kappa0 <- (rbar*q - rbar^3) / (1-rbar^2) # approximation
    
    MU <- muhat * kappa0
  }
  
  
  {
    B0 <- lm(Y ~ -1 + X, offset=tcrossprod(rep(1,n),MU)+Offset) %>% coef
    beta <- as.vector(rbind(MU, B0))
  }
  
  
  X1 <- cbind(1, X)
  
  l1 <- 1
  beta_old <- beta
  beta_new <- beta + 1
  beta.list <- loglik.list <- crit.list <- NULL
  crit <- 1
  while( crit > eps & l1 <= maxit ){
    
    beta.list[[l1]] <- beta
    
    beta_old <- beta
    
    Xt.list <- apply(X1, 1, function(Xi) kronecker(diag(1,q), Xi), simplify=FALSE)
    Xt <- do.call("cbind", Xt.list)
    yt <- as.vector(t(Y))
    
    b1.list <- lapply(1:n, function(i) b1.vMF( Offset[i,] + crossprod(Xt.list[[i]], beta) ))
    eta <- do.call("rbind", b1.list)
    
    b2.list <- lapply(1:n, function(i) b2.vMF( Offset[i,] + crossprod(Xt.list[[i]], beta) ))
    W <- do.call("adiag", b2.list) # magic::adiag
    
    Z <- crossprod(Xt, beta) + solve(W) %*% (yt - eta)
    
    beta <- solve(Xt %*% W %*% t(Xt) + diag(lambda,q*(p+1),q*(p+1))) %*% Xt %*% W %*% Z
    
    # beta[1:3] <- beta[1:3] %>% {./norm(.,"2")}
    
    beta_new <- beta
    
    loglik <- Reduce("+", lapply(1:n, function(i){
      theta <- Offset[i,] + crossprod(Xt.list[[i]], beta)
      crossprod(theta, Y[i,]) + log( Cq(theta) )
    }))
    
    crit <- norm(beta_old-beta_new, "2")
    
    loglik.list[l1] <- loglik
    crit.list[l1] <- crit
    
    l1 <- l1 + 1
  }
  
  # print( matrix(beta,p,q) )
  # print( l1 )
  
  
  # plot(loglik.list, type="l")
  
  beta.list <- lapply(beta.list, function(bb) matrix(bb, p+1, q, byrow=FALSE) )
  beta_new <- matrix(beta_new, p+1, q, byrow=FALSE)
  
  if(standardize){
    beta.list <- lapply(beta.list, function(b) magic::adiag(1, sdx.inv) %*% b)
    beta_new <- magic::adiag(1, sdx.inv) %*% beta_new
  }
  
  
  
  list(mu=beta_new[1,], 
       beta=beta_new,
       beta.list=beta.list,
       offset=Offset,
       lambda=lambda,
       # beta=matrix(beta,p,q,byrow=TRUE),
       # beta.list=beta.list %>% lapply(function(x) matrix(x,p,q,byrow=TRUE)),
       # beta=matrix(beta,p+1,q,byrow=TRUE),
       # beta.list=beta.list %>% lapply(function(x) matrix(x,p+1,q,byrow=TRUE)),
       loglik.list=loglik.list[-1], crit.list=crit.list[-1])
  
}
