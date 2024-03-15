#' @export sim.sphere
sim.sphere <- function(n=50, p=2, r=1, mu=c(0,0,1), snr=10, s=.2, s0=0){
  # typeA = "sparse", typeB = "all", n = 50, p = 10, q = 10, 
  #  d = 3, rvec = NULL, nuA = 0.2, nuB = 0.5, d0 = 3, es = "1", 
  #  es.B = 1, snr = 1, simplify = TRUE, sigma = NULL, rho_X = 0.5, 
  #  rho_E = 0
  
  if(FALSE){
    with(as.data.frame(sim.sphere(r=1, snr=5, s=.5)), png.sphere(V1,V2,V3))
  }
  if(FALSE){
    n=10; p=2; r=2; mu=c(0,0,1); snr=5; s=.5; s0=0; zero.prop=0.5; seed=1;
    seed.U=seed; seed.V=seed; verbose=FALSE
  }
  
  # mu <- rep(0,p)
  
  
  V.tmp <- do.call("cbind", lapply(1:r, function(x) rnorm(p)))
  V <- qr.Q(qr( V.tmp ))
  # D <- sapply(1:r, function(k) (s/k + s0))
  D <- sapply(1:r, function(k) (s/1 + s0))
  
  
  U.tmp <- do.call("cbind", lapply(1:r, function(x) rnorm(n)))
  U <- qr.Q(qr( U.tmp )) %*% diag(D,r,r)
  # U <- matrix(0,n,r)
  # for( k in 1:r ){
  #   U[,k] <- rnorm(n, mean=0, sd=D[k])
  # }
  
  X0 <- tcrossprod(U, V)
  # X0 <- tcrossprod(rep(1,n), mu) + tcrossprod(U,V)
  
  E <- matrix(rnorm(n*p,0,1),n,p)
  sigma <- sqrt(sum(as.numeric(X0)^2)/sum(as.numeric(E)^2)/snr)
  E <- E * sigma
  
  X <- Expmu(mu, cbind( scale(X0 + E, scale=F), rep(0,n)))
  
  
  result <- list(mu=mu, U=U, V=V, X=X, X0=X0, E=E)
  
  result
}




#' @export sim.sphere.runif
sim.sphere.runif <- function(n){
  aa <- seq(0, 360, length.out = n)
  bb <- seq(0, 180, length.out = n)
  grid <- expand.grid(theta = aa, phi = bb)
  
  x <- sin(grid$theta * pi / 180) * cos(grid$phi * pi / 180)
  y <- sin(grid$theta * pi / 180) * sin(grid$phi * pi / 180)
  z <- cos(grid$theta * pi / 180)
  
  list(x=x,y=y,z=z)
}


#' @export sim.sphere.circle
sim.sphere.circle <- function(n){
  
  ndata = n
  theta = seq(from=-0.99,to=0.99,length.out=ndata)*pi
  tmpx  = cos(theta) + rnorm(ndata,sd=0.1)
  tmpy  = sin(theta) + rnorm(ndata,sd=0.1)
  
  data <- cbind.data.frame(x=tmpx, y=tmpy, z=1)
  data <- t( apply(data, 1, function(x) x/norm(x,"2")) )
  
  # data = riemfactory(data, name="sphere")
  
  return(data)
}















#' @export sim.glm_vmf
sim.glm_vmf <- function( n=100, p=1, q=3, sd=1, mu=c(0,0,1), orthogonal=TRUE ){
  
  if(FALSE){
    n=50; p=1; q=3; sd=1; mu=c(0,0,1); orthogonal=TRUE
  }
  
  X <- matrix(rnorm(n*p, sd=sd), n, p)
  if(orthogonal){
    B <- matrix(runif(p*q, -1, 1), p, q) %>% 
      # svd() %>% { diag(.$d, length(.$d), length(.$d)) %*% t(.$v) }
      apply(1,function(v){
      qr.Q(qr( cbind(mu, v) ))[,-1]
    }) %>% t()
  } else {
    B <- matrix(runif(p*q, -1, 1), p, q)
  }
  
  
  
  theta <- tcrossprod(rep(1,n), mu) + X %*% B
  
  Y <- t(apply(theta, 1, function(theta.i) Rfast::rvmf(1, theta.i/norm(theta.i,"2"), k=norm(theta.i,"2")) ))
  
  list(X=X, Y=Y, mu=mu, B=B, theta=theta)  
  
}








