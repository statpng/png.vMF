


png.FrechetMean <- function(X, lr=1, maxit=1000, tol=1e-10, verbose=FALSE){
  
  if(FALSE){
    set.seed(1)
    X=sim.sphere(n=10); lr=1e-2; maxit=1000; tol=1e-10; verbose=T
    
    FrechetMean(X, tol=1e-10)
    png.FrechetMean(X, tol=1e-10)
  }
  
  y.init <- X[1,,drop=F]
  
  yold <- y.init
  for( i in 1:maxit ){
    tmp <- lr*rowMeans(apply(X, 1, function(x) LogmapSphere(mu=yold, x)))
    ynew <- Expmu(mu=yold, matrix(tmp,nrow=1))
    
    # crit <- GeodRegr::geo_dist("sphere", yold, ynew)
    crit <- acos( sum(yold*ynew) )
    # crit <- norm( Logmu(yold,ynew), "2" )
    
    
    if(verbose) print(crit)
    if(crit < tol) break
    
    yold <- ynew
  }
  
  matrix(ynew, ncol=1)
}

FrechetMean <- function(X, ...){
  manifold::frechetMean(manifold::createM("Sphere"), t(X), ...)
}



#' @export vMF.MuKappa
vMF.MuKappa <- function(Y, kappa.type=4){
  n <- nrow(Y)
  q <- ncol(Y)
  
  muhat <- colSums(Y) %>% {./norm(.,"2")} # mle of mean direction in vMF
  rbar <- colSums(Y) %>% {norm(.,"2")/n}
  
  if(is.null(kappa.type)){
    kappa0 <- 1
  } else {
    kappa0 <- list(
      (q-1) / 2*(1-rbar),
      q*rbar*( 1+q/(q+2)*rbar^2 + q^2*(q+8)/((q+2)^2*(q+4))*rbar^4 ),
      (rbar*q) / (1-rbar^2),
      (rbar*q - rbar^3) / (1-rbar^2)
    )[[kappa.type]]
  }
  
  muhat * kappa0
  
}


#


rotation.x2y <- function(x, y) {
  if(FALSE){
    x <- rnorm(5) %>% {./norm(.,"2")}
    y <- rnorm(5) %>% {./norm(.,"2")}
    rotation.x2y(x, y) %*% x
    y
  }
  # 반사시키는 벡터 v 계산
  v <- x - y
  # 벡터 v를 정규화
  v <- v / sqrt(sum(v * v))
  # 단위 행렬 생성
  I <- diag(length(x))
  # 하우스홀더 변환을 사용하여 회전 행렬 Q 계산
  Q <- I - 2 * (v %*% t(v))
  return(Q)
}


some.useful.functions <- function(x,y){
  library(rotations)
  
  rotations::project.SO3()
  
  
  #
  rotations::rot.dist()
  #
  r <- rotations::rfisher(10000, kappa = 0.01)
  #
  r <- rotations::rvmises(10000, kappa = 0.01)
  range(r)
  
  #
  rotations::genR(10)
  #
  
  
  library(rotasym)
  
}





#' @export LogmapSphere
LogmapSphere <- function(mu, v) {
  # 두 점 사이의 각도 계산
  theta <- acos(sum(mu * v))
  
  # 0으로 나누는 것을 방지
  if (theta < .Machine$double.eps) {
    return(rep(0, length(mu)))
  }
  
  # 구면 로그 맵 계산
  return(theta / sin(theta) * (v - cos(theta) * mu))
  
  # apply(data$X[[2]] %>% head, 1, function(x){
  #   print(log_map('sphere', c(0, 0, 1), x))
  #   print(LogmapSphere(c(0, 0, 1), x))
  # })
  # 
  # log_map('sphere', c(0, 0, 1), c(0, 1/sqrt(2), 1/sqrt(2)))
  # LogmapSphere(c(0, 0, 1), c(0, 1/sqrt(2), 1/sqrt(2)))
  
}



#' @export Expmu
Expmu <- function(mu, V) {
  n <- dim(V)[1]
  # lv <- mapply(function(x1, x2, x3) norm(c(x1, x2, x3), type="2"), V[, 1], V[, 2], V[, 3])
  lv <- as.vector( apply(V, 1, norm, "2") )
  cos_lv <- diag(cos(lv), n, n)
  sin_lv <- diag(sin(lv) / lv, n, n)
  im <- cos_lv %*% kronecker(matrix(mu,nrow=1), matrix(1,n,1)) + sin_lv %*% V
  return(im)
}


#' @export Logmu
Logmu <- function(mu, X){
  # library(GeodRegr)
  # t(apply(X, 1, function(x){
  #   GeodRegr::log_map("sphere", mu, x)
  # }))
  
  t(apply(X, 1, function(x){
    LogmapSphere(mu, x)
  }))
  
}


#' @export sphere.Frechet
sphere.Frechet <- function(df){
  
  # df: n x 3 data.frame (x,y,z)
  
  lapply(1:nrow(df), function(i) as.vector(df[i,])) %>% 
    riemfactory(., name="sphere") %>% 
    rbase.mean(.)
  
}





if(FALSE){
  library(GeodRegr)
  
}



if(FALSE){
  library(RiemBase)
  
  ndata = 100
  theta = seq(from=-0.99,to=0.99,length.out=ndata)*pi
  tmpx  = cos(theta) + rnorm(ndata,sd=0.1)
  tmpy  = sin(theta) + rnorm(ndata,sd=0.1)
  
  ### Wrap it as 'riemdata' class
  data  = list()
  for (i in 1:ndata){
    tgt = c(tmpx[i],tmpy[i],1)
    data[[i]] = tgt/sqrt(sum(tgt^2)) # project onto Sphere
  }
  data = riemfactory(data, name="sphere")
  
  ### Compute Fréchet Mean
  out1 = rbase.mean(data)
  out2 = rbase.mean(data,parallel=TRUE) # test parallel implementation
  
  
}





if(FALSE){
  library(manifold)
  detach("package:manifold", unload=TRUE)
  
  manifold::runifSphere(100, 3) %>% t %>% png.sphere
  
  
  manifold::frechetMean(createM("Sphere"), manifold::runifSphere(100, 3))
  
}

