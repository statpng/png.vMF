function(){

  library(dplyr)  
  devtools::load_all()
  
  simdata <- sim.sphere(n=50, p=2, r=1, mu=c(0,0,1), snr=10, s=2)
  
  simdata$X %>% png.sphere()
  
  
  simdata2 <- sim.glm_vmf(n=50, p=1, q=3, sd=50, mu=c(0,0,50), orthogonal=T)
  
  simdata2$Y %>% png.sphere()
  
  
  fit <- glm_vmf_FixedMean_Offset(X=simdata2$X, Y=simdata2$Y, lambda=0)
  
  fit$beta
  
}







function(){
  
  devtools::load_all()
  
  
  set.seed(1)
  out.beta <- array(NA, c(3, 4, 5, 100))
  out <- array(NA, c(1, 3, 5, 100))
  for( j in 1:dim(out)[3] ){
    print(j)
    
    for( i in 1:dim(out)[4] ){
      
      tau <- c(5,10,20,50,100)[j]
      simdata <- sim.glm_vmf(n=50, p=1, q=3, sd=100, mu=c(0,0,100), orthogonal=TRUE)  # png.sphere(simdata$Y)
      
      fit.GeodRegr <- with(simdata, GeodRegr::geo_reg("sphere", X, t(Y), estimator="l2", max_iter=1e+2))
      fit.vmf <- with(simdata, glm_vmf_FixedMean_Offset(X=X, Y=Y, MU=c(0,0,tau), lambda=1e-2, maxit=100, eps=1e-8, standardize=TRUE))
      
      fit.vmf <- glm_vmf_FixedMean_Offset2(X=simdata$X, Y=simdata$Y, lambda=1e-12, maxit=100, eps=1e-16, standardize=TRUE)
      
      
      fit.vmf$beta
      simdata$mu
      simdata$B
      
      # install.packages("Directional")
      # library(Directional)
      # y <- rvmf(150, rnorm(3), 5)
      # b1 <- vmf.reg(y, iris[, 4])
      # b2 <- glm_vmf_FixedMean(iris[,4,drop=F], y)
      # b1$beta
      # rbind(b2$mu, b2$beta)
      # GeodRegr::geo_reg(manifold="sphere", x=as.matrix(iris[,4,drop=F]), y=t(y), estimator="l2") %>% {
      #   t(cbind(.$p, .$V))
      # }
      
      
      
      result <- list(simdata$B, 
                     t(fit.GeodRegr$V), 
                     fit.vmf$beta) %>% 
        lapply(t) %>% lapply(function(x) x/norm(x,"F"))
      
      result
      
      
      out[,1,j,i] <- png.utils::png.angle(result[[1]], result[[2]])$max
      out[,2,j,i] <- png.utils::png.angle(result[[1]], result[[3]])$max
      # out[,3,j,i] <- png.utils::png.angle(result[[1]], result[[4]])$max
      
      out.beta[,1,j,i] <- simdata$B
      out.beta[,2,j,i] <- t(fit.GeodRegr$V)
      out.beta[,3,j,i] <- fit1$beta
      # out.beta[,4,j,i] <- fit2$beta[-1,]
    }
  }
  
  if(FALSE){
    save(out, out.beta, file="out.glm_vmf1.RData")
  }
  
  load("out.glm_vmf1.RData")
  
  
  apply(out,2:3,mean)
  #          [,1]     [,2]     [,3]     [,4]     [,5]
  # [1,] 10.68148 5.796712 3.926280 4.011802 4.258179
  # [2,] 12.87231 7.023528 4.696356 3.490202 2.928996
  
  
  
  for( idx.sd in 1:dim(out.beta)[3] ){
    pdf(file=paste0("glm_vmf1-ScatterPlot", idx.sd, ".pdf"), width=5, height=5)
    
    par(mfrow=c(2,2), omi=c(.2,.4,.2,.2), mai=c(.8,.6,.2,.2))
    cbind.data.frame(true=out.beta[1,1,idx.sd,], est=out.beta[1,2,idx.sd,]) %>% {plot(., pch=18); cor(.) %>% {.[2]} %>% round(4) %>% mtext(., side=3, line=-1, adj=.01)}
    mtext("GeodRegr", 3, line=.5)
    mtext("Y1", 2, line=4.5, outer=F)
    cbind.data.frame(true=out.beta[1,1,idx.sd,], est=out.beta[1,3,idx.sd,]) %>% {plot(., pch=18); cor(.) %>% {.[2]} %>% round(4) %>% mtext(., side=3, line=-1, adj=.01)}
    mtext("IRLS with vMF", 3, line=.5)
    cbind.data.frame(true=out.beta[1,1,idx.sd,], est=out.beta[1,2,idx.sd,]) %>% {plot(., pch=18); cor(.) %>% {.[2]} %>% round(4) %>% mtext(., side=3, line=-1, adj=.01)}
    mtext("Y2", 2, line=4.5, outer=F)
    cbind.data.frame(true=out.beta[1,1,idx.sd,], est=out.beta[1,3,idx.sd,]) %>% {plot(., pch=18); cor(.) %>% {.[2]} %>% round(4) %>% mtext(., side=3, line=-1, adj=.01)}
    
    dev.off()
  }
  
  
  
  
  
  
}