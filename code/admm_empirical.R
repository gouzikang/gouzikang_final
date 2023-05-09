library(readr)
#library(abind)
library(glmnet)
library(methods)
library(rqPen)
library(doParallel)
library(RSpectra)
library(mvtnorm)
library(MASS)






#m <- sum(ryn<0)#累计收益率为负值的股票数

# #单只股票的夏普率
# rf <- rep(0,p)
# sd <- rep(0,p)
# sr1 <- rep(0,p)
# 
# #计算方式1
# for (i in 1:p){
#   rf[i] <- mean(R[,i])#日回报率的均值
#   sd[i] <- sd(R[,i])#日回报率的标准差
#   sr1[i] <- rf[i]/sd[i]*sqrt(252)
# }
# 
# #计算方式2
# sr2 <- rep(0,p)
# for (i in 1:p){
#   sr2[i] <- ryn[i]*100/sd(R[,i])/sqrt(252)
# }


p=477
# zhongshu <- function(x){
#   return(as.numeric(names(table(x))[table(x)==max(table(x))]))
# }

l2 <- function(x){
  sum(x^2)
}
shrink <- function(x,y){
  sign(x)*pmax(0, (abs(x)-y))
}

lpd_admm_constrained <- function(R, lambda,
                                 max_iter, eabs, erel,
                                 theta, rho, r){
  
  Sigma <- cov(R)
  DiffMean <- apply(R,2,mean)
  DiffMean <- matrix(DiffMean,p,1)
  
  omega <- rep(1/p,p) 
  z <- rep(0,p)
  gamma <- rep(0,p)
  tao <- rho*eigs_sym(t(Sigma)%*%Sigma, 1)$values
  phi <- 0
  Alpha <- rep(1,p)
  
  
  Sigmap <- rbind(Sigma,Alpha)
  DiffMeanp <- rbind(DiffMean,1/r)
  
  Z <- rbind(matrix(z,p,1),0)
  Z.1 <- Z
  Gamma <- rbind(matrix(gamma,p,1),phi)
  #omega.m<-matrix(0, nrow=p/K, ncol=K)
  #  d.3 <- l2(as.vector(beta.m)-beta.0)
  dualres = sqrt(l2(rho*t(Sigmap)%*%(Z-Z.1)))
  primres = sqrt(l2(Sigmap%*%omega-Z-r*DiffMeanp))
  #dualres <- sqrt(l2(phi*t(X)%*%(z.1-z)))
  epri = sqrt(p+1)*eabs + erel*max(sqrt(l2(Sigmap%*%omega)), sqrt(l2(Z)), sqrt(l2(r*DiffMeanp)))
  edual = sqrt(p)*eabs + erel*sqrt(l2(t(Sigmap)%*%Gamma))
  
  
  k <- 1
  while ((primres > epri | dualres > edual ) & k < max_iter){
    
    z.1 <- z
    Z.1 <- rbind(matrix(z.1,p,1),0)
    
    omega <- shrink(omega-rho/tao*(t(Sigma)%*%(Sigma%*%omega-r*DiffMean-z+1/rho*gamma)+phi/rho*Alpha),1/tao)
    
    z <- pmin(pmax(Sigma%*%omega-r*DiffMean + 1/rho*gamma,-lambda*Alpha), lambda*Alpha)
    Z <- rbind(matrix(z,p,1),0)
    
    gamma <- gamma+theta*rho*(Sigma%*%omega-r*DiffMean-z)
    
    phi <- phi +theta*rho*(sum(omega)-1)
    
    Gamma <- rbind(matrix(gamma,p,1),phi)
    
    dualres = sqrt(l2(rho*t(Sigmap)%*%(Z-Z.1)))
    primres = sqrt(l2(Sigmap%*%omega-Z-r*DiffMeanp))
    epri = sqrt(p+1)*eabs + erel*max(sqrt(l2(Sigmap%*%omega)), sqrt(l2(Z)), sqrt(l2(r*DiffMeanp)))
    edual = sqrt(p)*eabs + erel*sqrt(l2(t(Sigmap)%*%Gamma))
    
    k <- k+1 
    
    #   print(c(k, l2(as.vector(beta.m)-beta.0), sum(beta.m!=0)))
  }
  print(k)
  return(omega)
}





R_in <- read.csv("/Users/gouzikang/Desktop/博士/博二/2023春/高级应用统计/final-pre/insample.csv")
R_out <- read.csv("/Users/gouzikang/Desktop/博士/博二/2023春/高级应用统计/final-pre/outofsample.csv")
R_in <- as.matrix(R_in)*100
R_out <- as.matrix(R_out)*100


  mu_in <- apply(R_in,2,mean)
  mu_in <- as.matrix(mu_in)
  Sigma_in <- cov(R_in)
  
  mu_out <- apply(R_out,2,mean)
  mu_out <- as.matrix(mu_out)
  Sigma_out <- cov(R_out)
  
  # DiffMean_out <- apply(R_out,2,mean)
  # DiffMean_out <- matrix(DiffMean_out,p,1)
  # Sigma_out <- cov(R_out)
  
  # # ##in_sample
  # #单只股票的累计收益率
  # rtn <- rep(1,p)
  # ryn <- rep(0,p)
  # for (i in (1:p)){
  #   for (j in (1:n1)){
  #     rtn[i] <- rtn[i]*(1+R[j,i]/100)
  #   }
  #   rtn[i] <- rtn[i]-1#累计收益率
  # }
  # 
  # T <- n1/252
  # ryn <- (rtn+1)**(1/T)-1#年累计收益率
  
  
  res_admm_constrained=lpd_admm_constrained(R_in, lambda=0.1,
                                            max_iter=10000,
                                            eabs=0.001,
                                            erel=0.001,
                                            theta=1.618,
                                            rho=0.25,
                                            r=1/2)
  omega <- matrix(res_admm_constrained,p,1)
  return_in <- t(omega)%*%mu_in
  risk_in <- t(omega)%*%Sigma_in%*%omega
  #sharpe_ratio1 <- Ry1/Sd/sqrt(252)
  sharpe_ratio_in <- return_in/sqrt(risk_in)*sqrt(252)
  target_in <- 1/2*risk_in-1/2*return_in
  
  
  max<- max(omega)
  min<- min(omega)
  long <- sum(omega>0)
  short <- sum(omega<0)
  
  #out_of_sample
  return_out <- t(omega)%*%mu_out
  risk_out<- t(omega)%*%Sigma_out%*%(omega)
  #sharpe_ratio1 <- Ry1/Sd/sqrt(252)
  sharpe_ratio_out <- return_out/sqrt(risk_out)*sqrt(252)
  target_out <- 1/2*risk_out-1/2*return_out



save.image(file = "ADMM.Rdata")

