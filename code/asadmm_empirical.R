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

lpd_asadmm_constrained <- function(R, K,lambda,
                                   max_iter, eabs, erel,
                                   theta, rho, r){
  
  
  Sigma <- cov(R)
  DiffMean <- apply(R,2,mean)
  DiffMean <- matrix(DiffMean,p,1)
  
  
  Omega<-matrix(1/p, nrow=p/K, ncol=K) #给omega分block
  z <- rep(0,p)
  phi<- 0
  alpha <- rep(1,p/K)
  Alpha <- rep(1,p)
  
  for (i in 1:K){
    assign(paste("Sigma",i,sep=""), Sigma[,seq((1+p/K*(i-1)),p/K*i)]) #给Sigma分块
  }
  
  abeta <- matrix(NA, nrow=p, ncol=K)
  for (i in 1:K){
    abeta[,i]<-get(paste("Sigma",i,sep=""))%*%Omega[,i]
  }#储存Sigma_i✖️omega_i
  
  
  beta <- matrix(0, nrow=p, ncol=K-1)
  for (i in 1:(K-1)){
    beta[,i] <- 1/K*(r*DiffMean+z+K*abeta[,i+1]-rowSums(abeta))
  }
  
  gamma <- matrix(0, ncol = K, nrow = p)
  
  
  OmegaLasso <- matrix(NA, nrow=p/K, ncol=K) 
  OmegaLasso[,1] <- t(get(paste("Sigma",1,sep="")))%*%(abeta[,1]+rowSums(beta) - z - r*DiffMean + gamma[,1]/rho) + phi/rho*alpha
  registerDoParallel(10)  # use multicore, set to the number of our cores
  OmegaLasso[,2:K] = foreach (i=2:K, .combine=cbind) %dopar% {
    t(get(paste("Sigma",i,sep="")))%*%(abeta[,i]-beta[,(i-1)]+gamma[,i]/rho) + phi/rho*alpha
  }
  
  eta = rep(NA, K)
  registerDoParallel(10)  # use multicore, set to the number of our cores
  eta = foreach (i=1:K, .combine=c) %dopar% {
    rho*eigs_sym(t(get(paste("Sigma",i,sep="")))%*%get(paste("Sigma",i,sep="")), 1)$values
  }
  
  
  z.1 <- z
  
  Xz = matrix(0, nrow = p/K, ncol = K)
  Xg = matrix(0, nrow = p/K, ncol = K)
  
  for (i in 1:K){
    Xz[,i]<-t(get(paste("Sigma",i,sep="")))%*%((z.1-z))
  }
  dualres = sqrt(l2(rho*c(Xz)))
  primres = sqrt(l2(abeta[,1]-z+rowSums(beta) - r*DiffMean)) + sqrt(l2(as.vector(beta-abeta[, 2:K]))) +abs(sum(Omega)-1)
  
  epri = sqrt(p)*eabs + erel*max(sqrt(l2(rowSums(abeta))), sqrt(l2(z)), sqrt(l2(r*DiffMean)))
  for (i in 1:K){
    Xg[,i]<-t(get(paste("Sigma",i,sep="")))%*%(gamma[,i])
  }
  edual = sqrt(p)*eabs + erel*sqrt(l2(c(Xg)))
  
  k <- 1
  
  #===  update === #
  
  while ((primres > epri | dualres > edual ) & k < max_iter){
    
    z.1 <- z
    
    registerDoParallel(10)  # use multicore, set to the number of our cores
    Omega[,1:K] = foreach (i=1:K, .combine=cbind) %dopar% {
      shrink(Omega[,i]-rho/eta[i]*OmegaLasso[,i], 1/eta[i])
    }
    
    for (i in 1:K){
      abeta[,i]<-get(paste("Sigma",i,sep=""))%*%Omega[,i]
    }
    
    for (i in 1:(K-1)){
      beta[,i] <- 1/K*(r*DiffMean+z+K*abeta[,i+1]-rowSums(abeta))
    }
    
    z <- pmin(pmax(abeta[,1]+rowSums(beta)-r*DiffMean + gamma[,1]/rho,-lambda*Alpha), lambda*Alpha)
    
    for (i in 1:(K-1)){
      beta[,i] <- 1/K*(r*DiffMean + z + K*abeta[,i+1]-rowSums(abeta))
    }
    
    gamma[,1] <- gamma[,1] + theta*rho*(abeta[,1]+rowSums(beta)-z-r*DiffMean)
    
    for (i in 2:K){
      gamma[,i] <- gamma[,i] + theta*rho*(abeta[,i]-beta[,i-1])
    }
    
    phi <- phi + theta*rho*(sum(Omega)-1)
    
    OmegaLasso[,1] <- t(get(paste("Sigma",1,sep="")))%*%(abeta[,1]+rowSums(beta) - z - r*DiffMean + gamma[,1]/rho) + phi/rho*alpha
    registerDoParallel(10)  # use multicore, set to the number of our cores
    OmegaLasso[,2:K] = foreach (i=2:K, .combine=cbind) %dopar% {
      t(get(paste("Sigma",i,sep="")))%*%(abeta[,i]-beta[,(i-1)]+gamma[,i]/rho) + phi/rho*alpha
    }
    
    
    for (i in 1:K){
      Xz[,i]<-t(get(paste("Sigma",i,sep="")))%*%((z.1-z))
    }
    dualres = sqrt(l2(rho*c(Xz)))
    primres = sqrt(l2(abeta[,1]-z+rowSums(beta) - r*DiffMean)) + sqrt(l2(as.vector(beta-abeta[, 2:K]))) +abs(sum(Omega)-1)
    
    epri = sqrt(p)*eabs + erel*max(sqrt(l2(rowSums(abeta))), sqrt(l2(z)), sqrt(l2(r*DiffMean)))
    for (i in 1:K){
      Xg[,i]<-t(get(paste("Sigma",i,sep="")))%*%(gamma[,i])
    }
    edual = sqrt(p)*eabs + erel*sqrt(l2(c(Xg)))
    
    # if (sqrt(l2(primres)) > 10*sqrt(l2(dualres))){
    #    phi = phi*2
    #  }
    #  if (sqrt(l2(dualres)) > 10*sqrt(l2(primres))){
    #    phi = phi/2
    #  }
    k <- k+1 
    
    
    #   print(c(k, l2(as.vector(beta.m)-beta.0), sum(beta.m!=0)))
  }
  print(k)
  
  return(Omega)
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


res_asadmm_constrained=lpd_asadmm_constrained(R_in, K=9,lambda=0.1,
                                              max_iter=10000, eabs=0.001, erel=0.001,
                                              theta=1.618, rho=1.1, r=0.5)
omega <- matrix(res_asadmm_constrained,p,1)
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



save.image(file = "ASADMM.Rdata")

