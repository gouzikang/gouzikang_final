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

lpd_qp_constrained <- function (R,r,lambda,rho,theta,max_iter){
  
  Sigma <- cov(R)
  DiffMean <- apply(R,2,mean)
  DiffMean <- matrix(DiffMean,p,1)
  
  omega <- rep(1/p,p) 
  eta <- eigs_sym(t(Sigma)%*%Sigma, 1)$values
  phi <- 0
  Alpha <- rep(1,p)
  k <- 1
  omega.1 <- omega
  epri <- sum(abs(omega.1-omega))+0.1
  epsilon <- 1e-4
  #===  update === #
  
  while (epri > epsilon & k < max_iter){
    
    omega.1 <- omega
    
    omega <- shrink(omega-1/eta*Sigma%*%omega-phi/eta*Alpha+r/eta*DiffMean,lambda/eta)
    
    phi <- phi +theta*rho*(sum(omega)-1)
    
    epri <- sum(abs(omega.1-omega))
    
    k <- k+1 
  }
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


res_qp_constrained=lpd_qp_constrained(R_in,
                                      r=1/2,
                                      lambda=0.1,
                                      rho=0.25,
                                      theta=1.618,
                                      max_iter=10000)

omega <- matrix(res_qp_constrained,p,1)
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



save.image(file = "lpsolve.Rdata")

