library(readr)
#library(abind)
library(glmnet)
library(methods)
library(rqPen)
library(doParallel)
library(RSpectra)
library(mvtnorm)
library(MASS)

n=756
n1=504
n2=n-n1
p=1000

mub <- c(0.78282,0.51803,0.41003)
covb <- matrix(c(0.029145,0.023873,0.010184,0.023873,0.053951,-0.006967,0.010184,-0.006967,0.086856),3,3)
muf <- c(0.023558,0.012989,0.020714)
covf <- matrix(c(1.2507,-0.034999,-0.20419,-0.034999,0.31564,-0.0022526,-0.20419,-0.0022526,0.19303),3,3)


set.seed(2020)
B <- mvrnorm(p, mub, covb) #factor loading
std <- rgamma(p,shape=3.3586,scale=0.1876)
mue <- rep(0,p)
cove <- diag(std**2)




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



T=5
return_in <- rep(0,T)
risk_in <- rep(0,T)
target_in <- rep(0,T)
sharpe_ratio_in <- rep(0,T)

OMEGA <- matrix(0,p,T)


return_out <- rep(0,T)
risk_out <- rep(0,T)
target_out <- rep(0,T)
sharpe_ratio_out <- rep(0,T)

long <- rep(0,T)
short <- rep(0,T)
max<- rep(0,T)
min<- rep(0,T)

for (t in 1:T){
  set.seed(2019+t)
  F <- mvrnorm(n,muf,covf)
  F <-t(F)
  
  e <- mvrnorm(n,mue,cove)
  e <- t(e)
  R <- B%*%F+e
  
  R_in <- R[,1:n1] #in_sample
  R_out<- R[,(n1+1):n] #out_of_sample
  # 
  R_in <- t(R_in)
  R_out <- t(R_out)
  
  mu_in <- apply(R_in,2,mean)
  mu_in <- matrix(mu_in,p,1)
  Sigma_in <- cov(R_in)
  
  mu_out <- apply(R_out,2,mean)
  mu_out <- matrix(mu_out,p,1)
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
                                        lambda=0.04,
                                        rho=0.25,
                                        theta=1.618,
                                        max_iter=10000)
  
  omega <- matrix(res_qp_constrained,p,1)
  return_in[t] <- t(omega)%*%mu_in
  risk_in[t] <- t(omega)%*%Sigma_in%*%omega
  #sharpe_ratio1 <- Ry1/Sd/sqrt(252)
  sharpe_ratio_in[t] <- return_in[t]/sqrt(risk_in[t])*sqrt(252)
  target_in[t] <- 1/2*risk_in[t]-1/2*return_in[t]
  OMEGA[,t] <- omega
  
  max[t]<- max(OMEGA[,t])
  min[t]<- min(OMEGA[,t])
  long[t] <- sum(OMEGA[,t]>0)
  short[t] <- sum(OMEGA[,t]<0)
  
  #out_of_sample
  return_out[t] <- t(OMEGA[,t])%*%mu_out
  risk_out[t] <- t(OMEGA[,t])%*%Sigma_out%*%(OMEGA[,t])
  #sharpe_ratio1 <- Ry1/Sd/sqrt(252)
  sharpe_ratio_out[t] <- return_out[t]/sqrt(risk_out[t])*sqrt(252)
  target_out[t] <- 1/2*risk_out[t]-1/2*return_out[t]
}


save.image(file = "qp.Rdata")

