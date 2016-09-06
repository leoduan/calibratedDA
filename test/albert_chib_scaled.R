require(msm)
rm(list=ls())
N <- 1E2

X0 <- 1
X1 <- rnorm(N, 1, 1)
X <- cbind(X0, X1)
beta <- c(0, 4)


theta<- pnorm(X%*%beta)

y<- runif(N)< theta

vecInf<- rep(Inf, N)

X2inv<- solve(t(X)%*%X)
cholX2inv<- t(chol(X2inv))
lb<- -vecInf
ub<- vecInf
lb[y==1]<- 0
ub[y==0]<- 0

 beta<- rnorm(2)

# r<- 1E5

trace_beta<- numeric()
trace_r<- numeric()

# r<- 1.01

r<-1


for(i in 1:500){
  
  xbeta<- X%*%beta
  Z<- rtnorm(N, xbeta, 1, lb, ub )
  
  beta_hat<- X2inv%*% (t(X)%*% (Z))
  
  RSS = sum((Z- X%*%beta_hat)^2)
  
  # r<- sqrt(RSS/ rchisq(1, N))

  
  m<-  beta_hat/r
  
  beta <- cholX2inv%*%rnorm(2) + m 

  if(i>100){
    trace_beta<- rbind(trace_beta, t(beta))
    trace_r<- c(trace_r, r)
  }
  print(i)
}

beta
 
ts.plot(trace_beta[,2])
acf(trace_beta[,2],lag.max = 40)
# acf(trace_beta[,2],lag.max = 40)

# ts.plot(trace_r)
# ts.plot(trace_beta[,2])
# acf(trace_beta[,1],lag.max = 40)


