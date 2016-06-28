require(msm)

N <- 1E3

X0 <- 1#rnorm(N, 1, 1)
X1 <- rnorm(N, 1, 1)
X <- cbind(X0, X1)
beta <- c(0, -3)

theta<- pnorm(X%*%beta)

y<- as.numeric( runif(N)<theta)
sum(y)

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
trace_r2<- numeric()

r<- 3
for(i in 1:300){
  
  Xbeta  <- X%*%beta
  theta<- Xbeta
  Z<- rtnorm(N, theta, 1, lb, ub )
  m<- X2inv%*% (t(X)%*% (Z))
  
  # r2<- 1 / rgamma(1,N/2, rate = sum((Z*r-X%*%beta*r)^2)/2)
  
  # r<- sqrt(r2)
  
  
  beta <- cholX2inv%*%rnorm(2) + m /r

  if(i>100){
    trace_beta<- rbind(trace_beta, t(beta))
    trace_r2<- c(trace_r2, r2)
  }
  print(i)
}

beta
 
ts.plot(trace_beta[,2])
acf(trace_beta[,1],lag.max = 40)
# acf(trace_beta[,2],lag.max = 40)

# ts.plot(trace_r2)
# ts.plot(trace_beta[,2])
# acf(trace_beta[,1],lag.max = 40)
