require(msm)

r<- 1E2

N <- 1E3

X0 <- 1#rnorm(N, 1, 1)
X1 <- rnorm(N, 1, 1)
X <- cbind(X0, X1)
beta <- c(-3, 0)


theta<- pnorm(X%*%beta)

# y<- as.numeric( runif(N)<theta)
# sum(y)
y<- rep(0,N)
y[1]<-1

vecInf<- rep(Inf, N)

X2inv<- solve(t(X)%*%X)
cholX2inv<- t(chol(X2inv))
lb<- -vecInf
ub<- vecInf
lb[y==1]<- 0
ub[y==0]<- 0

# beta<- rnorm(2)

# r<- 1E5

trace_beta<- numeric()
trace_r2<- numeric()

for(i in 1:300){
  
  theta<- X%*%beta *r
  Z<- rtnorm(N, theta, r, lb, ub )
  m<- X2inv%*% (t(X)%*% (Z)) /r
  
  
  
  beta <- cholX2inv%*%rnorm(2) + m 

  if(i>100){
    trace_beta<- rbind(trace_beta, t(beta))
    trace_r2<- c(trace_r2, r2)
  }
  print(i)
}

beta
 
ts.plot(trace_beta[,1])
acf(trace_beta[,1],lag.max = 40)
# acf(trace_beta[,2],lag.max = 40)

# ts.plot(trace_r2)
# ts.plot(trace_beta[,2])
# acf(trace_beta[,1],lag.max = 40)
