
n <- 1E4

X0 <- 1#rnorm(N, 1, 1)
# X1 <- rnorm(n, 1, 1)
# X <- cbind(X0, X1)
# beta <- c(-1, 1)
X<- matrix(rep(1,n))

beta<- -3
p<- length(beta)
# hist(X %*% beta)
Z <- rnorm(n,X %*% beta,1)


y<- (Z>0)*1

table(y)

probit_r<- function(r,steps=1000, ver =1, back_iter=50){
  
  trace_beta<- numeric()
  trace_g<- numeric()
  
  
  lb<- rep(-Inf,n)
  ub<- rep(Inf,n)
  lb[y==1]<- 0
  ub[y==0]<- 0
  
  for(i in 1:steps){
    
    
    if(i>back_iter){
      beta<- trace_beta[i- ceiling(runif(1,back_iter/2,back_iter)),]
    }
    
    Xbeta <- X%*%beta
    Z <- rtruncnorm(n, mean= Xbeta, a = lb, b = ub)

    XRX <- t(X)%*% X
    XRZ <- t(X) %*% Z
    m<- solve(XRX, XRZ)
    cholvar<- t(chol(solve(XRX)))
    beta<- cholvar%*% rnorm(p)+m
    
    # if(i> back_iter){
      trace_beta<- rbind(trace_beta,t(beta))
    # }

    print(i)
  }
  
  list('beta'= trace_beta[back_iter:steps,], 'g'=trace_g)
}
# fit<- probit_r(1,10)

    
fit<-probit_r(1,3000,back_iter = 100)

acf(fit$beta[c(600:1900)], lag.max = 300)

ts.plot(fit$beta)

