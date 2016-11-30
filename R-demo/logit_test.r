N<- 1E4

X0<- 1#rnorm(N, 1, 1)
X1<- rnorm(N, 1, 1)
X<- cbind(X0,X1)
beta<- c(-5,1)


Xbeta<- X%*%beta
theta<-  exp(Xbeta)/(1+exp(Xbeta))


y<- as.numeric(runif(N)<theta)

sum(y)
n<- N


fit<- logitCDA(y,X, r_ini=0.1,burnin=100, run=100 ,fixR = FALSE)
  

ts.plot(fit$beta)
acf(fit$beta)
