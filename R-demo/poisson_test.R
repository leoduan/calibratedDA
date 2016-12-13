setwd("~/git/ImbalancedPG/R-demo/")
N<- 1E4

X0<- 1#rnorm(N, 1, 1)
X1<- rt(N, 5)
X<- cbind(X0,X1)
beta<- c(-5,1)


Xbeta<- X%*%beta
eta0<-  rnorm(N,0,sd= 1)
theta<-  exp(Xbeta + eta0)


y<- rpois(N,theta)

hist(y)

sum(y)
n<- N

source("poisson.r")

fit<- poissonCDA(y,X, r_ini= 1,burnin=100, run=1000 ,fixR = F,randomEffect = T,sigma2 = 2)

fittedXbeta<-X%*%colMeans(fit$beta)
# plot(fittedXbeta,log(theta))
# plot(fittedXbeta+ fit$eta,log(theta))

ts.plot(fit$beta)
ts.plot(fit$proposal)
acf(fit$beta,lag.max = 40)
acf(fit$proposal,lag.max = 40)
fit$accept_rate

ts.plot(fit$sigma2)
ts.plot(fit$eta)

