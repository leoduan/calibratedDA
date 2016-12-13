require("scalableDA")

N<- 1E4


X0<- 1#rnorm(N, 1, 1)
X1<- rnorm(N, 1,1)
X<- cbind(X0,X1)
beta<- c(-5,1)


Xbeta<- X%*%beta
eta0<-  rnorm(N,0,sd= 1)
eta0[1]<- 0
theta<-  exp(Xbeta + eta0)


y<- rpois(N,theta)

hist(y)

sum(y)
n<- N

setwd("~/git/ImbalancedPG/poisson/")
source("poisson.r")
fit2<- poissonCDA(y = y,X = X, r_ini= 1,tune = 100,burnin=100, run=1000,fixR = F,randomEffect =T,sigma2 = 1, bigR=100)
ts.plot(fit2$beta)
print(fit2$accept_rate)
acf(fit2$beta,lag.max = 40)
ts.plot(fit2$sigma2)

