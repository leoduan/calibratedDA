require("scalableDA")

N<- 1E4


X0<- 1#rnorm(N, 1, 1)
X1<- rnorm(N, 1,1)
X<- cbind(X0,X1)
beta<- c(-5,1)


Xbeta<- X%*%beta
# eta0<-  rnorm(N,0,sd= 1)
theta<-  exp(Xbeta)


y<- rpois(N,theta)

hist(y)

sum(y)
n<- N

# setwd("~/git/ImbalancedPG/R-demo/")
# source("poisson.r")
# fit<- poissonCDA(y = y,X = X, r_ini= 1,tune = 100,burnin=500, run=1000 ,fixR = F,randomEffect =T,sigma2 = 1)
# ts.plot(fit$beta)
# fit$accept_rate
# acf(fit$beta,lag.max = 40)

fit2<- scalableDA::poisson_reg(y , X,tune = 100,burnin = 100,run = 1000)
ts.plot(fit2$beta)
acf(fit2$beta, lag.max = 40)
