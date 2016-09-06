require("ImbalancedPG")

N<- 1E3

X1<- rnorm(N, 0, 1)
X2<- rnorm(N, 0, 1)

X<- cbind(1,X1,X2)
beta<- c(-3, 0, 0)

theta<-  exp(X%*%beta)
# hist(X%*%beta)

y<- rpois(N,theta)
# hist(y)
min(y)

B<- diag(1000,3,3)
b<- rep(0,3)

fit<- ImbalancedPG::poisson_reg(y , X,b,B, r0ini =  10,c = 2,burnin = 1000,run = 1000)

fit$w[1000,]

ts.plot(fit$beta[500:1000,])


acf(fit$beta[,1])
acf(fit$beta[,2])
acf(fit$beta[,3])

ts.plot(fit$r[,1:50])

hist(fit$r)

sd(fit$beta[,1])
sd(fit$beta[,2])
sd(fit$beta[,3])

hist(fit$beta[,1])
hist(fit$beta[,2])
hist(fit$beta[,3])

.rs.restartR()
rm(list=ls())

