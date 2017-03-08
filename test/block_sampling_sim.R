require("scalableDA")
require("coda")
# .rs.restartR()
setwd("~/git/ImbalancedPG/poisson/")

n<- 1E3
X<- cbind(1,rnorm(n),rnorm(n))

beta<- c(-5,1,-1)

y<-rpois(n,exp(X%*%beta + rnorm(n,sd = 1)))

fit3<- scalableDA::poisson_reg_block_random(y = y,X = as.matrix(X[,-1]), r_ini= 1,tune = 1,burnin=500, run=500,fixR = F,MH = T,C=1E9,c_ini = 1,adaptC = F,nu_ini = 1, fixNu=F,sigma2 = 100,centeredRanEff = T)

fit4<- scalableDA::poisson_reg_block_random(y = y,X = as.matrix(X), r_ini= 1,tune = 1,burnin=500, run=500,fixR = F,MH = T,C=1E9,c_ini = 10,adaptC = F,nu_ini = 1, fixNu=F,sigma2 = 100,centeredRanEff = F)

acf(fit3$beta[101:500,],lag.max = 40)
acf(fit4$beta[101:500,2:3],lag.max = 40)

ts.plot(fit3$beta)
ts.plot(fit4$beta[,2:3])

ts.plot(fit3$nu)
ts.plot(fit3$tau[,1])
ts.plot(fit3$tau[,2])
ts.plot(fit3$tau[,3])
ts.plot(fit3$tau[,4])

ts.plot(fit3$eta0)



acf(fit3$nu,lag.max = 100)

# fit3$r
