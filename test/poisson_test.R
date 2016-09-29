# .rs.restartR()

require("ImbalancedPG")

N<- 1E4

X1<- rnorm(N, 0, 1)
X2<- rnorm(N, 0, 1)

X<- cbind(1,X1,X2)
beta<- c(-2,1, 1)

theta<-  exp(X%*%beta + rnorm(N,sd=0.1))

y<- rpois(N,theta)
hist(y)
min(y)


B<- diag(1000,3,3)
b<- rep(0,3)

fit<- ImbalancedPG::poisson_reg(y , X,b,B, r0ini =  10,c = 1,burnin = 1000,run = 1000)

fit<- ImbalancedPG::poisson_reg_random_effect(y , X, r0ini =  10,c = 1,burnin = 1000,run = 1000,da_ver = 1,update_sigma2 = T)

ts.plot(fit$r[,2])
ts.plot(fit$beta)
ts.plot(fit$sigma2)

acf(fit$beta[,1],lag.max = 40)
acf(fit$beta[,2],lag.max = 40)
acf(fit$beta[,3],lag.max = 40)


sd(fit$beta[,1])
sd(fit$beta[,2])
sd(fit$beta[,3])

hist(fit$beta[,1])
hist(fit$beta[,2])
hist(fit$beta[,3])



plot(exp(fit$theta[1000,]),y)

acf(fit$theta[,1],lag.max = 40)
