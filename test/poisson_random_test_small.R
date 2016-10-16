# .rs.restartR()

require("scalableDA")

N<- 1E3

X1<- rnorm(N, 0, 1)

X<- cbind(1,X1)
beta<- c(-5,1)

theta<-  exp(X%*%beta + rnorm(N,sd= 1))

y<- rpois(N,theta)
hist(y)
min(y)


B<- diag(1000,3,3)
b<- rep(0,3)

# fit1<- poisson_reg(y , X,b,B, r0ini =  10,c = 1,burnin = 1000,run = 1000)

fit2<- poisson_reg_random_effect(y , X, tau =  5,c = 1,burnin = 1000,run = 500,da_ver = 1,max_r = 1000,update_sigma = T)

acf(fit2$beta[1:500,1],lag.max = 40)
acf(fit2$beta[1:500,2],lag.max = 40)

# mean(fit1$r)
mean(fit2$r)

# ts.plot(fit$r[,2])
# ts.plot(fit1$beta)
ts.plot(fit2$beta)

ts.plot(fit2$sigma2)



sd(fit$beta[,1])
sd(fit$beta[,2])
sd(fit$beta[,3])

hist(fit$beta[,1])
hist(fit$beta[,2])
hist(fit$beta[,3])



plot(exp(fit$theta[1000,]),y)

acf(fit$theta[,1],lag.max = 40)
