require("ImbalancedPG")


N<- 10000

X1<- rnorm(N, 1, 1)
X<- cbind(1,X1)
beta<- c(-8,3)

theta<-  exp(X%*%beta)

y<- rpois(N,theta)
hist(y)

B<- diag(1000,2,2)
b<- rep(0,2)

fit<- ImbalancedPG::poisson_reg(y , X,b,B,r0 = 10)

ts.plot(fit$beta)

acf(fit$beta[,1])
acf(fit$beta[,2])

sd(fit$beta[,1])
sd(fit$beta[,2])
