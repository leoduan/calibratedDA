require("ImbalancedPG")

N<- 1E4

X0<- 1#rnorm(N, 1, 1)
X1<- rnorm(N, 0, 1)
X<- cbind(X0,X1)
beta<- c(1,-5)

Xbeta<- X%*%beta
theta<-  exp(Xbeta)/(1+exp(Xbeta))


y<- as.numeric(runif(N)<theta)
sum(y)
B<- diag(1000,2,2)
b<- rep(0,2)

fit<- ImbalancedPG::logit_reg_simple(y , X,b,B,r0 = 1.5)

fit$beta
acf(fit$beta[,1])
acf(fit$beta[,2])

# hist(fit$beta[,1])
# hist(fit$beta[,2])

sd(fit$beta[,1])
sd(fit$beta[,2])

ts.plot(fit$beta)


require("BayesLogit")
fit2<- logit(y,X) 

acf(fit2$beta[,1])
acf(fit2$beta[,2])

sd(fit2$beta[,1])
sd(fit2$beta[,2])

ts.plot(fit2$beta)

hist(fit$w,xlim=c(0,1))
hist(fit2$w,xlim=c(0,1))

sd(fit$w[,3])
sd(fit2$w[,3])


Xbeta <-fit2$beta %*% t(X)
Xbeta2 <-fit2$beta %*% t(X)

median(1/fit$w[,3])
median(1/fit2$w[,3])

sd(1/fit2$w[,1])
fit$beta

#binomial
fit3<- logit(1,n = 1, X=1)
acf(fit3$w)

