require("ImbalancedPG")


N<- 1E3

X1<- rnorm(N, 3, 1)
X2<- rnorm(N, 3, 1)

X<- cbind(1,X1,X2)
beta<- c(-5,1 , 1)

theta<-  exp(X%*%beta)
# hist(X%*%beta)

y<- rpois(N,theta)
# hist(y)
min(y)

B<- diag(1000,3,3)
b<- rep(2,3)

fit<- ImbalancedPG::poisson_reg(y , X,b,B, r0= 10,c = 0.1)

ts.plot(fit$beta)

acf(fit$beta[,1])
acf(fit$beta[,2])
acf(fit$beta[,3])

sd(fit$beta[,1])
sd(fit$beta[,2])
sd(fit$beta[,3])

hist(fit$beta[,1])
hist(fit$beta[,2])
hist(fit$beta[,3])


logDensity1 <- dpois(y, exp(X%*%beta),log = T)


theta<- 10
r<- 100
y<- 10
dpois(y,theta ,log = F) 
dnbinom(y, r, r/(r+theta), log=F)

dpois(2,theta ,log = T) - dnbinom(2, r, r/(r+theta), log=T)

##
approx_den<- function(y,theta,r){
 lgamma(r+1) - lgamma(y+1) - lgamma(r-y+1) + y*log(exp(theta)/r) - r* log(1+ exp(theta)/r)
}

logDensity2<- approx_den(y, X%*%beta,  exp(X%*%beta*2)*10)

sum(logDensity1 - logDensity2)

sum(abs(exp(X%*%beta)/(1+0.001) - exp(X%*%beta)))
median(abs(exp(X%*%beta)/(1+0.1)/(1+0.1) - exp(X%*%beta)))


y<- 2

theta<- 3

comb = function(n, x) {
  return(factorial(n) / (factorial(x) * factorial(n-x)))
}

y<- 7
theta<- 2
r<- exp(theta)*10

p<- dpois(y, exp(theta))/ppois(r,exp(theta))

log(p) - log(q)

y<-100
r<- 100
- lgamma(y+r+1) + y*log(r) + lgamma(r+1)

gamma(y+r+1)/ r^y/gamma(r+1)

# (exp(theta))^y / (1+exp(theta)/1E2)^1E2/factorial(y)


comb(1E2*0.947563, y)*(exp(theta)/1E2)^y / (1+exp(theta)/1E2)^1E2


r<- theta^2*10

sum( lgamma(r+1)-log(r)*y -lgamma(r-y+1) - r*log(1+(theta)/r)) + sum(theta)

