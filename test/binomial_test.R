require("ImbalancedPG")

N<- 1E6
n<- rep(N,2)
theta<- 1
p<- 1/(1+exp(-theta))

y<- rbinom(2,N,p)

fit<- ImbalancedPG::binomial_simple(y,n, r0 = N,burnin = 500,run = 500)

fit$theta

acf(fit$theta[,1])
acf(fit$theta[,2])

ts.plot(fit$theta)
pM<- colMeans(fit$theta)

pM
