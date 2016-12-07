require("monomvn")

N<- 20
p<- 100

X<- matrix(rnorm(N*p),N)

beta<- c(1,rep(0,p-1))

y<- rnorm(N,X%*%beta, 0.01)


fit<- blasso(X,y,T=1000,RJ=FALSE,rao.s2=F)
mean(fit$beta[500:1000,1])
mean(fit$beta[500:1000,2])
mean(fit$beta[500:1000,3])
