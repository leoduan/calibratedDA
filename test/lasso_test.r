require('monomvn')

N<-10000
p<- 1000
p1<- 3

X<- matrix(rnorm(N*p),N)

beta<- c(rnorm(p1,0,1),rep(0,p-p1))

y<- rnorm(N,X%*%beta,1)

bl<-bhs(X,y,rao.s2=FALSE,icept = F,thin = 1)


acf(bl$lambda2)


ts.plot(bl$beta[501:1000,1])

ts.plot(bl$beta[,2])
ts.plot(bl$beta[,3])

ts.plot(bl$beta[1000,1:3])
beta


s<- 1/rgamma(N,6/2,6/2)
tsample<- rnorm(N,0,sqrt(s))
hist(tsample)

sd(tsample)
sd(rt(10000,6))


sd(rchisq(10000,2))
sd(rgamma(10000,2/2,rate=1/2))
