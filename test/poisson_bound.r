require("ggplot2")

theta<- 10

r_i<- ceiling(exp(theta*1.5)*10)

p<- exp(theta)/r_i / (1+exp(theta)/r_i)

y_poi <- rpois(10000,lambda = exp(theta))
y_bin<- rbinom(10000, r_i,  p)
y_norm<- rnorm(10000, exp(theta),  sqrt(exp(theta)))



mean(y_poi)
mean(y_bin)
mean(y_norm)

sd(y_poi)
sd(y_bin)
sd(y_norm)

hist(y_poi)
hist(y_bin)

dataHist<- data.frame("Dis"=rep(c("Poi","Bin"),each=10000),
                      "RV"= c(y_poi,y_bin))

ggplot(dataHist, aes(RV, fill = as.factor(Dis))) + geom_histogram(alpha = 0.2,bins = 100,  position="identity")



dataHist<- data.frame("Dis"=rep(c("Poi","Normal"),each=10000),
                      "RV"= c(y_poi,y_norm))

ggplot(dataHist, aes(RV, fill = as.factor(Dis))) + geom_histogram(alpha = 0.2,bins = 100,  position="identity")


exp(theta)
r_i * log(1+exp(theta)/r_i)


maxDist<- function(x){
  exp( exp(x)/10 -exp(x))-exp( -exp(x))
}

maxDist(c(1:100))

dist<- function(x){
  ( dpois(x,exp(5)) - dnorm(x,exp(5),exp(5/2)))
}


max(abs(dist(1:(exp(5)*2))))






pg_mean<- function(r, psi){
  r/(abs(psi)-log(r))*tanh(abs(psi)-log(r))
}

pg_mean(0.01,-0.1)
pg_mean(1,-0.1)



gamma(4+1)/(gamma(2+1)*gamma(2+1))
