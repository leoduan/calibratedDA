require("ImbalancedPG")


N<- 1000

X1<- rnorm(N, 1, 1)
X<- cbind(1,X1)
beta<- c(-5,2)

theta<-  exp(X%*%beta)

y<- rpois(N,theta)
hist(y)
# sum(y==0)

require("rjags")


setwd("~/git/ImbalancedPG/test/")

jags <- jags.model('zip.jags',
                   data = list('X' = X1,
                               'y' = y,
                               'N' = N),
                   n.chains = 1,
                   n.adapt = 200)

update(jags, 1000)

out<- coda.samples(jags,
             c('a', 'b', 'c','d', 'p'),
             1000)

acf(out[[1]][,1])
acf(out[[1]][,2])
# acf(out[[1]][,3])
# acf(out[[1]][,4])
acf(out[[1]][,5])

