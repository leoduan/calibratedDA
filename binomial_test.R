library(devtools)

setwd("~/git/")
build('ImbalancedPG/')
install.packages("ImbalancedPG_1.0.tar.gz", repos = NULL, type = "source")

unload("ImbalancedPG")
require("ImbalancedPG")

y<- 12
n<- 1E8
r<- 10
x<- ImbalancedPG::binomial(rep(n,100),rep(y,100), r_ratio = r)
acf(x$theta[,2],lag.max = 20)


theta_mean<- colMeans(x$theta)
hist( exp((theta_mean))/(1+exp(theta_mean)))

sum(exp(x$theta)/ (r*y/n)> 1)



#logit regression

N<- 10000
pos<- 50


X1<- c(c(rnorm(N, 1, 0.2), rnorm(pos, 1.5, 0.2)))
X1<- (X1 - mean(X1))/sd(X1)

X<- cbind(1,X1)


beta0<- as.vector(c(-8, 1))

theta<- (X%*%beta0)
p<- 1/(1+ exp(-theta))

y0<- as.numeric(runif(N+pos)< p)

b<- rep(0,2)
B<- diag(1000,2)
test<- ImbalancedPG::logit_reg(y0, X, b, B, mc_draws=1E2, r_ratio = 2)

sum(test$eta>0)
nrow(test$test)
