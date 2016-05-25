library(devtools)

setwd("~/git/")
build('ImbalancedPG/')
install.packages("ImbalancedPG_1.0.tar.gz", repos = NULL, type = "source")

unload("ImbalancedPG")
require("ImbalancedPG")


#logit regression

N<- 10000


X1<- rnorm(N, 1, 0.2)
X1<- (X1 - mean(X1))/sd(X1)

X<- cbind(1,X1)


beta0<- as.vector(c(-8, 2))

theta<- (X%*%beta0)
p<- 1/(1+ exp(-theta))

y0<- as.numeric(runif(N)< p)
sum(y0)

b<- rep(0,2)
B<- diag(1000,2)
test<- ImbalancedPG::logit_reg(y0, X, b, B, mc_draws=1, r_ratio = 10 )

#trace of the percentage of y=0's ignored
ts.plot(test$filter_count)

acf(test$beta[,1])
acf(test$beta[,2])

hist(test$beta[,1])
hist(test$beta[,2])
