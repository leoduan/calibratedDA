library(devtools)

setwd("~/git/")
build('ImbalancedPG/')
install.packages("ImbalancedPG_1.0.tar.gz", repos = NULL, type = "source")

require(ImbalancedPG)

y<- 1234
n<- 1E8
r<- 10
x<- ImbalancedPG::multinomial(rep(n,100),rep(y,100), r_ratio = r)
acf(x$theta[,1],lag.max = 100)

theta_mean<- colMeans(x$theta)
mean( exp((theta_mean))/(1+exp(theta_mean)))

sum(exp(x$theta)/ (r*y/n)> 1)

