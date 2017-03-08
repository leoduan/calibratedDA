setwd("c:/Users/leo/git/calibratedDA/code-for-submission/")
source("probitCDA.r")
require("coda")

#generate data using Beta= (-5,1,-1)
N <- 1E4

X0 <- 1
X1 <- rnorm(N, 1, 1)
X2 <- rnorm(N, 1, 1)
X <- cbind(X0, X1,X2)
beta <- c(-5, 1,-1)
theta <- pnorm(X %*% beta)

#generate Bernoulli outcomes
y <- as.numeric(runif(N) < theta)
#total number of positive outcomes
sum(y)


#############################################################
#Estimation using the original Albert-Chib (1993) algorithm
#(equivalent to fixing r to 1 and b to 0 in CDA)
#############################################################
fit.AlbertChib <-probitCDA(y,X, burnin = 500, run = 1000, fixR = T,r_ini = 1)

#traceplot
par(mfrow=c(1,1))
ts.plot(fit.AlbertChib$beta)

#autocorrelation plot
par(mfrow=c(1,3))
autocorr.plot(fit.AlbertChib$beta[,1],lag.max =  40,auto.layout=F)
autocorr.plot(fit.AlbertChib$beta[,2],lag.max =  40,auto.layout=F)
autocorr.plot(fit.AlbertChib$beta[,3],lag.max =  40,auto.layout=F)


#############################################################
#Estimation using the calibrated data augmentation (CDA)
#############################################################

fit.CDA <-probitCDA(y,X, burnin = 500, run = 1000)

#traceplot
par(mfrow=c(1,1))
ts.plot(fit.CDA$beta)

#autocorrelation plot
par(mfrow=c(1,3))
autocorr.plot(fit.CDA$beta[,1],lag.max =  40,auto.layout=F)
autocorr.plot(fit.CDA$beta[,2],lag.max =  40,auto.layout=F)
autocorr.plot(fit.CDA$beta[,3],lag.max =  40,auto.layout=F)
