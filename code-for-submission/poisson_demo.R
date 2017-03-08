setwd("c:/Users/leo/git/calibratedDA/code-for-submission/")

require("coda")

N<- 1E4


#generate data using Beta= (-8,1)
X0<- 1
X1<- rnorm(N, 1, 1)
X<- cbind(X0,X1)
beta<- c(-6,3)

Xbeta<- X%*%beta
theta<-  exp(Xbeta)

#generate Poisson outcomes
y <- rpois(N,theta)
hist(y, breaks = 30)
hist(log(y+1),breaks = 30)

source("poissonCDA.r")

#########################################################################################################################
#Estimation with approximate MCMC using Polya-Gamma algorithm proposed by Zhou et al. (2015)
#########################################################################################################################
fit.PolyaGamma <-poissonCDA(y,X, burnin = 200, run = 1000, fixR = T,r_ini = 1, lambda = 1000,MH = F)

#traceplot
par(mfrow=c(1,1))
ts.plot(fit.PolyaGamma$beta)

#autocorrelation plot
par(mfrow=c(1,2))
autocorr.plot(fit.PolyaGamma$beta[,1],lag.max =  100,auto.layout=F)
autocorr.plot(fit.PolyaGamma$beta[,2],lag.max =  100,auto.layout=F)



#########################################################################################################################
#Estimation using the calibrated data augmentation (CDA)
##########################################################################################################################


fit.CDA <- poissonCDA(y,X, burnin = 200, run = 1000)

#traceplot
par(mfrow=c(1,1))
ts.plot(fit.CDA$beta)

#autocorrelation plot
par(mfrow=c(1,2))
autocorr.plot(fit.CDA$beta[,1],lag.max =  100,auto.layout=F)
autocorr.plot(fit.CDA$beta[,2],lag.max =  100,auto.layout=F)
