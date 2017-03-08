require("coda")

N<- 1E4


#generate data using Beta= (-8,1)
X0<- 1
X1<- rnorm(N, 1, 1)
X<- cbind(X0,X1)
beta<- c(-8,1)

Xbeta<- X%*%beta
theta<-  exp(Xbeta)/(1+exp(Xbeta))

#generate Bernoulli outcomes
y <- as.numeric(runif(N) < theta)
#total number of positive outcomes
sum(y)

source("logitCDA.r")

#############################################################
#Estimation using the original Polya-Gamma (2013) algorithm
#(equivalent to fixing r to 1 and b to 0 in CDA)
#############################################################
fit.PolyaGamma <-logitCDA(y,X, burnin = 200, run = 500, fixR = T,r_ini = 1)

#traceplot
par(mfrow=c(1,1))
ts.plot(fit.PolyaGamma$beta)

#autocorrelation plot
par(mfrow=c(1,2))
autocorr.plot(fit.PolyaGamma$beta[,1],lag.max =  40,auto.layout=F)
autocorr.plot(fit.PolyaGamma$beta[,2],lag.max =  40,auto.layout=F)


#########################################################################################################################
#Estimation using the calibrated data augmentation (CDA)
#(note: for easier testing on the reviewer's computer, we use the R version of the Polya-Gamma generator,
#which is slightly slower for r_i<1 cases due to the R implementation.
#The C++ implementation does not have this issue.)
##########################################################################################################################


fit.CDA <-logitCDA(y,X, burnin = 200, run = 500)

#traceplot
par(mfrow=c(1,1))
ts.plot(fit.CDA$beta)

#autocorrelation plot
par(mfrow=c(1,2))
autocorr.plot(fit.CDA$beta[,1],lag.max =  40,auto.layout=F)
autocorr.plot(fit.CDA$beta[,2],lag.max =  40,auto.layout=F)
