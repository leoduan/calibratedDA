# .rs.restartR()

require("scalableDA")
require("coda")
setwd("~/git/ImbalancedPG/poisson/")
source("maxpoint_data.r")

source("poisson.r")

X<-X

# fit2<- scalableDA::poisson_reg(y = y,X = X, r_ini= 1,tune = 100,burnin=100, run=100,fixR = F,MH = T,C=1E9,c_ini = 5,adaptC = F,randomEff = F,nu = 10)

fit3<- scalableDA::poisson_reg_block_random(y = y,X = as.matrix(X[,-1]), r_ini= 1,tune = 50,burnin=50, run=500,fixR = F,MH = T,C=1E9,c_ini = 50, c_ini2=1,adaptC = F,nu_ini = 3, fixNu=F,sigma2 = 100,centeredRanEff = T)

fit4<- scalableDA::poisson_reg_block_random(y = y,X = as.matrix(X[,-1]), r_ini= 1,tune = 1,burnin=1, run=10,fixR = T,MH = F,C=1E4,c_ini = 1,adaptC = F,nu_ini = 3, fixNu=F,sigma2 = 100,centeredRanEff = T)
  
# fit3<- scalableDA::poisson_reg_block_random(y = y,X = X, r_ini= 1,tune = 1,burnin=1, run=100,fixR = T,MH = T,C=1E9,c_ini = 0.1,adaptC = F,nu_ini = 0.00001, fixNu=T,sigma2 = 100)


acf(fit3$beta[,1],lag.max = 40)
acf(fit4$beta[,1],lag.max = 40)

acf(fit3$beta,lag.max = 40)
ts.plot(fit4$beta)

fit3$c
ts.plot(fit3$beta[300:500,])

hist(X%*%fit3$beta[100,]+ fit3$tau[100,])

# plot(exp(X%*%colMeans(fit2$beta) +colMeans(fit2$tau)),y)
plot(exp(X%*%colMeans(fit3$beta) +colMeans(fit3$tau)),y)

# acf(fit2$beta[,1:5],lag.max = 40)
# acf(fit3$beta[300:500,1:5],lag.max = 40)

ts.plot(fit3$nu)
acf(fit3$nu)

fit3$c
fit3$acceptance_rate

acf(fit3$tau[,1:5],lag.max = 100)

hist(fit3$r)


fit3$c


DIC<- function(x){
  xbeta<- X%*%x
  -2*sum( y*xbeta - exp(xbeta))
}

mean(apply(fit3$beta[1:500,], 1,DIC))
mean(apply(fit3$beta[201:500,], 1,DIC))

essFit3<-effectiveSize(fit3$beta[300:500,])/201
mean(essFit3)
quantile(essFit3,c(0.025,0.975))

13*60/(mean(essFit3)*201)


fit3$acceptance_rate
ts.plot(fit3$beta)

# fit3$beta<- fit3$beta[c(1:500)*2,]



load("../../CDA_data/poissonCDA.RDa")
load("../../CDA_data/poissonDA.RDa")

mean(rowMeans(fit3$beta))
quantile(rowMeans(fit3$beta),c(0.025,0.975))


mean(rowMeans(fit3$tau))
quantile(rowMeans(fit3$tau),c(0.025,0.975))

mean(rowMeans(fit3$tau^2))
quantile(rowMeans(fit3$tau^2),c(0.025,0.975))


mean(rowMeans(fit3$beta^2))
quantile(rowMeans(fit3$beta^2),c(0.025,0.975))

mean(rowMeans(fit4$beta))
quantile(rowMeans(fit4$beta),c(0.025,0.975))


mean(rowMeans(fit4$beta^2))
quantile(rowMeans(fit4$beta^2),c(0.025,0.975))

mean(rowMeans(fit4$tau))
quantile(rowMeans(fit4$tau),c(0.025,0.975))

mean(rowMeans(fit4$tau^2))
quantile(rowMeans(fit4$tau^2),c(0.025,0.975))



mean(rowMeans(fit3$tau))
quantile(rowMeans(fit3$tau),c(0.025,0.975))


mean(rowMeans(fit3$tau^2))
quantile(rowMeans(fit3$tau^2),c(0.025,0.975))
# save(fit3,file="../../CDA_data/poissonCDA.RDa")
# save(fit4,file="../../CDA_data/poissonDA.RDa")

require("reshape")
require("ggplot2")

ACFfit3<-apply(fit3$beta[300:500,], MARGIN = 2, function(x){c(c(acf(x,plot = F,lag.max = 40))$acf[,,1])})
ts.plot(ACFfit3)

CDAacf<- t(ACFfit3)
p<- nrow(CDAacf)

colnames(CDAacf)<-  c(0:40)

dfCDAacf <- melt(CDAacf)
colnames(dfCDAacf)<- c("Variable","Lag","ACF")
dfCDAacf$Lag<- as.factor(dfCDAacf$Lag)

pdf("./poisson_acf_cda.pdf",5,3)
ggplot(data = dfCDAacf, aes(x=Lag, y=ACF)) + geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(-0.2,1))+ scale_x_discrete(breaks = c(0:8)*5)
dev.off()


# fit4<- scalableDA::poisson_reg(y = y,X = X, r_ini= 1,tune = 200,burnin=200, run=1000,fixR = T,MH = T,C=1E4, c=1,)

fit4<- scalableDA::poisson_reg_block_random(y = y,X = as.matrix(X[,-1]), r_ini= 1,tune = 1,burnin=1, run=500,fixR = T,MH = T,C=1E9,c_ini = 30,adaptC = F,nu_ini = 3, fixNu=F,sigma2 = 100,centeredRanEff = T)


essFit4<-effectiveSize(fit4$beta)/1000
mean(essFit4)
quantile(essFit4,c(0.025,0.975))

fit4$acceptance_rate
ts.plot(fit4$beta)

mean(apply(fit4$beta[1:500,], 1,DIC))


ACFfit4<-apply(fit4$beta, MARGIN = 2, function(x){c(c(acf(x,plot = F,lag.max = 100))$acf[,,1])})
ts.plot(ACFfit4,ylim=c(0,1))



CDAacf<- t(ACFfit4)[,c(0:20)*5+1]
p<- nrow(CDAacf)

colnames(CDAacf)<- c(0:20)*5

dfCDAacf <- melt(CDAacf)
colnames(dfCDAacf)<- c("Variable","Lag","ACF")
dfCDAacf$Lag<- as.factor(dfCDAacf$Lag)

pdf("./poisson_acf_da.pdf",5,3)
ggplot(data = dfCDAacf, aes(x=Lag, y=ACF)) + geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(-0.2,1))+ scale_x_discrete(breaks = c(0:20)*5)
dev.off()


