require("scalableDA")

N<- 1E4


X0<- 1#rnorm(N, 1, 1)
X1<- rnorm(N, 1,1)
X<- cbind(X0,X1)
beta<- c(-3,1)


Xbeta<- X%*%beta
theta<-  exp(Xbeta)


y<- rpois(N,theta)

hist(y)

sum(y)
n<- N
# 


setwd("~/git/ImbalancedPG/poisson/")
source("poisson.r")
fit2<- poissonCDA(y = y,X = X, r_ini= 1,tune = 300,burnin=100, run=1000,fixR = F,MH = T,C=1E9, c=2)
fit2$accept_rate
# fit2$r
ts.plot(fit2$beta)
# print(fit2$accept_rate)
acf(fit2$beta,lag.max = 40)
# ts.plot(fit2$sigma2)


fit3<- scalableDA::poisson_reg(y = y,X = X, r_ini= 1,tune = 100,burnin=100, run=500,fixR = F,MH = T,C=1E9,c_ini = 2)

DIC<- function(x){
  xbeta<- X%*%x
  -2*sum( y*xbeta - exp(xbeta))
}

mean(apply(fit3$beta[1:500,], 1,DIC))
mean(apply(fit3$beta[201:500,], 1,DIC))



fit3$acceptance_rate
ts.plot(fit3$beta)
# fit3$beta<- fit3$beta[c(1:500)*2,]

ACFfit3<-apply(fit3$beta, MARGIN = 2, function(x){c(c(acf(x,plot = F,lag.max = 40))$acf[,,1])})
ts.plot(ACFfit3)


require("reshape")
require("ggplot2")

CDAacf<- t(ACFfit3)
p<- nrow(CDAacf)

colnames(CDAacf)<- c(0:40)

dfCDAacf <- melt(CDAacf)
colnames(dfCDAacf)<- c("Variable","Lag","ACF")
dfCDAacf$Lag<- as.factor(dfCDAacf$Lag)

pdf("./poisson_acf_cda.pdf",5,3)
ggplot(data = dfCDAacf, aes(x=Lag, y=ACF)) + geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(-0.2,1))+ scale_x_discrete(breaks = c(0:8)*5)
dev.off()


fit4<- scalableDA::poisson_reg(y = y,X = X, r_ini= 1,tune = 200,burnin=200, run=1000,fixR = T,MH = T,C=1E9, c=1)

fit4$acceptance_rate
ts.plot(fit4$beta)

mean(apply(fit4$beta[1:500,], 1,DIC))


ACFfit4<-apply(fit4$beta, MARGIN = 2, function(x){c(c(acf(x,plot = F,lag.max = 40))$acf[,,1])})
ts.plot(ACFfit4,ylim=c(0,1))

CDAacf<- t(ACFfit4)
p<- nrow(CDAacf)

colnames(CDAacf)<- c(0:100)

dfCDAacf <- melt(CDAacf)
colnames(dfCDAacf)<- c("Variable","Lag","ACF")
dfCDAacf$Lag<- as.factor(dfCDAacf$Lag)

pdf("./poisson_acf_da.pdf",5,3)
ggplot(data = dfCDAacf, aes(x=Lag, y=ACF)) + geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(-0.2,1))+ scale_x_discrete(breaks = c(0:8)*5)
dev.off()


