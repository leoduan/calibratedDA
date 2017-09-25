require('R.matlab')
setwd("~/git/ImbalancedPG/poisson/")
mp<-readMat("./convdat.mat")

nonzero<-  c(mp$Ntr!=0) & c(mp$Nts!=0)
# 
Ytr<- mp$Ytr[nonzero,]
Ntr<- mp$Ntr[nonzero]
# 
y<- Ytr[,46]
N<- Ntr

setwd("~/git/ImbalancedPG/R-demo/")
source("logit.r")
source("logitRanIntercept.r")

DIC<- function(x){
 -2*sum(x*y - N*log(1+exp(x)))
}

#MH
fit2<- logitRanIntercept(y, N,r_ini= 1,tune = 100,burnin=100, run=1000 ,fixR = F,MH = T,c=2,priorMean = -12,priorVar = 49)
fit2$accept_rate

# save(fit2,file="../../CDA_data/binomialCDA.RDa")

load("../../CDA_data/binomialCDA.RDa")
require("coda")

fit2$beta<- fit2$beta[200:1000,]
fit2$beta0<- fit2$beta0[200:1000]

mean(rowMeans((fit2$beta)))
quantile(rowMeans((fit2$beta)),c(0.025,0.975))

mean(rowMeans((fit2$beta^2)))
quantile(rowMeans((fit2$beta^2)),c(0.025,0.975))

mean(fit2$beta0)
quantile(fit2$beta0,c(0.025,0.975))

mean(fit2$sigma0)
quantile(fit2$sigma0,c(0.025,0.975))

essFit2<-effectiveSize(fit2$beta[,1:2000])/800
# essFit2<- essFit2*1000/800
mean(essFit2)
quantile(essFit2,c(0.025,0.975))


# acf(fit2$beta[,1])
# 
# mean(fit2$sigma0)
# 
# mean(apply(fit2$beta[1:500,], 1,DIC))
# mean(apply(fit2$beta[500:1000,], 1,DIC))


ACFfit2<-apply(fit2$proposal, MARGIN = 2, function(x){c(c(acf(x,plot = F,lag.max = 40))$acf[,,1])})
# ts.plot(ACFfit2[,1:2000])


require("reshape")
require("ggplot2")

CDAacf<- t(ACFfit2)
p<- nrow(CDAacf)

colnames(CDAacf)<- c(0:40)

dfCDAacf <- melt(CDAacf)
colnames(dfCDAacf)<- c("Variable","Lag","ACF")
dfCDAacf$Lag<- as.factor(dfCDAacf$Lag)

pdf("./binomial_random_acf_cda.pdf",5,3)
ggplot(data = dfCDAacf, aes(x=Lag, y=ACF)) + geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(-0.2,1))+ scale_x_discrete(breaks = c(0:8)*5)
dev.off()

acf(fit2$beta0,lag.max = 100)
acf(fit2$sigma0,lag.max = 100)

mean(fit2$beta0[500:1000])
sd(fit2$beta0[500:1000])

mean(fit2$sigma0[500:1000])
sd(fit2$sigma0[500:1000])

#r fixed to 1
fit3<- logitRanIntercept(y, N,r_ini= 1,tune = 200,burnin=2000, run=1000 ,fixR = T,MH = F,c=5,priorMean = -12,priorVar = 49)

# save(fit3,file="../../CDA_data/binomialDA.RDa")
load("../../CDA_data/binomialDA.RDa")


fit3$accept_rate
ACFfit3<-apply(fit3$beta, MARGIN = 2, function(x){c(c(acf(x,plot = F,lag.max = 40))$acf[,,1])})

# ts.plot(fit3$beta[,1])

mean(rowMeans((fit3$beta)))
quantile(rowMeans((fit3$beta)),c(0.025,0.975))

mean(rowMeans((fit3$beta^2)))
quantile(rowMeans((fit3$beta^2)),c(0.025,0.975))



mean(colMeans((fit3$beta)))
quantile(colMeans((fit3$beta)),c(0.025,0.975))

mean(fit3$beta0)
quantile(fit3$beta0,c(0.025,0.975))

mean(fit3$sigma0)
quantile(fit3$sigma0,c(0.025,0.975))

essFit3<-effectiveSize(fit3$beta[,1:2000])/1000

20/ (mean(essFit3)*1000)

mean(essFit3)
quantile(essFit3,c(0.025,0.975))

mean(apply(fit3$beta, 1,DIC))


CDAacf<- t(ACFfit3)
p<- nrow(CDAacf)

colnames(CDAacf)<- c(0:40)

dfCDAacf <- melt(CDAacf)
colnames(dfCDAacf)<- c("Variable","Lag","ACF")
dfCDAacf$Lag<- as.factor(dfCDAacf$Lag)

pdf("./binomial_random_acf_da.pdf",5,3)
ggplot(data = dfCDAacf, aes(x=Lag, y=ACF)) + geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(-0.2,1))+ scale_x_discrete(breaks = c(0:8)*5)
dev.off()

mean(fit3$beta0[500:1000])
sd(fit3$beta0[500:1000])

mean(fit3$sigma0[500:1000])
sd(fit3$sigma0[500:1000])


#No MH
fit4<- logitRanIntercept(y, N,r_ini= 1,tune = 300,burnin=300, run=1000 ,fixR = F,MH = F,c=1,priorMean = -12,priorVar = 49)
ACFfit4<-apply(fit4$beta, MARGIN = 2, function(x){c(c(acf(x,plot = F,lag.max = 40))$acf[,,1])})
ts.plot(ACFfit4)

ts.plot(fit4$beta[,ACFfit4[20,]>0.6])
ts.plot(fit2$beta[,ACFfit4[20,]>0.6])

acf(fit4$beta0,lag.max = 40)
acf(fit4$sigma0,lag.max = 40)


ts.plot(fit2$beta[,ACFfit4[20,]>0.6])


#r fixed to 1
fit5<- logitRanIntercept(y, N,r_ini= 1E-5,tune = 200,burnin=300, run=1000 ,fixR = T,MH = F,c=5,priorMean = -10,priorVar = 10)
fit3$accept_rate
ACFfit3<-apply(fit3$beta, MARGIN = 2, function(x){c(c(acf(x,plot = F,lag.max = 40))$acf[,,1])})
ts.plot(ACFfit3)
acf(fit3$beta0)
acf(fit3$sigma0)



############
fit4$b
max(fit4$r)

ncol(
  
  ACFfit4[,N==max(N)]

ts.plot(fit2$beta[,3])

fit$c

hist(fit$r)


acf(fit$beta,lag.max = 40)
ts.plot(fit$beta[,1])
ts.plot(fit$beta[,2])

ts.plot(fit$proposal[,2])
hist(fit$r)




fitDA<- logitCDA(y,X, N,r_ini= 1,burnin=100, run=500 ,fixR = T)

acf(fitDA$beta,lag.max = 40)
ts.plot(fitDA$beta)

# fitMVN<- logitMVN(y,X, r_ini=1,burnin=2000, run=500,fixR = FALSE)
# ts.plot(fitMVN$beta)



acf(fitMVN$beta,lag.max = 40)

acf_all <- apply(cbind(fitMVN$beta[,1], fitDA$beta[,1], fit$beta[,1]),2,function(x){ acf(x, lag.max = 40, plot=F)$acf})


df<- data.frame("ACF"=c(acf_all),"Method"=rep(c("MH-MVN","DA","CDA"),each=41),"Lag"=rep(c(0:40),3))
df$Method<- ordered(df$Method, levels = c("DA","MH-MVN","CDA"))

pdf("./logit_demo_acf.pdf",5,3)
ggplot(data=df, aes(x=Lag, y=ACF,linetype=Method))+ geom_line(size=.75)
dev.off()


df<- data.frame("Value"=c(fitMVN$beta[,1], fitDA$beta[,1], fit$beta[,1]),"Method"=rep(c("MH-MVN","DA","CDA"),each=500),"Iteration"=rep(c(1:500),3))
df$Method<- ordered(df$Method, levels = c("DA","MH-MVN","CDA"))

pdf("./logit_demo_trace_plot.pdf",5,3)
ggplot(data=df, aes(x=Iteration, y=Value))+ geom_line(size=.75) + facet_grid(Method~.)#, scales = "free_y")
dev.off()



# fitMVN$accept_rate






df<- data.frame("ACF"=c(acf_all),"r"=rep(c("1","10","100","1000","5000"),each=41),"Lag"=rep(c(0:40),5))
require("ggplot2")
pdf("./probit_demo_acf.pdf",5,3)
ggplot(data=df, aes(x=Lag, y=ACF,linetype=r))+ geom_line(size=.75)+ylim(-0.2,1)
dev.off()



acf_all_prop <- apply(cbind(fit1$proposal,fit2$proposal,fit3$proposal,fit4$proposal,fit5$proposal),2,function(x){ acf(x, lag.max = 40, plot=F)$acf})

df<- data.frame("ACF"=c(acf_all_prop),"r"=rep(c("1","10","100","1000","5000"),each=41),"Lag"=rep(c(0:40),5))
pdf("./probit_demo_acf_prop.pdf",5,3)
ggplot(data=df, aes(x=Lag, y=ACF,linetype=r))+ geom_line(size=.75)+ylim(-0.2,1)
dev.off()
