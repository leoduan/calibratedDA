setwd("~/git/ImbalancedPG/R-demo/")
n<- 1E4

X0<- 1#rnorm(N, 1, 1)
X1<- rnorm(n, 1, 1)
X<- cbind(X0,X1)
# beta<- c(-20,1)


# Xbeta<- X%*%beta
Xbeta<- rep(-19,n)
theta<-  exp(Xbeta)/(1+exp(Xbeta))

N<- rep(1E6,n)

y<- rbinom(n,N,theta)
# y<- as.numeric(runif(n)<theta)
sum(y)

fit<- logitCDA(y,as.matrix(X), N,r_ini= 1,tune = 200,burnin=300, run=1000 ,fixR = F,MH = T,c=8)

setwd("~/git/ImbalancedPG/R-demo/")
source("logit.r")
source("logitRanIntercept.r")

#
fit2<- logitRanIntercept(y, N,r_ini= 1,tune = 200,burnin=300, run=1000 ,fixR = F,MH = T,c=10,priorMean = -10,priorVar = 10)
fit2$accept_rate

ACFfit2<-apply(fit2$beta, MARGIN = 2, function(x){c(c(acf(x,plot = F,lag.max = 40))$acf[,,1])})

ts.plot(ACFfit2)

#r fixed to 1
fit3<- logitRanIntercept(y, N,r_ini= 1,tune = 200,burnin=300, run=1000 ,fixR = T,MH = F,c=5,priorMean = -10,priorVar = 10)
fit3$accept_rate
ACFfit3<-apply(fit3$beta, MARGIN = 2, function(x){c(c(acf(x,plot = F,lag.max = 40))$acf[,,1])})
ts.plot(ACFfit3)

#No MH
fit4<- logitRanIntercept(y, N,r_ini= 1,tune = 200,burnin=300, run=1000 ,fixR = F,MH = F,c=2,priorMean = -10,priorVar = 10)
# hist(fit4$r)
# fit4$b
ACFfit4<-apply(fit4$beta, MARGIN = 2, function(x){c(c(acf(x,plot = F,lag.max = 40))$acf[,,1])})
ts.plot(ACFfit4)

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
