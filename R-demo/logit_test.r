setwd("~/git/ImbalancedPG/R-demo/")
N<- 1E4

X0<- 1#rnorm(N, 1, 1)
X1<- rnorm(N, 1, 1)
X<- cbind(X0,X1)
beta<- c(-1,1)


Xbeta<- X%*%beta
theta<-  exp(Xbeta)/(1+exp(Xbeta))


y<- as.numeric(runif(N)<theta)

sum(y)
n<- N

source("logit.r")

fit<- logitCDA(y,X, r_ini= 1,burnin=100, run=500 ,fixR = FALSE)
fitDA<- logitCDA(y,X, r_ini= 1,burnin=100, run=500 ,fixR = T)
fitMVN<- logitMVN(y,X, r_ini=1,burnin=2000, run=500,fixR = FALSE)
ts.plot(fitMVN$beta)



acf(fit$beta,lag.max = 40)
acf(fitDA$beta,lag.max = 40)
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
