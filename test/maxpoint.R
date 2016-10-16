# .rs.restartR()

require("scalableDA")
setwd("~/git/ImbalancedPG/test/")

source("maxpoint_data.r")

# fit<- ImbalancedPG::poisson_reg_random_effect(y , X, r0ini =  10,c = 1,burnin = 1000,run = 1000,update_sigma2 = T)

# fit<- poisson_reg(y , X,b,B, r0ini =  10,c = 1.1,burnin = 1000,run = 1000,fixed_R = T)

# fit2<- poisson_reg(y , X,b,B, r0ini =  10,c = 1.1,burnin = 1000,run = 1000,fixed_R = F)

# fit<- poisson_reg_random_effect(y , X, tau =  10,c = 1,burnin = 1000,run = 1000,da_ver = 1)


load(file="maxpoint_fit1_random.rda")
load(file="maxpoint_fit2_random.rda")

ts.plot(fit1$beta[,1])
ts.plot(fit1$sigma2)

plot(exp(fit1$theta),y)

ts.plot(fit2$beta[,1])
ts.plot(fit2$sigma2)

acf(fit1$beta[,1],lag.max = 40)
acf(fit2$beta[,1],lag.max = 40)

ts.plot(fit1$sigma2)
ts.plot(fit2$sigma2)

# fit2_beta_trace<- fit2$beta
# save(fit2_beta_trace,file="fit_maxpoint_ada.rda")

# fit<- ImbalancedPG::poisson_reg_random_effect(y , X, r0ini =  10,c = 1,burnin = 10,run = 10,update_sigma2 = T)


load(file="~/fit_maxpoint_ada.rda")


sum(y>0)/N

beta<-colMeans(fit$beta)

plot(exp(X%*%beta)[y<1000],y[y<1000])

plot(exp(fit$theta[10,])[y<1000],y[y<1000])

plot(exp(X2%*%beta)[y2<1000],y2[y2<1000])

library(latex2exp)
library(reshape)
library(ggplot2)

pdf("./poisson_beta0.pdf",6,6)
par(mfrow=c(3,1))
ts.plot(fit2_beta_trace[,1],xlab="Iteration",ylab=TeX('$\\beta_0$'))
ts.plot(fit2_beta_trace[,2],xlab="Iteration",ylab=TeX('$\\beta_1$'))
ts.plot(fit2_beta_trace[,3],xlab="Iteration",ylab=TeX('$\\beta_2$'))
dev.off()



getAcf<-function(x){
  c(acf(x, lag.max = 40,plot = F)$acf)
}

acfADA<- apply(fit2_beta_trace, 2, getAcf)

df<- melt(acfADA)

df$X1<- df$X1-1
df$X2<- as.factor(df$X2)

names(df)<- c("Lag","ParIndex","ACF")

pdf("./poisson_ada_acf.pdf",3,3)
ggplot(data=df, aes(x=Lag, y=ACF,linetype=ParIndex))+ geom_line(size=.75)+   theme(legend.position="none")
dev.off()


plot(exp(fit$theta[1000,]),y,xlim=range(y))

beta<-(fit$beta[1000,])
plot(exp(X%*%beta),y,xlim=range(y))

beta<-colMeans(fit2_beta_trace)

plot((X%*%beta),log(y+0.1),xlim=c(0,10), ylim=c(0,10))
plot((X2%*%beta),log(y2+0.1),xlim=c(0,10), ylim=c(0,10))


plot((X2%*%beta),log(y2+1),xlim=c(0,10), ylim=c(0,10))


