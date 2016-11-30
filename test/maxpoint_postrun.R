# .rs.restartR()

require("scalableDA")
setwd("~/git/ImbalancedPG/test/")

source("maxpoint_data.r")

Xinv<-  solve(t(X)%*%X,t(X))

max(abs(Xinv))
sum(abs(Xinv))


load("maxpoint_fit1.rda")
load("maxpoint_fit2.rda")

ts.plot(fit1$beta[,1])
ts.plot(fit2$beta[,1])


library(latex2exp)
library(reshape)
library(ggplot2)



df1<-data.frame("Value"= c( fit1$beta[1001:2000,1], fit1$beta[1001:2000,9], fit1$beta[1001:2000,46] ),
                "Iteration" = rep(c(1:1000),3),
                "Parameter" = as.factor( rep(c("beta.0","beta.8","beta.45"),each=1000)))


pdf("./traceplot_poisson_da.pdf",6,6)
ggplot(data=df1, aes(x=Iteration, y=Value))+ geom_line(size=.75) + facet_grid(Parameter~., scale="free_y")
dev.off()

df2<-data.frame("Value"= c( fit2$beta[,1], fit2$beta[,9], fit2$beta[,46] ),
           "Iteration" = rep(c(1:1000),3),
           "Parameter" = as.factor( rep(c("beta.0","beta.8","beta.45"),each=1000)))


pdf("./traceplot_poisson_cda.pdf",6,6)
ggplot(data=df2, aes(x=Iteration, y=Value))+ geom_line(size=.75) + facet_grid(Parameter~., scale="free_y")
dev.off()



######acf

getAcf<-function(x){
  c(acf(x, lag.max = 40,plot = F)$acf)
}



acfADA<- apply(fit2$beta, 2, getAcf)

df<- melt(acfADA)

df$X1<- df$X1-1
df$X2<- as.factor(df$X2)

names(df)<- c("Lag","ParIndex","ACF")

pdf("./poisson_cda_acf.pdf",6,6)
ggplot(data=df, aes(x=Lag, y=ACF,col= as.factor(ParIndex)))+ geom_line(size=.75)+   theme(legend.position="none")+  scale_colour_manual(values = rep("black",100))
dev.off()




acfADA<- apply(fit1$beta, 2, getAcf)

df<- melt(acfADA)

df$X1<- df$X1-1
df$X2<- as.factor(df$X2)

names(df)<- c("Lag","ParIndex","ACF")

pdf("./poisson_da_acf.pdf",6,6)
ggplot(data=df, aes(x=Lag, y=ACF,col= as.factor(ParIndex)))+ geom_line(size=.75)+   theme(legend.position="none")+  scale_colour_manual(values = rep("black",100))
dev.off()



######

mean(fit1$beta[,1])
sd(fit1$beta[,1])

mean(fit2$beta[,1])
sd(fit2$beta[,1])



mean(rowSums(fit1$beta[,-1]))
sd(rowSums(fit1$beta[,-1]))

mean(rowSums(abs(fit1$beta)))
sd(rowSums(abs(fit1$beta)))

mean(rowSums(abs(fit2$beta)))
sd(rowSums(abs(fit2$beta)))


sqrt( mean( (y-exp(X%*%colMeans(fit1$beta)))^2))
sqrt( mean( (y-exp(X%*%colMeans(fit2$beta)))^2))

rmse<- function(y,mu){
  sqrt( mean( (y-mu)^2))
}

deviance<- function(y,mu){
  a<- y*log(y/mu)
  a[is.na(a)]<-0
  2*sum(a-(y-mu))
}

deviance(y, exp(X%*%colMeans(fit1$beta)))
deviance(y, exp(X%*%colMeans(fit2$beta)))


rmse(y2, exp(X2%*%colMeans(fit1$beta)))
rmse(y2, exp(X2%*%colMeans(fit2$beta)))




pdf("./poisson_fitting_da.pdf",4,4)
plot(exp(X%*%colMeans(fit1$beta)),y,xlim = c(0,800),ylim=c(0,800),xlab="Fitted",ylab="True Value")
abline(a=0,b=1,lty=2)
dev.off()

pdf("./poisson_fitting_ada.pdf",4,4)
plot(exp(X%*%colMeans(fit2$beta)),y,xlim = c(0,800),ylim=c(0,800),xlab="Fitted",ylab="True Value")
abline(a=0,b=1,lty=2)
dev.off()


pdf("./poisson_cv_da.pdf",4,4)
plot(exp(X2%*%colMeans(fit1$beta)),y2,xlim = c(0,800),ylim=c(0,800),xlab="Prediction",ylab="True Value")
abline(a=0,b=1,lty=2)
dev.off()

pdf("./poisson_cv_ada.pdf",4,4)
plot(exp(X2%*%colMeans(fit2$beta)),y2,xlim = c(0,800),ylim=c(0,800),xlab="Prediction",ylab="True Value")
abline(a=0,b=1,lty=2)
dev.off()



hist(abs(colMeans(fit2$beta)[-1]))

ts.plot((fit2$beta[,1]))
ts.plot(rowSums(fit2$beta[,-1]))

mean(fit2$beta[,1])
sd(fit2$beta[,1])

mean(rowSums(fit2$beta[,-1]))
sd(rowSums(fit2$beta[,-1]))

acf(fit2_beta_trace[,1],lag.max = 40)



acf(rowSums(fit2$beta[,-1]),lag.max = 40)



sum(y>0)/N

beta<-colMeans(fit$beta)

plot(exp(X%*%beta)[y<1000],y[y<1000])

plot(exp(fit$theta[10,])[y<1000],y[y<1000])

plot(exp(X2%*%beta)[y2<1000],y2[y2<1000])



pdf("./poisson_beta0.pdf",6,6)
par(mfrow=c(3,1))
ts.plot(fit2_beta_trace[,1],xlab="Iteration",ylab=TeX('$\\beta_0$'))
ts.plot(fit2_beta_trace[,2],xlab="Iteration",ylab=TeX('$\\beta_1$'))
ts.plot(fit2_beta_trace[,3],xlab="Iteration",ylab=TeX('$\\beta_2$'))
dev.off()




plot(exp(fit$theta[1000,]),y,xlim=range(y))

beta<-(fit$beta[1000,])
plot(exp(X%*%beta),y,xlim=range(y))

beta<-colMeans(fit2_beta_trace)

plot((X%*%beta),log(y+0.1),xlim=c(0,10), ylim=c(0,10))
plot((X2%*%beta),log(y2+0.1),xlim=c(0,10), ylim=c(0,10))


plot((X2%*%beta),log(y2+1),xlim=c(0,10), ylim=c(0,10))



load(file="stan_fit_poisson.Rda")
hmc_beta<- matrix(0,2000,p)

for(i in 1:p){
  hmc_beta[,i]<-  eval( parse(text= paste("fit@sim$samples[[1]]$`beta[", i,"]`" ,sep="")))
}
hmc_beta<- hmc_beta[1001:2000,]


cor(fit2$beta)

hmc_beta

mean(rowSums(abs(hmc_beta)))
sd(rowSums(abs(hmc_beta)))

plot(colMeans(hmc_beta[,-1]), colMeans(fit2$beta[,-1]))


sqrt(mean((colMeans(hmc_beta[,-1])- colMeans(fit2$beta[,-1]))^2))




colSD<- function(x){
  apply(x, 2, function(t){sd(t)})
}

sqrt(mean((colSD(hmc_beta[,-1])- colSD(fit2$beta[,-1]))^2))



df3<-data.frame("Value"= c( hmc_beta[,1], hmc_beta[,9], hmc_beta[,46] ),
                "Iteration" = rep(c(1:1000),3),
                "Parameter" = as.factor( rep(c("beta.0","beta.8","beta.45"),each=1000)))


pdf("./traceplot_poisson_hmc.pdf",6,6)
ggplot(data=df3, aes(x=Iteration, y=Value))+ geom_line(size=.75) + facet_grid(Parameter~., scale="free_y")
dev.off()