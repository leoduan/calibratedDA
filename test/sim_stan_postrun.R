setwd("~/git/ImbalancedPG/test/")

library(rstan) # observe startup messages
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("maxpoint_data.r")

#p<-3
#N<-1000

X<- as.matrix(X[1:N,1:p])
y<-y[1:N]
data = list('X' = X,
            'y' = y,
            'N' = N,
            'p' = p
)
# 
# X2<- t(X)%*%X
# diag(X2)<- diag(X2)+1E-3
# b0<- solve(X2, t(X)%*%log(y+1))
# save(fit,file="stan_fit_poisson.Rda")
load(file="stan_fit_poisson.Rda")


for(i in 1:N){
  eval( parse(text= paste("fit@sim$samples[[1]]$`mu[", i,"]`<-NULL" ,sep="")))
  print(i)
}

hmc_beta<- matrix(0,2000,p)

for(i in 1:p){
  hmc_beta[,i]<-  eval( parse(text= paste("fit@sim$samples[[1]]$`beta[", i,"]`" ,sep="")))
}



library(latex2exp)
library(reshape)
library(ggplot2)

pdf("./traceplot_poisson_hmc.pdf",6,6)
par(mfrow=c(3,1))
ts.plot((hmc_beta[1001:2000,1]),xlab="Iteration",ylab=TeX('$\\beta_0$'))
ts.plot((hmc_beta[1001:2000,9]),xlab="Iteration",ylab=TeX('$\\beta_8$'))
ts.plot((hmc_beta[1001:2000,46]),xlab="Iteration",ylab=TeX('$\\beta_45$'))
dev.off()



getAcf<-function(x){
  c(acf(x, lag.max = 40,plot = F)$acf)
}



acfADA<- apply(hmc_beta[1001:2000,], 2, getAcf)

df<- melt(acfADA)

df$X1<- df$X1-1
df$X2<- as.factor(df$X2)

names(df)<- c("Lag","ParIndex","ACF")

pdf("./poisson_hmc_acf.pdf",3,3)
ggplot(data=df, aes(x=Lag, y=ACF,col= as.factor(ParIndex)))+ geom_line(size=.75)+   theme(legend.position="none")+  scale_colour_manual(values = rep("black",100))
dev.off()


load("maxpoint_fit2.rda")


hmc_beta<- hmc_beta[1001:2000,]

mean(hmc_beta[,1])
sd(hmc_beta[,1])


mean(rowSums(hmc_beta[,-1]))
sd(rowSums(hmc_beta[,-1]))



rmse<- function(y,mu){
  sqrt( mean( (y-mu)^2))
}

deviance<- function(y,mu){
  a<- y*log(y/mu)
  a[is.na(a)]<-0
  2*sum(a-(y-mu))
}


rmse(y, exp(X%*%colMeans(hmc_beta)))


deviance(y,  exp(X%*%colMeans(hmc_beta)))

rmse(y2, exp(X2%*%colMeans(hmc_beta)))

fit2_thin<-fit2$beta[c(1:1000),]
rmse(y2, exp(X2%*%colMeans(fit2_thin)))




plot( colMeans(fit2$beta[,-1]), colMeans(hmc_beta[,-1]),xlim=c(-0.2,0.5),ylim=c(-0.2,0.5))

pdf("ADAvsHMC_mean.pdf",6,6)
plot( colMeans(fit2$beta[,-1]), colMeans(hmc_beta[,-1]),
      xlab="ADA",ylab="HMC"
      )
dev.off()


pdf("ADAvsHMC_sd.pdf",6,6)
plot( apply(fit2$beta[,-1],MARGIN = 2,sd), apply(hmc_beta[,-1],MARGIN = 2,sd)    ,xlim=c(0,0.04),ylim=c(0,0.04),
      xlab="ADA",ylab="HMC"
      )
dev.off()



rmse( colMeans(fit2$beta[,-1]), colMeans(hmc_beta[,-1])) /  mean(colMeans(fit2$beta[,-1]))


load("maxpoint_fit1.rda")

fit1_thin<-fit1$beta[(100+c(1:100))*10,]

rmse(y, exp(X%*%colMeans(fit1_thin)))

ts.plot(fit1_thin)
