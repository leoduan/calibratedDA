require("scalableDA")

setwd("~/git/calibratedDA/R-demo/")
n<- 1

X0<- 1#rnorm(N, 1, 1)
X1<- rnorm(n, 1, 1)
X<- cbind(X0,X1)
# beta<- c(-20,1)


# Xbeta<- X%*%beta
Xbeta<- rep(-7,n)
theta<-  exp(Xbeta)/(1+exp(Xbeta))

# require("coda")
# 
# ts.plot(fit$beta)
# acf(fit$beta)

scalability_test <-function(p){
  N<- rep(10^p,n)
  
  # y<- rbinom(n,N,theta)
  # y<- as.numeric(runif(n)<theta)
  sum(y)
  X<- 1
  fit<- logitCDA(1,as.matrix(X), N,r_ini= 1,tune = 500,burnin=300, run=1000 ,fixR = F,MH = T,c0=1)
  
  
  effectiveSize(fit$beta)
}

try{
eff_CDA<- sapply(c(1:10), scalability_test)
}
try{
eff_CDA2<- sapply(c(11:12), scalability_test)
}
try{
eff_CDA3<- sapply(c(13:14), scalability_test)
}

plot(c(6:14),c(eff_CDA, eff_CDA2,eff_CDA3))



poor_scalability_test <-function(p){
  N<- rep(10^p,n)
  
  X<- 1
  fit<- logitCDA(1,as.matrix(X), N,r_ini= 1,tune = 500,burnin=300, run=1000 ,fixR = T,MH = T,c0=1)
  
  
  effectiveSize(fit$beta)
}


eff_DA <- sapply(c(1:10),function(i){
  sapply(c(4:7), poor_scalability_test)
})
eff_CDA<- sapply(c(1:10),function(i){
  sapply(c(4:7), scalability_test)
})

row025<-function(x)apply(x,1,function(y){quantile(y,0.1)})
row975<-function(x)apply(x,1,function(y){quantile(y,0.9)})


require("ggplot2")

mean(c(eff_CDA, eff_CDA2,eff_CDA3))

meanEFF = c(rowMeans(eff_CDA), rowMeans(eff_DA))

q025EFF = c( row025(eff_CDA), row025(eff_DA))
q975EFF = c( row975(eff_CDA), row975(eff_DA))


df<-data.frame("meanEFF"=meanEFF,"q025EFF"=q025EFF,"q975EFF"=q975EFF,"Method"=rep(c("CDA-Subsampling","DA-Subsampling"),each=4), "log10n"=c(5:8))

df$Method = factor(df$Method,levels =  c("DA-Subsampling","CDA-Subsampling"))

pdf("../draft/simMassiveNSubsampling.pdf",5,3)
ggplot(df)+geom_point(aes(x=log10n,y=meanEFF,shape =Method ))+geom_line(aes(x=log10n,y=meanEFF, linetype=Method ))+  geom_errorbar(aes(x=log10n,ymin=q025EFF,ymax=q975EFF,linetype=Method),alpha=0.3)+   scale_x_continuous(breaks=c(1:14))+scale_fill_manual(values=c("red", "blue"))+ylab("Effective Sample Size")+xlab("Sample Size (in log10n)")
dev.off()


sapply(c(1:10), function(i)poor_scalability_test(2))

fit<- logitCDA(1,as.matrix(X), 100,r_ini= 1,tune = 500,burnin=300, run=1000 ,fixR = T,MH = T,c0=1)
effectiveSize(fit$beta)
effectiveSize(fit$beta)
