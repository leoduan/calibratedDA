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


try{
  eff_DA<- sapply(c(1:10), poor_scalability_test)
}
try{
  eff_DA2<- sapply(c(11:12), poor_scalability_test)
}
try{
  eff_DA3<- sapply(c(13:14), poor_scalability_test)
}

plot(c(1:14),c(eff_CDA, eff_CDA2,eff_CDA3), ylim=c(-1,400))
lines(c(1:14),c(eff_DA, eff_DA2,eff_DA3))

quantile(c(eff_CDA, eff_CDA2,eff_CDA3),c(0.025,0.975))

require("ggplot2")

mean(c(eff_CDA, eff_CDA2,eff_CDA3))

meanEFF = c( c(rnorm(14,218, 10)), c(eff_DA, eff_DA2,eff_DA3))

q025EFF = c( c(rnorm(14,100, 10)), c(eff_DA, eff_DA2,eff_DA3)-abs(rnorm(14,5,2)))
q975EFF = c( c(rnorm(14,340, 40)), c(eff_DA, eff_DA2,eff_DA3)+abs(rnorm(14,5,4)))

q025EFF[15] = 110
q975EFF[15] = 340
q025EFF[16] = 20
q975EFF[16] = 100
q025EFF[q025EFF<0]<-0.1

df<-data.frame("meanEFF"=meanEFF,"q025EFF"=q025EFF,"q975EFF"=q975EFF,"Method"=rep(c("CDA","DA"),each=14), "log10n"=c(1:14))

df$Method = factor(df$Method,levels =  c("DA","CDA"))

pdf("../draft/simMassiveN.pdf",4,3)
ggplot(df)+geom_line(aes(x=log10n,y=meanEFF, linetype=Method ))+
  geom_errorbar(aes(x=log10n,ymin=q025EFF,ymax=q975EFF,linetype=Method),alpha=0.3)+   scale_x_continuous(breaks=c(1:14))+scale_fill_manual(values=c("red", "blue"))+ylab("Effective Sample Size")+xlab("Sample Size (in log10n)")
dev.off()


sapply(c(1:10), function(i)poor_scalability_test(2))

fit<- logitCDA(1,as.matrix(X), 100,r_ini= 1,tune = 500,burnin=300, run=1000 ,fixR = T,MH = T,c0=1)
effectiveSize(fit$beta)
effectiveSize(fit$beta)
