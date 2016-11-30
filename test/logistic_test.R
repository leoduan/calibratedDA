# .rs.restartR()
require("scalableDA")
setwd("~/git/ImbalancedPG/test/")

N<- 1E5

X0<- 1#rnorm(N, 1, 1)
X1<- rnorm(N, 1, 1)
X<- cbind(X0,X1)
# beta<- c(-11,1)
beta<- c(-5,1)


Xbeta<- X%*%beta
theta<-  exp(Xbeta)/(1+exp(Xbeta))

a<-numeric()
for(i in 1:10){
  y<- as.numeric(runif(N)<theta)
  a<- c(a, sum(y))
}

mean(a)/N
sd(a)/N

B<- diag(1000,2,2)
b<- rep(0,2)

fit<- logistic_reg(y , X,b,B,r0 = 10,burnin = 1000,run = 1000)
fit2<- ImbalancedPG::logit_reg_simple(y , X,b,B,r0 = N,burnin = 1000,run = 1000)

colMeans(fit2$beta)
apply(fit2$beta, 2, sd)


colMeans(fit$beta)
apply(fit$beta, 2, sd)




require("ggplot2")
require("reshape")

a<- fit$r

df<- melt(a[,1:20])

pdf("./r_adaptation.pdf",3,3)
ggplot(data=df, aes(x=X1, y=value,col= factor(X2)))+ geom_line(size=0.5)+  scale_colour_manual(values = rep("black",30))+ylab("r")+xlab("Iteration")+  theme(legend.position="none")
dev.off()



df<- data.frame("z"=c((fit$w[1000,]), (fit2$w[1000,])),"Method"=rep(c("CDA","DA"),each=N))

pdf("./logit_z.pdf",4,3)
ggplot(df, aes(x=z, fill=Method, col=Method)) + geom_histogram(alpha=0.2, position="identity", bins=50)+xlim(-1,0.3)+ylim(0,15000)+ scale_fill_grey(start = 0,end = 0.5)+ scale_color_grey(start = 0,end = 0.5)
dev.off()




ada<- c(acf(fit$beta[, 1], lag.max = 40, main="DA",plot = F)$acf)
da<- c(acf(fit2$beta[, 1], lag.max = 40, main="DA",plot = F)$acf)

df<- data.frame("ACF"=c(da,ada),"Method"=rep(c("DA","CDA"),each=41),"Lag"=rep(c(0:40),2))
df$Method<- ordered(df$Method, levels = c("DA","CDA"))



pdf("./logit_11.pdf",4,3)
ggplot(data=df, aes(x=Lag, y=ACF,linetype=Method))+ geom_line(size=.75)
dev.off()

pdf("./logit_15.pdf",4,3)
ggplot(data=df, aes(x=Lag, y=ACF,linetype=Method))+ geom_line(size=.75)
dev.off()

pdf("./logit_19.pdf",4,3)
ggplot(data=df, aes(x=Lag, y=ACF,linetype=Method))+ geom_line(size=.75)
dev.off()





fit$beta
acf(fit$beta[,1])
acf(fit$beta[,2])

# hist(fit$beta[,1])
# hist(fit$beta[,2])

sd(fit$beta[,1])
sd(fit$beta[,2])

ts.plot(fit$beta)


require("BayesLogit")
fit2<- logit(y,X) 

acf(fit2$beta[,1])
acf(fit2$beta[,2])

sd(fit2$beta[,1])
sd(fit2$beta[,2])

ts.plot(fit2$beta)

hist(fit$w,xlim=c(0,1))
hist(fit2$w,xlim=c(0,1))

sd(fit$w[,3])
sd(fit2$w[,3])


Xbeta <-fit2$beta %*% t(X)
Xbeta2 <-fit2$beta %*% t(X)

median(1/fit$w[,3])
median(1/fit2$w[,3])

sd(1/fit2$w[,1])
fit$beta

#binomial
fit3<- logit(1,n = 1, X=1)
acf(fit3$w)

