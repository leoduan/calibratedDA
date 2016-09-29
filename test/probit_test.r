# .rs.restartR()
require("ImbalancedPG")
setwd("~/git/ImbalancedPG/test/")

N <- 1E4

X0 <- 1#rnorm(N, 1, 1)
X1 <- rnorm(N, 1, 1)
X <- cbind(X0, X1)
beta <- c(-5, 1)

# hist(X %*% beta)
theta <- pnorm(X %*% beta)

a<-numeric()
for(i in 1:10){
  y <- as.numeric(runif(N) < theta)
  a<- c(a,sum(y))
}

mean(a)
sd(a)
     
     b<- rep(0,2)
B<- diag(1000,2)


# fit<- ImbalancedPG::probit_reg_px(y,X,b,B, r0= 1, burnin = 1000, run = 500,nu0 = 0.1)
# 
fit<- ImbalancedPG::probit_reg_simple(y,X,b,B, r0= 100, burnin = 1000, run = 1000)
fit2<- ImbalancedPG::probit_reg_simple(y,X,b,B, r0= 1, burnin = 1000, run = 1000)
fit3<- ImbalancedPG::probit_reg_px(y,X,b,B, r0= 1, burnin = 100, run = 1000,nu0 = 0.1)

colMeans(fit2$beta)
apply(fit2$beta, 2, sd)

colMeans(fit3$beta)
apply(fit3$beta, 2, sd)

colMeans(fit$beta)
apply(fit$beta, 2, sd)

# hist(fit$trace_r0)
# 
# acf(fit$beta[fit$trace_r0>1,1])
# acf(fit$beta[fit$trace_r0>1,2])

hist(fit$r)

library(latex2exp)

pdf("./probit_ada_r.pdf",8,6)
plot(X%*%beta,fit$r, xlab=TeX('$\\beta_0+\\x\\beta_1$'),ylab="r")
dev.off()



colMeans(fit$beta)

sd(fit$beta[,1])
sd(fit$beta[,2])


# sd(fit$beta[,1])

  



table(fit2$r)
sd(fit2$beta[,1])
sd(fit2$beta[,2])



##
par(mfrow=c(3,1))

pdf("./probit_14ts.pdf",8,6)
par(mfrow=c(3,1))
ts.plot(fit2$beta[,1], main="DA",xlab="Iteration",ylab="",lwd=2)
ts.plot(fit3$beta[,1], main="PX-DA",xlab="Iteration",ylab="",lwd=2)
ts.plot(fit$beta_prop[,1]-0.2, main="ADA",xlab="Iteration",ylab="",lwd=2)
dev.off()

# acf(fit$beta[, 1], lag.max = 40)

require("ggplot2")

da<- c(acf(fit2$beta_prop[, 2], lag.max = 40, main="DA",plot = F)$acf)
pxda<-c(acf(fit3$beta[, 2], lag.max = 40, main="DA",plot = F)$acf)
ada<- c(acf(fit$beta[, 2], lag.max = 40, main="DA",plot = F)$acf)

df<- data.frame("ACF"=c(da,pxda,ada),"Method"=rep(c("DA","PX-DA","ADA"),each=41),"Lag"=rep(c(0:40),3))



pdf("./probit_13.pdf",4,3)
ggplot(data=df, aes(x=Lag, y=ACF,linetype=Method))+ geom_line(size=.75)
dev.off()


pdf("./probit_15.pdf",4,3)
ggplot(data=df, aes(x=Lag, y=ACF,linetype=Method))+ geom_line(size=.75)
dev.off()

pdf("./probit_11.pdf",4,3)
ggplot(data=df, aes(x=Lag, y=ACF,linetype=Method))+ geom_line(size=.75)
dev.off()

par(mfrow=c(1,3))

dev.off()

acf(fit$beta_prop[, 2], lag.max = 40)
acf(fit2$beta[, 2], lag.max = 40)
acf(fit3$beta_prop[, 2], lag.max = 40)


acf(fit$beta_prop[, 1]+fit$beta_prop[, 2], lag.max = 40)
acf(fit2$beta[, 1]+fit2$beta[, 2], lag.max = 40)


