# .rs.restartR()
require("scalableDA")
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
     
a     
b<- rep(0,2)
B<- diag(1000,2)


# fit<- ImbalancedPG::probit_reg_px(y,X,b,B, r0= 1, burnin = 1000, run = 500,nu0 = 0.1)
# 
fit<- probit_reg_simple2(y,X,b,B, r0= 200, burnin = 5000, run = 10)
fit2<- ImbalancedPG::probit_reg_simple(y,X,b,B, r0= 1, burnin = 1000, run = 1000)
fit3<- probit_reg_px(y,X,b,B, r0= 1, burnin = 100, run = 1000,nu0 = 0.1)

plot(X%*%beta,fit$r)


colMeans(fit$beta)
acf(fit$beta)



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


setwd("~/work/Dropbox/scalable_da/ms/")

pdf("./probit_cda_r.pdf",5,5)
plot(X%*%beta,fit$r, xlab=TeX('Posterior mean $x_i^T\\beta$'),ylab=TeX('Adapted $r_i'))
dev.off()

pdf("./probit_cda_b.pdf",5,5)
plot(  (sqrt(fit$r)-1)*X%*%colMeans(fit$beta),fit$alpha, xlab=TeX('Posterior mean $( \\sqrt{r_i}-1 ) x_i^T\\beta$'),ylab=TeX('Adapted $b_i'),xlim=c(-80,0),ylim=c(-80,0))
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
ts.plot(fit$beta_prop[,1]-0.2, main="CDA",xlab="Iteration",ylab="",lwd=2)
dev.off()

# acf(fit$beta[, 1], lag.max = 40)

require("ggplot2")

da<- c(acf(fit2$beta[, 2], lag.max = 40, main="DA",plot = F)$acf)
pxda<-c(acf(fit3$beta[, 2], lag.max = 40, main="DA",plot = F)$acf)
cda<- c(acf(fit$beta[, 2], lag.max = 40, main="DA",plot = F)$acf)

df<- data.frame("ACF"=c(da,pxda,cda),"Method"=rep(c("DA (Albert-Chib)","PX-DA","CDA"),each=41),"Lag"=rep(c(0:40),3))

df$Method<- ordered(df$Method, levels = c("DA (Albert-Chib)","PX-DA","CDA"))


pdf("./probit15_acf.pdf",5,3)
ggplot(data=df, aes(x=Lag, y=ACF,linetype=Method))+ geom_line(size=.75)
dev.off()




df<- data.frame("Value"=c(fit2$beta[501:1000, 1],fit3$beta[501:1000, 1], fit$beta[501:1000,1]),"Method"=rep(c("DA (Albert-Chib)","PX-DA","CDA"),each=500),"Iteration"=rep(c(1:500),3))
df$Method<- ordered(df$Method, levels = c("DA (Albert-Chib)","PX-DA","CDA"))

pdf("./probit15_trace_plot.pdf",6,3.6)
ggplot(data=df, aes(x=Iteration, y=Value))+ geom_line(size=.75) + facet_grid(Method~.)
dev.off()
