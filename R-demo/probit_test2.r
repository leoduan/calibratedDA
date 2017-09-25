setwd("~/git/calibratedDA/R-demo/")
require("ggplot2")

#generate data
N <- 1E4


#intercept only

y<- rep(0,N)
y[1]<-1

X<- matrix(1,N)
# sum(y)


source("probitPW.r")

#generate data
N <- 1E4

X0 <- 1
X1 <- rnorm(N, 1, 1)
X2 <- rnorm(N, 1, 1)
X <- cbind(X0, X1,X2)
beta <- c(-5, 1,-1)
theta <- pnorm(X %*% beta)

y <- as.numeric(runif(N) < theta)

sum(y)

require("coda")

fit0 <-probitCDA(y,X, burnin = 100, run = 500, fixR = F, MH=T,r_ini = 10,pw_precision = 0)
fit1 <-probitCDA(y,X, burnin = 100, run = 500, fixR = F, MH=T,r_ini = 10,pw_precision = 1)
fit2 <-probitCDA(y,X, burnin = 100, run = 500, fixR = F, MH=T,r_ini = 10,pw_precision = 2)
fit3 <-probitCDA(y,X, burnin = 100, run = 500, fixR = F, MH=T,r_ini = 10,pw_precision = 3)


effectiveSize(fit0$beta[,1])
effectiveSize(fit1$beta[,1])
effectiveSize(fit2$beta[,1])
effectiveSize(fit3$beta[,1])

hist(fit0$trace_accept_ind[,500])
hist(fit1$trace_accept_ind[,500])
hist(fit2$trace_accept_ind[,500])

ar0<-apply(fit0$trace_accept_ind,1,function(x)prod(x))
ar1<-apply(fit1$trace_accept_ind,1,function(x)prod(x))
ar2<-apply(fit2$trace_accept_ind,1,function(x)prod(x))
ar3<-apply(fit3$trace_accept_ind,1,function(x)prod(x))

Xbeta = X%*%fit3$beta[500,]
plot( Xbeta, (Xbeta + fit3$b)/sqrt(fit3$r))


source("probitMAP.r")
fit5 <-probitCDA(y,X, burnin = 200, run = 500, fixR = F, MH=T,r_ini = 10)


Xbeta = X%*%fit5$beta[500,]

pdf("../draft/adaptDiag.pdf",6,6)
plot( Xbeta, (Xbeta + fit5$b)/sqrt(fit5$r),xlab= expression(paste( frac(x * beta + b, sqrt(r)))) ,ylab= expression(paste( x * beta))    ,lwd=0.5,xlim=c(-12,2),
      ylim=c(-12,2))
dev.off()


rmse<- function(a,b){
  sqrt(mean((a-b)^2))
}

plot(Xbeta, (Xbeta + fit5$b)/sqrt(fit5$r))
rmse(Xbeta, (Xbeta + fit5$b)/sqrt(fit5$r))


trace_r = fit5$trace_r
# trace_r[log(trace_r)>30]<- NA

ts.plot((log(trace_r[1:500,99:100])),ylab=expression(paste(log(r) )),ylim=c(0,60))


pdf("../draft/adaptTraceR.pdf",6,6)
ts.plot((log(trace_r[1:400,51:98])),ylab=expression(paste(log(r) )))
dev.off()
# plot( Xbeta, (Xbeta + fit5$b)/sqrt(fit5$r))


mean(ar0)
mean(ar1)

ar0[ar0>1]<-1
ar1[ar1>1]<-1
ar2[ar2>1]<-1
ar3[ar3>1]<-1

mean(ar0)
hist(ar0,breaks = 100)

mean(ar1)
hist(ar1,breaks = 100)

mean(ar2)
hist(ar2,breaks = 100)

mean(ar3)
hist(ar3,breaks = 100)


fit6$accept_rate
ts.plot(fit6$beta)
acf(fit6$beta,lag.max = 100)
colMeans(fit6$beta)


ts.plot(log(fit6$trace_r[,1:50]))



{
  pw_precision=0
  theta_list = seq(1/10^pw_precision,15,length.out = 15* (10^pw_precision))
  theta_list = c(rev(-theta_list),theta_list)
  prob <- pnorm(theta_list)
  dprob<- dnorm(theta_list)
  r_list <-    c(prob*(1 - prob)/dprob^2)
  b_list <- (sqrt(r_list)-1) * theta_list
}

plot(X%*%colMeans(fit6$beta),log(fit6$r))
lines(theta_list,log(r_list),col="red")

ts.plot(log(fit6$trace_r[,1:100]))



plot(X%*%colMeans(fit6$beta),sign(fit6$b)*log(abs(fit6$b)))


colMeans(fit7$beta^2)
sum(exp(fit7$importance)* fit7$beta[,1]^2) /sum(exp(fit7$importance))
sum(exp(fit7$importance)* fit7$beta[,2]^2) /sum(exp(fit7$importance))
sum(exp(fit7$importance)* fit7$beta[,3]^2) /sum(exp(fit7$importance))

# fit6$r[fit6$r>10000]<-NA

plot(X%*%beta,log(fit7$r))


Xbeta = X%*%beta
prob <- pnorm(Xbeta)

dprob<- dnorm(Xbeta)
r<-    c(prob*(1 - prob)/dprob^2)
r[r<1]<- 1





# acf(fit6$beta)
ts.plot(fit6$beta[,1])

mean(fit6$r)

hist(X%*%beta, breaks=100)

mean(X%*%beta)
sd(X%*%beta)^2


ts.plot(fit6$beta[,1])



library(latex2exp)


setwd("~/work/Dropbox/scalable_da/ms/")

pdf("./probit_cda_r.pdf",5,5)
plot(X%*%beta,fit6$r, xlab=TeX('Posterior mean $x_i^T\\beta$'),ylab=TeX('Adapted $r_i'))
dev.off()

pdf("./probit_cda_b.pdf",5,5)
plot(  (sqrt(fit6$r)-1)*X%*%colMeans(fit6$beta),fit6$b, xlab=TeX('Posterior mean $( \\sqrt{r_i}-1 ) x_i^T\\beta$'),ylab=TeX('Adapted $b_i'),xlim=c(-80,0),ylim=c(-80,0))
dev.off()


source("probitNaive.r")
fitNaive <-probitCDA(y,X, burnin = 100, run = 500, fixR = F, MH=T,r_ini = 10,delta = 5)
ts.plot((fitNaive$beta[,1]))
acf((fitNaive$beta[,1]),lag.max = 40)
