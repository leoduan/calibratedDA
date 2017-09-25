setwd("~/git/calibratedDA/R-demo/")
require("ggplot2")

#generate data
N <- 1E4


#intercept only

y<- rep(0,N)
y[1]<-1

X<- matrix(1,N)
# sum(y)

source("probit.r")

#run cDA
fit10<-probitCDA(y,X, burnin = 1000, run = 1000,r_ini = 10,fixR = T)
fit100<-probitCDA(y,X, burnin = 200, run = 1000,r_ini = 100,fixR = T)
fit1000<-probitCDA(y,X, burnin = 200, run = 1000,r_ini = 1000,fixR = T)
fit5000<-probitCDA(y,X, burnin = 200, run = 1000,r_ini = 5000,fixR = T)
fit10000<-probitCDA(y,X, burnin = 200, run = 1000,r_ini = 10000,fixR = T)




fit1 <-probitCDA(y,X, r_ini = 1,burnin = 500, run = 1000, fixR = T)
fit2 <-probitCDA(y,X, r_ini = 10,burnin = 500, run = 1000, fixR = T)
fit3 <-probitCDA(y,X, r_ini = 100,burnin = 500, run = 1000, fixR = T)
fit4 <-probitCDA(y,X, r_ini = 1000,burnin = 500, run = 1000, fixR = T)
fit5 <-probitCDA(y,X, r_ini = 5000,burnin = 500, run = 1000, fixR = T)



posterior<- data.frame(theta= c( fit1$proposal[,1]
                                 ,fit2$proposal[,1]
                                 ,fit3$proposal[,1]
                                 ,fit4$proposal[,1]
                                 ,fit5$proposal[,1]
), r= as.factor(rep(c(1,10,100,1000,5000), each=1000)))



pdf("./density_probit.pdf",5,3)
ggplot(posterior, aes(theta, linetype = r))+  geom_density(alpha = 0.2,bw = 0.2,lwd=2) + xlim(-8,0)
dev.off()




acf_all <- apply(cbind(fit1$beta,fit2$beta,fit3$beta,fit4$beta,fit5$beta),2,function(x){ acf(x, lag.max = 40, plot=F)$acf})

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


fit6 <-probitCDA(y,X, burnin = 500, run = 500, fixR = F, MH=T)
fit6$accept_rate

ts.plot(fit6$beta)
acf(fit6$beta)
colMeans(fit6$beta)

fit7 <-probitCDA(y,X, burnin = 500, run = 500, fixR = F,MH = F)
ts.plot(fit7$beta)
acf(fit7$beta)



colMeans(fit7$beta^2)
sum(exp(fit7$importance)* fit7$beta[,1]^2) /sum(exp(fit7$importance))
sum(exp(fit7$importance)* fit7$beta[,2]^2) /sum(exp(fit7$importance))
sum(exp(fit7$importance)* fit7$beta[,3]^2) /sum(exp(fit7$importance))

# fit6$r[fit6$r>10000]<-NA

plot(X%*%beta,log(fit6$r))

acf(fit6$beta)
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

