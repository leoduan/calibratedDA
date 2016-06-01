library(devtools)

setwd("~/git/")
build('ImbalancedPG/')
install.packages("ImbalancedPG_1.0.tar.gz", repos = NULL, type = "source")

unload("ImbalancedPG")
require("ImbalancedPG")
require("ggplot2")
require("reshape")

#logit regression

N<- 1E5

X1<- rnorm(N, 0, 0.5)
# X1<- (X1 - mean(X1))/sd(X1)

X<- cbind(1,X1)

beta0<- as.vector(c(-8, 1))

theta<- (X%*%beta0)
p<- 1/(1+ exp(-theta))

y0<- as.numeric(runif(N)< p)
sum(y0)

b<- rep(0,2)
B<- diag(1000,2)

system.time(
test<- ImbalancedPG::logit_reg(y0, X, b, B, mc_draws=5, r0_ratio = 10 ,r1= 2,burnin = 1000,run = 1000, downsampling = TRUE)
)

#trace of the percentage of y=0's ignored
mean(test$filter_count)

plot_acf<- function(x){
  bacf <- acf(x, plot = FALSE,lag.max = 40)
  bacfdf <- with(bacf, data.frame(lag, acf))
  
  q <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0)) +  theme_bw() 
  q
}



pdf("acf_reg_imbal_beta0.pdf",4,4)
plot_acf(test$beta[,1])
dev.off()

pdf("acf_reg_imbal_beta1.pdf",4,4)
plot_acf(test$beta[,2])
dev.off()

hist(test$beta[,1])
hist(test$beta[,2])
colMeans(test$beta)
sd(test$beta[,1])
sd(test$beta[,2])
acf(test$beta[,1])
acf(test$beta[,2])



require(BayesLogit)
fit2 <- logit(y0,X,samp = 2000,burn = 1000)
ts.plot(fit2$beta[1500:2000,1])
ts.plot(fit2$beta[1500:2000,2])

colMeans(fit2$beta[1500:2000,])
sd(fit2$beta[1000:2000,1])
sd(fit2$beta[1000:2000,2])

pdf("acf_reg_pg_beta0.pdf",4,4)
plot_acf(fit2$beta[,1])
dev.off()

pdf("acf_reg_pg_beta1.pdf",4,4)
plot_acf(fit2$beta[,2])
dev.off()

hist(fit2$beta[,1])
hist(fit2$beta[,2])

ts.plot(fit2$beta[,2])
