# .rs.restartR()
require("ImbalancedPG")
setwd("~/git/ImbalancedPG/test/")

N <- 1E4

X0 <- 1#rnorm(N, 1, 1)
X1 <- rnorm(N, 1, 1)
X <- cbind(X0, X1)
beta <- c(-4, 1)

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
fit<- ImbalancedPG::probit_reg_simple2(y,X,b,B, r0= 100, burnin = 3000, run = 1000)
# fit2<- ImbalancedPG::probit_reg_simple(y,X,b,B, r0= 200, burnin = 1000, run = 1000)
acf(fit$beta,lag.max = 40)
# acf(fit2$beta,lag.max = 40)

colMeans(fit$beta)
apply(fit$beta, 2, sd)
ts.plot(fit$beta)

plot(X%*%beta,fit$r)

hist(fit$beta[,2])

