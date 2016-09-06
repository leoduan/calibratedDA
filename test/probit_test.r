.rs.restartR()
require("ImbalancedPG")

N <- 1E4

X0 <- 1#rnorm(N, 1, 1)
X1 <- rnorm(N, 1, 1)
X <- cbind(X0, X1)
beta <- c(-2, 0)

# hist(X %*% beta)
theta <- pnorm(X %*% beta)

y <- as.numeric(runif(N) < theta)
sum(y)


b<- rep(0,2)
B<- diag(1000,2)


fit<- ImbalancedPG::probit_reg_px(y,X,b,B, r0= 1, burnin = 1000, run = 500,nu0 = 0.1)

ts.plot(fit$beta)

hist(fit$trace_r0)

# acf(fit$beta[fit$trace_r0>1,1])
# acf(fit$beta[fit$trace_r0>1,2])

acf(fit$beta[, 1], lag.max = 40)
acf(fit$beta[, 2], lag.max = 40)

sd(fit$beta[,1])

fit2<- ImbalancedPG::probit_reg_simple(y,X,b,B, r0= 2, burnin = 1000, run = 1000)

ts.plot(fit2$beta)

acf(fit2$beta[, 1], lag.max = 40)
acf(fit2$beta[, 2], lag.max = 40)

sd(fit2$beta[,1])

table(fit2$r)

