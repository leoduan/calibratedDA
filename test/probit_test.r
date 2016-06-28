require("ImbalancedPG")

N <- 1E4

X0 <- 1#rnorm(N, 1, 1)
X1 <- rnorm(N, 1, 1)
X <- cbind(X0, X1)
beta <- c(-5, 1)

theta <- pnorm(X %*% beta)

y <- as.numeric(runif(N) < theta)
sum(y)

b<- rep(0,2)
B<- diag(1000,2)

fit<- ImbalancedPG::probit_reg_simple(y,X,b,B, r0= N/10, burnin = 1000, run = 500)

acf(fit$beta[, 1], lag.max = 40)
acf(fit$beta[, 2], lag.max = 40)

ts.plot(fit$beta[, 1])
ts.plot(fit$beta[, 2])

# acf(trace_beta[,2],lag.max = 40)

# ts.plot(trace_beta[, 2])
# ts.plot(trace_beta[,2])

hist(X%*%beta)
