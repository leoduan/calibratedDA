require("ImbalancedPG")

N <- 1E4

X0 <- 1#rnorm(N, 1, 1)
X1 <- rnorm(N, 1, 1)
X <- cbind(X0, X1)
beta <- c(-3, 0)

# hist(X %*% beta)
theta <- pnorm(X %*% beta)

y <- as.numeric(runif(N) < theta)
sum(y)

b<- rep(0,2)
B<- diag(1000,2)

fit<- ImbalancedPG::probit_reg_simple(y,X,b,B, r0= 40, burnin = 100, run = 200)

acf(fit$beta[, 1], lag.max = 40)
acf(fit$beta[, 2], lag.max = 40)

ts.plot(fit$beta)

# acf(trace_beta[,2],lag.max = 40)

ts.plot(fit$beta[, 1])
ts.plot(fit$beta[, 2])

sd(fit$beta[, 1])


plot(pnorm(X%*%colMeans(fit$beta)), pnorm(X%*%(beta)), xlim=c(0,0.1),ylim=c(0,0.1))
