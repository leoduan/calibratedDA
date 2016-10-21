
#generate data
N <- 1E4

X0 <- 1
X1 <- rnorm(N, 1, 1)
X <- cbind(X0, X1)
beta <- c(-5, 1)

theta <- pnorm(X %*% beta)

a<-numeric()
y <- as.numeric(runif(N) < theta)

#run cDA
fit<-probitCDA(y,X, burnin = 1000, run = 1000)

ts.plot(fit$beta)

acf(fit$beta, lag.max = 40)


#posterior mean
beta<- colMeans(fit$beta)



#run the original Albert-Chib
fit_AC<-probitCDA(y,X, r_ini = 1,burnin = 1000, run = 1000)


ts.plot(fit_AC$beta)
acf(fit_AC$beta, lag.max = 40)
