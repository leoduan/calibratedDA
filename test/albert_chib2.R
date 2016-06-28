require(msm)

N <- 1E3
y=0
while(sum(y)==0){


X0 <- 1#rnorm(N, 1, 1)
X1 <- rnorm(N, 1, 1)
X <- cbind(X0, X1)
beta <- c(-3.8, 0)

Xbeta<- X%*%beta
theta <- pnorm(X %*% beta)

y <- as.numeric(runif(N) < theta)
}

vecInf <- rep(Inf, N)

X2inv <- solve(t(X) %*% X)
cholX2inv <- t(chol(X2inv))
lb <- -vecInf
ub <- vecInf


# beta <- rnorm(2)
# C <- Xbeta

r <- 1.1

# r<-1

trace_beta <- numeric()


C <- (1-r) * rep(-3.578, N)

for (i in 1:300) {
  
  Xbeta  <- X%*%beta
  theta <- Xbeta
  
  C <- (1-r) * rep(Xbeta, N)
  
  
  lb[y == 1] <- C[y == 1]
  ub[y == 0] <- C[y == 0]
  
  Z <- rtnorm(N, theta, sd = r, lb, ub)
  m <- X2inv %*% (t(X) %*% (Z)/r)
  beta <- cholX2inv %*% rnorm(2) * r  + m
  
  # r <- rep(1,N)

  if (i > 100){
    trace_beta <- rbind(trace_beta, t(beta))
  }
  
  print(i)
}

ts.plot(trace_beta[,1])
ts.plot(trace_beta[,2])

ts.plot(trace_beta[,1]/trace_beta[,2])

acf(trace_beta[,1])
acf(trace_beta[,2])

# hist(trace_beta[,1])

xbeta<-X%*%colMeans(trace_beta)

1-pnorm(3.578, xbeta, r)
1-pnorm(0,-3.8,1)

1-pnorm(0, xbeta, 1)

hist(pnorm(xbeta))
hist(pnorm(xbeta, C, r))

