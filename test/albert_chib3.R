require(msm)

N <- 1E4
y=0
# while(sum(y)==0){


X0 <- 1#rnorm(N, 1, 1)
X1 <- rnorm(N, 1, 1)
X <- cbind(X0, X1)
beta <- c(-3, 0)

Xbeta<- X%*%beta
theta <- pnorm(X %*% beta)

y <- as.numeric(runif(N) < theta)
# }

sum(y)
vecInf <- rep(Inf, N)

X2inv <- solve(t(X) %*% X)
cholX2inv <- t(chol(X2inv))
lb <- -vecInf
ub <- vecInf


beta <- rnorm(2)
# C <- Xbeta

r0 <- 15
r <- rep(r0,N)

# r<-1

trace_beta <- numeric()

Z<- y
C <- rep(0,N)# (1-r) * rep(-3.578, N)

for (i in 1:300) {
  
  Xbeta  <- X%*%beta
  theta <- Xbeta
  
  C <-  (1-r) * Xbeta

  r<- rep(r0,N)
  r[((Z <= C) & (y==1))] <- 1
  r[((Z > C) & (y==0))] <- 1
  
  C <-  (1-r) * Xbeta
  
  
  lb[y == 1] <- C[y == 1]
  ub[y == 0] <- C[y == 0]
  
  Z <- rtnorm(N, theta, sd = r, lb, ub)
  X2<- t(X)%*% (1/r/r*X)
  X2inv<- solve(X2)
  cholX2inv <- t(chol(X2inv))
  m <-  X2inv %*% (t(X/r/r) %*% (Z))
  beta <- cholX2inv %*% rnorm(2)  + m
  
  
  
  # r <- rep(1,N)

  if (i > 100){
    trace_beta <- rbind(trace_beta, t(beta))
  }
  
  print(i)
}

ts.plot(trace_beta)
# ts.plot(trace_beta[,2])

acf(trace_beta[,1], lag.max = 40)
# acf(trace_beta[,2])

# hist(trace_beta[,1])
