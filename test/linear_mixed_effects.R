
n <- 1E4

X0 <- 1#rnorm(N, 1, 1)
X1 <- rnorm(n, 1, 1)
X <- cbind(X0, X1)
beta <- c(-1, 1)


# X<- matrix(rep(1,n))
# beta<- 1

p<- length(beta)
# hist(X %*% beta)

sigma1<-  1^2
sigma2<-   0.1^2


Z <- rnorm(n,X %*% beta,sqrt(sigma1))

y<- rnorm(n,Z,sd = sqrt(sigma2))


X2<- t(X)%*%X
trace_beta<- numeric()

for(i in 1:1000){
  Xbeta<- X%*%beta
  var<- 1/(1/sigma1 + 1/sigma2)
  m<- var* (Xbeta/sigma1 + y/sigma2)
  Z<-  rnorm(n,m, sqrt(var))
  eps<- Z - Xbeta
  
  m<- solve(X2,t(X)%*% Z )
  var<- solve(X2)*sigma1
  
  # m<- solve(X2,t(X)%*% (y-eps) )
  # var<- solve(X2)*sigma2
  
  beta <- t(chol(var))%*%rnorm(2)+m
  trace_beta<- rbind(trace_beta, t(beta))
}

acf(trace_beta[,1])
acf(trace_beta[,2])
ts.plot(trace_beta)



#######
#logit#
#######

require("BayesLogit")

n <- 1E4

X0 <- 1#rnorm(N, 1, 1)
X1 <- rnorm(n, 1, 1)
X <- cbind(X0, X1)


X2<- t(X)%*%X
X2inv<- solve(X2)

# X<- matrix(rep(1,n))
# beta<- 1

p<- length(beta)
# hist(X %*% beta)

sigma1<-  1^2

beta <- c(-10, 0.1)
Z <- rnorm(n,X %*% beta,sqrt(sigma1))

y<- (runif(n)< 1/(1+exp(-Z)))*1

table(y)

ver=1

logit_random_effect<- function(ver=1){
  k<- y - 0.5
  beta<- c(-1,1)
  eps<- Z - X%*%beta
  trace_beta<- numeric()
  trace_sigma1<- numeric()
  
  for(i in 1:1000){
    
    if(ver==1){
      m<- X2inv %*% (t(X)%*% Z )
      var<- X2inv*sigma1
    }
    if(ver==2){
      var <- solve( t(X) %*% (w* X))
      m<- var %*% (t(X)%*% (k - w* eps))
    }
    beta <- t(chol(var))%*%rnorm(2)+m
    beta
    
    var<- 1/(w + 1/ sigma1)
    Xbeta<- X%*%beta
    m<-   var*(k + Xbeta/sigma1)
    Z<- sqrt(var)*rnorm(n)+m
    eps<-  Z - Xbeta
    
    w <- rpg(n, 1, Z)
    
    sigma1<- rgamma(1, n/2, rate= sum(eps^2)/2)
    
    
    trace_beta<- rbind(trace_beta, t(beta))
    trace_sigma1<- c(trace_sigma1, sigma1)
    
  }
  list(
    "beta"= trace_beta,
    "sigma1"=trace_sigma1
  )
}

fit1<- logit_random_effect(ver=1)
fit2<- logit_random_effect(ver=2)

ts.plot(fit1$sigma1)
ts.plot(fit2$sigma1)



#######
#probit#
#######

n <- 1E4

X0 <- 1#rnorm(N, 1, 1)
X1 <- rnorm(n, 1, 1)
X <- cbind(X0, X1)


X2<- t(X)%*%X
X2inv<- solve(X2)

# X<- matrix(rep(1,n))
# beta<- 1

p<- length(beta)
# hist(X %*% beta)

sigma1<-  5^2

beta <- c(-5, 1)
xbetaEps <- rnorm(n,X %*% beta,sqrt(sigma1))

y<- (rnorm(n,xbetaEps,1)>0)*1

table(y)

ver=1



probit_random_effect<- function(ver=1){
  beta<- c(-1,1)
  eps<- Z - X%*%beta
  trace_beta<- numeric()
  trace_sigma1<- numeric()
  Xbeta<- X%*%beta
  
  
  lb<- rep(-Inf,n)
  ub<- rep(Inf,n)
  lb[y==1]<- 0
  ub[y==0]<- 0
  
  for(i in 1:1000){
    
    Z<- rtruncnorm(n,lb,ub,xbetaEps,1)
    
    
    Xbeta<- X%*%beta
    var<- 1/(1/sigma1 + 1/1)
    m<- var* (Xbeta/sigma1 + Z/1)
    xbetaEps <-  rnorm(n,m, sqrt(var))
    eps<- xbetaEps - Xbeta
    
    if(ver==1){
      m<- solve(X2,t(X)%*% xbetaEps )
      var<- solve(X2)*sigma1
    }
    if(ver==2){
      m<- solve(X2,t(X)%*% (Z-eps) )
      var<- solve(X2)*1
    }
    
    
    beta <- t(chol(var))%*%rnorm(2)+m
    beta
    
    sigma1<- 1/rgamma(1, n/2, rate= sum(eps^2)/2)
    
    if(i>100){
      trace_beta<- rbind(trace_beta, t(beta))
      trace_sigma1<- c(trace_sigma1, sigma1)
    }
    
  }
  list(
    "beta"= trace_beta,
    "sigma1"=trace_sigma1
  )
}

fit1<- probit_random_effect(ver=1)
fit2<- probit_random_effect(ver=2)

ts.plot(fit1$sigma1)
ts.plot(fit2$sigma1)

# acf(fit1$sigma1)
# acf(fit2$sigma1)

ts.plot(fit1$beta)
ts.plot(fit2$beta)

acf(fit1$beta)
acf(fit2$beta)


