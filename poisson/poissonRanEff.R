setwd("~/git/ImbalancedPG/poisson/")

library(rstan) # observe startup messages
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("maxpoint_data.r")

# N<- 1E3
# X0<- 1#rnorm(N, 1, 1)
# X1<- rnorm(N, 1,1)
# X<- cbind(X0,X1)
# beta<- c(-5,1)
# p<- ncol(X)
# 
# Xbeta<- X%*%beta
# eta0<-  rnorm(N,0,sd= 1)
# theta<-  exp(Xbeta + eta0)
# y<- rpois(N,theta)
# 
# sum(y)
# n<- N


X<- as.matrix(X[1:N,1:p])
y<-y[1:N]
data = list('X' = X,
            'y' = y,
            'N' = N,
            'p' = p
)

fit <- stan(file = 'poissonRanEff.stan', data = data, iter = 2000, chains = 1)

save(fit,file="stan_fit_poisson.Rda")


# load(file="stan_fit_poisson.Rda")
# 
# acf(fit@sim$samples[[1]]$`beta[1]`,lag.max = 40)
# acf(fit@sim$samples[[1]]$`beta[2]`)
# 
# ts.plot(fit@sim$samples[[1]]$`beta[1]`)
# ts.plot(fit@sim$samples[[1]]$`beta[2]`)
# ts.plot(fit@sim$samples[[1]]$`sigma2`)
# 
# mean(fit@sim$samples[[1]]$sigma2)
# mean(fit@sim$samples[[1]]$`eta[1]`)
# 

# mean(1/rgamma(1000,2,1))
