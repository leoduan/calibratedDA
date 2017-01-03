setwd("~/git/ImbalancedPG/poisson/")

library(rstan) # observe startup messages
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# 

setwd("~/git/ImbalancedPG/poisson/")
source("maxpoint_data.r")

source("poisson.r")

X<-X

p<- ncol(X)

N<- nrow(X)

X<- as.matrix(X[1:N,1:p])
y<-y[1:N]



# n<- 1E4
# X<- cbind(1,rnorm(n))
# 
# beta<- c(-3,1)
# N<-n
# 
# y<-rpois(n,exp(X%*%beta + rnorm(n,sd = 1)))
# p<- 1

data = list('X' = X,
            'y' = y,
            'N' = N,
            'p' = p
)

fit <- stan(file = 'poissonRanEff.stan', data = data, iter = 2000, chains = 1, init="0")

save(fit,file="stan_fit_poisson.Rda")

# 
load(file="stan_fit_poisson.Rda")
# # # 
acf(fit@sim$samples[[1]]$`beta[1]`[1001:2000],lag.max = 40)
acf(fit@sim$samples[[1]]$`beta[100]`[1001:2000])
# acf(fit@sim$samples[[1]]$`nu`[1001:2000])
# acf(fit@sim$samples[[1]]$`eta0`[1001:2000])
# 
# # 
# # # 
# ts.plot(fit@sim$samples[[1]]$`beta[1]`)
# ts.plot(fit@sim$samples[[1]]$`beta[2]`)
# ts.plot(fit@sim$samples[[1]]$`nu`)
# ts.plot(fit@sim$samples[[1]]$`eta0`[1001:2000])


# 
# mean(fit@sim$samples[[1]]$sigma2)
# mean(fit@sim$samples[[1]]$`eta[1]`)
# 

# mean(1/rgamma(1000,2,1))
