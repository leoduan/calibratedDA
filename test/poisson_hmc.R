setwd("~/git/ImbalancedPG/test/")

library(rstan) # observe startup messages
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("maxpoint_data.r")

#p<-3
#N<-1000

X<- as.matrix(X[1:N,1:p])
y<-y[1:N]
data = list('X' = X,
            'y' = y,
            'N' = N,
            'p' = p
)

X2<- t(X)%*%X
diag(X2)<- diag(X2)+1E-3
b0<- solve(X2, t(X)%*%log(y+1))


fit <- stan(file = 'poisson.stan', data = data, iter = 2000, chains = 1)

save(fit,file="stan_fit_poisson.Rda")
