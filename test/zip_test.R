require("ImbalancedPG")

require("rjags")


setwd("~/git/ImbalancedPG/test/")

source("maxpoint_data.r")                      

N<- 1000

X<- X[1:N,1:p]
y<- y[1:N]

jags <- jags.model('zip.jags',
                   data = list('X' = as.matrix(X),
                               'y' = y,
                               'N' = N,
                               'P' = p
                               ),
                   n.chains = 1,
                   n.adapt = 200)

update(jags, 1000)

jags_zip_fit<- coda.samples(jags,
             c( 'b', 'p'),
             1000)


save(jags_zip_fit, file="jags_zip_fit.Rda")
