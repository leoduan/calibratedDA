# .rs.restartR()

require("ImbalancedPG")
setwd("~/git/ImbalancedPG/test/")

source("maxpoint_data.r")

fit1<- ImbalancedPG::poisson_reg(y , X,b,B, r0ini =  10,c = 1.1,burnin = 2000,run = 2000,fixed_R = T)


fit1$r<- NULL
fit1$w<- NULL


save(fit1,file="maxpoint_fit1.rda")
