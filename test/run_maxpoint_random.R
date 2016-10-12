# .rs.restartR()

require("scalableDA")
setwd("~/git/ImbalancedPG/test/")

source("maxpoint_data.r")



fit1<- poisson_reg_random_effect(y , X, tau =  10,c = 1,burnin = 2000,run = 2000,da_ver = 1)
fit2<- poisson_reg_random_effect(y , X, tau =  100,c = 1,burnin = 2000,run = 2000,da_ver = 1)



save(fit1,file="maxpoint_fit1_random.rda")
save(fit2,file="maxpoint_fit2_random.rda")
