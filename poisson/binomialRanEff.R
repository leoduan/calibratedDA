setwd("~/git/ImbalancedPG/poisson/")

library(rstan) # observe startup messages
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


require('R.matlab')
setwd("~/git/ImbalancedPG/poisson/")
mp<-readMat("./convdat.mat")

nonzero<-  c(mp$Ntr!=0) & c(mp$Nts!=0)
Ytr<- mp$Ytr[nonzero,]
Ntr<- mp$Ntr[nonzero]
y<- Ytr[,46]
N<- Ntr
n<- length(N)

data = list(
            'y' = y,
            'n' = n,
            'N' = N
)

fit <- stan(file = 'binomialRanEff.stan', data = data, iter = 2000, chains = 1, init="0")

save(fit,file="stan_fit_binomial.Rda")


load(file="stan_fit_binomial.Rda")
# 
# acf(fit@sim$samples[[1]]$`beta[1]`,lag.max = 40)
# acf(fit@sim$samples[[1]]$`beta[2]`)

p<- n
hmc_eta<- matrix(0,1000,p)

for(i in 1:p){
  hmc_eta[,i]<-  eval( parse(text= paste("fit@sim$samples[[1]]$`eta[", i,"]`[1001:2000]" ,sep="")))
}


require("reshape")
require("ggplot2")

ACFfit3<-apply(hmc_eta, MARGIN = 2, function(x){c(c(acf(x,plot = F,lag.max = 40))$acf[,,1])})
ts.plot(ACFfit3)

CDAacf<- t(ACFfit3)
p<- nrow(CDAacf)

colnames(CDAacf)<-  c(0:40)

dfCDAacf <- melt(CDAacf)
colnames(dfCDAacf)<- c("Variable","Lag","ACF")
dfCDAacf$Lag<- as.factor(dfCDAacf$Lag)

pdf("./binomial_acf_hmc.pdf",5,3)
ggplot(data = dfCDAacf, aes(x=Lag, y=ACF)) + geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(-0.2,1))+ scale_x_discrete(breaks = c(0:8)*5)
dev.off()



# 
# ts.plot(fit@sim$samples[[1]]$`beta[1]`)
# ts.plot(fit@sim$samples[[1]]$`beta[2]`)
# ts.plot(fit@sim$samples[[1]]$`sigma2`)
# 
# mean(fit@sim$samples[[1]]$sigma2)
# mean(fit@sim$samples[[1]]$`eta[1]`)
# 

# mean(1/rgamma(1000,2,1))
