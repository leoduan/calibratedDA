p<- 96

hmc_beta<- matrix(0,1000,p)

for(i in 1:p){
  hmc_beta[,i]<-  eval( parse(text= paste("fit@sim$samples[[1]]$`beta[", i,"]`[1001:2000]" ,sep="")))
}

p<- 1000
hmc_eta<- matrix(0,1000,p)



mean(colMeans(abs(hmc_beta)))


mean(fit@sim$samples[[1]]$`beta0`[1001:2000])
quantile(fit@sim$samples[[1]]$`beta0`[1001:2000],c(0.025,0.975))


mean(fit@sim$samples[[1]]$`sigma2`[1001:2000])
quantile(fit@sim$samples[[1]]$`sigma2`[1001:2000],c(0.025,0.975))


essFit2<-effectiveSize(hmc_beta[,1:96])/1000
mean(essFit2)
quantile(essFit2,c(0.025,0.975))
