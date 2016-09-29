require("ImbalancedPG")

n<- c(1E7,1E8,1E9,1E10,1E12,1E15)

# theta<- -25
# p<- 1/(1+exp(-theta))

y<- rep(10,6)#rbinom(2,N,p)

fit<- ImbalancedPG::binomial_simple(y,n, r0 = 10,burnin = 500,run = 500)

fit$theta

getACF<- function(x){
  c(acf(x, lag.max = 40,plot = F)$acf)
}


df<- data.frame("ACF"=c(apply(fit$theta, 2, getACF)),"p"= as.factor( rep(10/n,each=41)),"Lag"=rep(c(0:40),2))

pdf("./logit_big_n_test.pdf",4,3)
ggplot(data=df, aes(x=Lag, y=ACF,linetype=p,col=p))+ geom_line(size=.75)
dev.off()



acf(fit$theta[,1])
acf(fit$theta[,2])
acf(fit$theta[,6])


ts.plot(fit$theta)
pM<- colMeans(fit$theta)

pM
