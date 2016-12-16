
require('R.matlab')
setwd("~/git/ImbalancedPG/poisson/")
mp<-readMat("./convdat.mat")

# 
# dim(mp$Ntr)
# dim(mp$Ytr)
# 
# dim(mp$Nts)
# dim(mp$Yts)
# 
nonzero<-  c(mp$Ntr!=0) & c(mp$Nts!=0)
# 
Ytr<- mp$Ytr[nonzero,]
Ntr<- mp$Ntr[nonzero]
# 
y<- Ytr[1:100,46]
N<- Ntr[1:100]
# X<- matrix(rep(1, length(N)))
# N<-N[1:4000]
# y<-y[1:4000]
# X<-as.matrix(X[1:4000,])
# 
# source("../R-demo/logit.r")
# fit<- logitCDA(y,X, N,r_ini= 1,tune = 300,burnin=300, run=500 ,fixR = F,MH = F,c = 100)
# 
# mean(fit$beta)
# mean(fit$beta^2)
# sum(exp(fit$importance)*fit$beta)/sum(exp(fit$importance))
# sum(exp(fit$importance)*fit$beta^2)/sum(exp(fit$importance))
# ts.plot(fit$beta)
# 
# 
# fit2<- logitCDA(y,X, N,r_ini= 1,tune = 200,burnin=200, run=1000, fixR = F,MH = T,c = 10)
# 
# fit2$r
# mean(fit2$beta)
# mean(fit2$beta^2)
# fit2$accept_rate
# 
# acf(fit2$beta,lag.max = 40)
# ts.plot(fit2$beta)

# sRateR<- svd(rateR)
# max(rateR,na.rm = T)
# max(rateS,na.rm = T)

data<- as.matrix(mp$Ytr)
keep<- colSums(data)>0
data<- data[,keep]

# data2<- as.matrix(mp$Ytr)[30001:59792,]
data2<- as.matrix(mp$Yts)
data2<- data2[,keep]


colSums(data)

cor_Y<- cor(data)


highCorIdx<- c(1:101)[rowSums(cor_Y==max(cor_Y[cor_Y<0.9],na.rm = T),na.rm = T)==1]

# cor_Y[highCorIdx[2],]

pick<- highCorIdx[2]
# pick<- 1

# pick<- 10

N<- nrow(data)

# N<- 10000

o<- order(cor_Y[pick,],decreasing = T)
feature<-   o<= 100
feature[pick]<-FALSE



X<- data[,feature]
X2<- data2[,feature]

y<- data[,pick]
y2<- data[,pick]

keep <- apply(X, 2, sd)>0
X<- X[,keep]

X2<- X2[,keep]



X<- log(X+1)
X<- cbind(1,X)

X2<- log(X2+1)
X2<- cbind(1,X2)

X<- cbind(X, log(mp$Ntr+1))
X2<- cbind(X2, log(mp$Nts+1))


p<- ncol(X)


B<- diag(1000,p,p)
b<- rep(0,p)

sum(y>0)/N

