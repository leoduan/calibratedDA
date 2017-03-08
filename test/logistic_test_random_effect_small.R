# .rs.restartR()
# require("scalableDA")
require("BayesLogit")

setwd("~/git/ImbalancedPG/test/")

N<- 1E3

X0<- 1#rnorm(N, 1, 1)
X1<- rnorm(N, 1, 1)
X<- cbind(X0,X1)
beta<- c(-5,1)
p<- ncol(X)


Xbeta<- X%*%beta
sigma<- 3


XbetaEps<- rnorm(N,Xbeta,sd = sqrt(sigma))
theta<-  exp(XbetaEps)/(1+exp(XbetaEps))

a<-numeric()
for(i in 1:10){
  y<- as.numeric(runif(N)<theta)
  a<- c(a, sum(y))
}

a
# mean(a)/N
# sd(a)/N

beta<- c(-5,1)

Xbeta <- X%*%beta
trace_beta<- numeric()
trace_sigma<- numeric()
eps<- rnorm(N,0,sqrt(1))
w<- rep(1,N)
r<- pmin( exp(Xbeta +eps) *100,1)
logr<- log(r)

theta<- XbetaEps

r<- 1
logr<- log(r)
k = y - r / 2.0 + w * logr;

max_theta<- rep(-Inf, N)

for(i in 1:1000){
  
  if(i>100){
    max_theta <-  pmax(theta,max_theta)
    r<- pmin( exp(max_theta) *10,1)
    logr<- log(r)
  }
  
  #update w
  
  alpha1<- r
  alpha2<- abs( theta - logr)
  w <- rpg(N,alpha1,alpha2)
  fail<- is.na(w)
  if(length(fail)>1){
    w[fail]<- alpha1[fail]/2/ alpha2[fail] * tanh(alpha2[fail])
  }
  
  #update theta & eps
  var = 1 / (1 / sigma + w);
  m = var * (y- r/2 +  w*(logr- Xbeta))
  eps = sqrt(var) * rnorm(N) + m
  # m = var * (w* logr + Xbeta/sigma + y - r/2)
  # theta = sqrt(var) * rnorm(N) + m
  # eps = theta - Xbeta
  theta<- Xbeta+ eps
  
  #update sigma
  # a= N/2;
  # b = sum((theta-Xbeta)^2)/2;
  # sigma = 1/rgamma(1,a,rate=b)  ;
  
  
  #update beta
  
  # V<- solve(t(X*w)%*%X +diag(1E-3,nrow = 2))
  # k<- y- r/2 + w* (logr- eps)
  # M = V %*% (t(X)%*% k);
  
  v<- (w - w^2/(w+1/sigma))
  V<-  solve(t(X*v) %*%X +diag(1E-3,nrow = 2))
  M<-  V%*% (t(X)%*%( (y- r/2+ w*logr )*(1- w/(w+1/sigma))))
  
  # V<- solve(t(X)%*%X) * sigma
  # M<- V %*% (t(X)%*% (theta) /sigma );
  
  
  beta = t(chol(V)) %*% rnorm(p) + M;
  Xbeta = X%*% beta;
  
  
  print(i)
  print(beta)
  
  if(i>100){
    trace_beta<- rbind(trace_beta, t(beta))
    trace_sigma<- c(trace_sigma, sigma)
  }
}

ts.plot(trace_beta)

ts.plot(trace_sigma)
ts.plot(trace_sigma[1:100])

acf(trace_beta)
colMeans(trace_beta)
acf(trace_beta)




fit<- ImbalancedPG::logit_reg_random_effect(y , X,r0 = 10,burnin = 20,run = 20)

fit$beta
colMeans(fit2$beta)
apply(fit2$beta, 2, sd)


colMeans(fit$beta)
apply(fit$beta, 2, sd)




require("ggplot2")
require("reshape")

a<- fit$r

df<- melt(a[,1:20])

pdf("./r_adaptation.pdf",3,3)
ggplot(data=df, aes(x=X1, y=value,col= factor(X2)))+ geom_line(size=0.5)+  scale_colour_manual(values = rep("black",30))+ylab("r")+xlab("Iteration")+  theme(legend.position="none")
dev.off()



df<- data.frame("z"=c((fit$w[1000,]), (fit2$w[1000,])),"Method"=rep(c("CDA","DA"),each=N))

pdf("./logit_z.pdf",4,3)
ggplot(df, aes(x=z, fill=Method, col=Method)) + geom_histogram(alpha=0.2, position="identity", bins=50)+xlim(-1,0.3)+ylim(0,15000)+ scale_fill_grey(start = 0,end = 0.5)+ scale_color_grey(start = 0,end = 0.5)
dev.off()





pdf("./logit_15.pdf",4,3)
ggplot(data=df, aes(x=Lag, y=ACF,linetype=Method))+ geom_line(size=.75)
dev.off()

pdf("./logit_19.pdf",4,3)
ggplot(data=df, aes(x=Lag, y=ACF,linetype=Method))+ geom_line(size=.75)
dev.off()





fit$beta
acf(fit$beta[,1])
acf(fit$beta[,2])

# hist(fit$beta[,1])
# hist(fit$beta[,2])

sd(fit$beta[,1])
sd(fit$beta[,2])

ts.plot(fit$beta)


require("BayesLogit")
fit2<- logit(y,X) 

acf(fit2$beta[,1])
acf(fit2$beta[,2])

sd(fit2$beta[,1])
sd(fit2$beta[,2])

ts.plot(fit2$beta)

hist(fit$w,xlim=c(0,1))
hist(fit2$w,xlim=c(0,1))

sd(fit$w[,3])
sd(fit2$w[,3])


Xbeta <-fit2$beta %*% t(X)
Xbeta2 <-fit2$beta %*% t(X)

median(1/fit$w[,3])
median(1/fit2$w[,3])

sd(1/fit2$w[,1])
fit$beta

#binomial
fit3<- logit(1,n = 1, X=1)
acf(fit3$w)

