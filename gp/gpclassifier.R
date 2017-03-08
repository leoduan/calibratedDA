require("truncnorm")
require("reshape")
require("ggplot2")

n<- 100
w<- c(0:(n-1))/n*2*pi
w_folded<- (pi-abs(w - pi))

rho<- 50
G<- 100*exp(- 2* rho^2* w_folded^2)
Q<- mvfft(diag(1,n))/sqrt(n)

sum(G>1E-5)/2

M<- G[1:10]
M<- c(M, M[-1])
S1<- Re(Q[,1:10])
S2<- Im(Q[,2:10])
S<- cbind(S1,S2)
S<- S%*%diag(sqrt(M))

m<- 19
Yc<- rnorm(m)

priorMean<- -3.2
Mu0<- S%*% Yc + priorMean

if(spectral){
  Yc<- complex(n, real = rnorm(n),imaginary =  rnorm(n))
  Mu0<- Re(fft(G^0.5 *Yc)/sqrt(n)) -5
}
# plot(Mu0)

prob<- pnorm(Mu0)

y<- (runif(n)<prob)*1
y<- rep(0,n)
# y[5:6]<-1
# y=rep(0,n)
y[10]<-1
sum(y)/n

plot(y)

lines(prob)
plot(prob)


# SZ<- t(S)%*%Z
# r<- rep(1,n)
# r<- runif(n)
# ZrZ<- SZ/r


#Woodbury
# mat1<- solve(diag(r,n) + S%*%t(S))
# mat2=  diag(1/r) - diag(1/r) %*% S %*%
# solve ( diag(1,m) + t(S)%*% diag(1/r) %*% S  ) %*% t(S) %*% diag(1/r)



# SS<- S%*%t(S)
# SSr<- SS +diag(r)
# temp1<-SS%*%solve(SSr,Z)
# temp2<- S%*%mMu0
# temp1<-SS-SS%*%solve(SSr,t(SS))
# temp2<-S %*% vMu0%*%t(S)
# max(abs(temp1-temp2))

tlog<- function(x){
  x[x<1E-12]= 1E-12
  log(x)
}


runGP<-function(c0, fixR=F, MH=T, spectral=FALSE){
  trace_mu<-numeric()
  
  lb<- rep(-Inf,n)
  ub<- rep(Inf,n)
  lb[y==1]<-0
  ub[y==0]<-0
  
  Mu<- rnorm(n)
  r=rep(1,n)
  prob<- pnorm(Mu)
  loglik <- tlog(1-prob) * (y==0) + tlog(prob) * (y==1)
  b<- rep(0, n)
  accept<-0
  
  for(i in 1:2000){
    
    
    Z<- rtruncnorm(n,a= lb,b= ub,mean = Mu+b, sd = sqrt(r))
    
    if(spectral){
      m<- Re(fft( fft(Z-b)/sqrt(n) / (G+r)*G, inverse = T))/sqrt(n)
      qv<- G*r/(G+ r)
      Yc<- complex(n, real = rnorm(n),imaginary =  rnorm(n))
      new_Mu<- Re(fft(sqrt(qv)*Yc,inverse = T))/sqrt(n) + m
    }else{
      
      SrS<- t(S)%*% (S/r)
      SrZ<- t(S)%*% ((Z-b)/r)
      ISrS = SrS+diag(1,m)
      ISrsInv<- solve(ISrS)
      mMu0<- (SrZ - SrS%*% ISrsInv %*% SrZ)
      vMu0<- diag(1,m)-SrS + SrS %*% ISrsInv %*% t(SrS)
      new_Mu<- S%*%(t(chol(vMu0))%*% rnorm(m) + mMu0)
    }
    new_prob <- pnorm(new_Mu)
    new_loglik<-  tlog(1-new_prob) * (y==0) + tlog(new_prob) * (y==1)
    
    q_prob <- pnorm((Mu+b)/sqrt(r))
    q_loglik<-  tlog(1-q_prob) * (y==0) + tlog(q_prob) * (y==1)
    
    new_q_prob <- pnorm((new_Mu+b)/sqrt(r))
    new_q_loglik<-  tlog(1-new_q_prob) * (y==0) + tlog(new_q_prob) * (y==1)
    
    # inidividual likelihood ratio
    alpha<- (new_loglik + q_loglik  - loglik - new_q_loglik  )
    
    #metropolis-hastings
    if(runif(1)< exp(sum(alpha))|| !MH)
    { 
      Mu<- new_Mu
      prob<- new_prob
      loglik<- new_loglik
      accept<- accept +1
      # if(tuning) {
      #   tune_accept<- tune_accept+1
      # }
    }
    
    
    if(i<500){
      if(!fixR){
        dprob<- dnorm(Mu)
        r<-    (c(prob*(1 - prob)/dprob^2)) /c0
        # r=rep(c0,n)
        r[r<1]<-1
        b<- (sqrt(r)-1) * Mu
      }
    }
    print(i)
    print(accept/(i+1))
    
    if(i>500)
      trace_mu<- rbind(trace_mu, t(Mu))
  }
  list("mu"= trace_mu, "accept"=accept/2000, "r"=r)
}


acfBoxPlot<- function(x){
  
  ACFfit<-apply(x, MARGIN = 2, function(x){c(c(acf(x,plot = F,lag.max = 40))$acf[,,1])})
  
  CDAacf<- t(ACFfit)
  colnames(CDAacf)<- c(0:40)
  
  dfCDAacf <- melt(CDAacf)
  colnames(dfCDAacf)<- c("Variable","Lag","ACF")
  dfCDAacf$Lag<- as.factor(dfCDAacf$Lag)
  
  p<- ggplot(data = dfCDAacf, aes(x=Lag, y=ACF)) + geom_boxplot(outlier.shape = NA)+
    scale_y_continuous(limits = c(-0.2,1))+ scale_x_discrete(breaks = c(0:8)*5)
  p
}

acfMean<- function(x){
  
  ACFfit<-apply(x, MARGIN = 2, function(x){c(c(acf(x,plot = F,lag.max = 40))$acf[,,1])})
  
  rowMeans(ACFfit)
}

fit1<- runGP(10,fixR = F,MH = T,spectral = F)
fit1$accept
fit1$r

fit2<- runGP(1,fixR = T,MH=F,spectral = F)
fit2$r


fit3<- runGP(1,fixR = F,MH = T,spectral = F)
fit3$accept

par(mfrow=c(3,1))
acf(fit1$mu[,1],lag.max = 40)
acf(fit2$mu[,1],lag.max = 40)
acf(fit3$mu[,1],lag.max = 40)

# ts.plot(fit1$mu[,8])
# ts.plot(fit2$mu[,8])

acfBoxPlot(fit2$mu)
acfBoxPlot(fit3$mu)

par(mfrow=c(1,1))
plot(acfMean(fit1$mu),type="l",ylim=c(-0.2,1))
lines(acfMean(fit2$mu),col="red")
lines(acfMean(fit3$mu),col="blue")


# plot(y)
# plot(pnorm(colMeans(fit$mu)),col="red")
# plot(pnorm(colMeans(fit2$mu)))
# lines(pnorm(Mu0),col="blue")
