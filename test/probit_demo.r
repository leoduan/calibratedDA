
setwd("~/git/ImbalancedPG/test/")

x<- -5
exp(-x^2)

1/x
exp(-x^2/2)

exp(-x^2/2)/x = 1/n


logn<-  seq(2,10,by = 0.1)
n_series<- 10^logn
  
v<- sapply(n_series, function(n){var(qnorm(rbeta(10000,50,n-50)))})

plot(logn,v ,type="l")
lines(logn,1/n_series)

# lines(logn,1/logn*v[1]/2.5)

cv100<- log(n_series)



require("ggplot2")

df<- data.frame("log10 n"=rep(logn,3),"Distribution"= as.factor( rep(c("Marginal","Original Conditional r=1","Calibriated Conditional r=1000"),each=length(logn))),"Variance Ratio"=c(v/v,1/n_series/v, pmin(1000/n_series/v,1)))

levels(df$Distribution)

df$Distribution<- ordered(df$Distribution, levels = c("Marginal", "Calibriated Conditional r=1000","Original Conditional r=1"))

ggplot(data=df, aes(x=log10.n, y=Variance.Ratio,linetype=Distribution,length(n_series)))+ geom_line(size=c(.75))


pdf("./probit_rate.pdf",6,3)
ggplot(data=df, aes(x=log10.n, y=Variance.Ratio,linetype=Distribution,length(n_series)))+ geom_line(size=c(.75))
dev.off()




require(truncnorm)

theta<- -3

n<- 1E4
# r<- n/log(n)

r<- 1E10
prob<- pnorm(theta)
# y<- runif(n)<prob
# C<- rep(0,n)
y<- rep(0,n)
y[1]<-1



probit_r<- function(r){
  
  trace_theta<- numeric()
  
  trace_theta_mean<- numeric()
  trace_theta_lb<- numeric()
  trace_theta_ub<- numeric()
  for(i in 1:1000){
    
    lb<- rep(-Inf,n)
    ub<- rep(Inf,n)
    lb[y==1]<- (1-sqrt(r))*theta
    ub[y==0]<- (1-sqrt(r))*theta
    
    z<- rtruncnorm(n, lb,ub, theta, sd=sqrt(r))
    
    
    theta_lb<- max( (z/(1- sqrt(r) ))[y==1])
    theta_ub<-min( (z/(1- sqrt(r) ))[y==0])
    
    
    if(r>1)
      theta <- rtruncnorm(1, a=  theta_lb, b= theta_ub, mean = sum(z)/n, sd=sqrt(r/n))
    else
      theta<- rnorm(1,mean(z),sqrt(1/n))
    
    trace_theta<- c(trace_theta,theta)
    trace_theta_mean<- c(trace_theta_mean,mean(z))
    trace_theta_lb<-  c(trace_theta_lb, theta_lb)
    trace_theta_ub<-  c(trace_theta_ub, theta_ub)
    print(i)
  }
  
  trace_theta
}

acf(trace_theta,lag.max = 40)
ts.plot(trace_theta)

sd(trace_theta)
sd(qnorm(rbeta(1E5,1,n-1)))

ts.plot(trace_theta)


hist(trace_theta_lb, breaks=30)
hist(trace_theta_ub, add=T, breaks=30,col="green")
hist(trace_theta_mean,add=T,col="red",breaks = 30)

ts.plot(trace_theta_ub- trace_theta_lb)

fit1 <-probit_r(1)
fit2 <-probit_r(10)
fit3 <-probit_r(100)
fit4 <-probit_r(1000)
fit5 <-probit_r(10000)
# fit6 <-probit_r(100000)

acf_all <- apply(cbind(fit1,fit2,fit3,fit4,fit5),2,function(x){ acf(x, lag.max = 40, plot=F)$acf})

df<- data.frame("ACF"=c(acf_all),"r value"=rep(c("1 (Albert-Chib 1993)","10","100","1000","10000"),each=41),"Lag"=rep(c(0:40),5))

pdf("./probit_demo_acf.pdf",5,3)
ggplot(data=df, aes(x=Lag, y=ACF,linetype=r.value))+ geom_line(size=.75)
dev.off()


