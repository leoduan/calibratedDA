require("monomvn")
require("statmod")

N<- 2
p<- 100

X<- matrix(rnorm(N*p),N)
corX<- matrix(1/5,p,p)
diag(corX)<- 1
X<- t(t(chol(corX))%*% t(X))


beta<- rep(0,p)
beta[1:floor(p/5)]<- rt(floor(p/5),0,df = 2)
lambda<- 1000


y<- X%*%beta + rt(N,df = 4)
y_star<-   y -mean(y)
Xy_star<- t(X)%*%y_star


X2<- t(X)%*%X

sigma2<- 1
beta<- matrix(1,p)

trace_beta<- numeric()
trace_sigma2<- numeric()


for(i in 1:1000){
  tau_inv <- rinvgauss(p, lambda*sqrt(sigma2)/abs(beta)
                       , dispersion= 1/lambda^2)
  
  sigma2<- 1/rgamma(1, (N+p-1)/2,  rate=sum((y-X%*%beta)^2)/2 + sum(beta* tau_inv*beta)/2)
  
  A<-  X2+ diag(tau_inv)
  invA<- solve(A)
  m<- invA %*%Xy_star
  beta<-  t(chol(invA)) %*% rnorm(p)+m
  if(i>100){
    trace_beta<- rbind(trace_beta, t(beta))
    trace_sigma2<- c(trace_sigma2,sigma2)
  }
  
}

acf(trace_beta[,2])
acf(trace_sigma2)

# plot(colMeans(trace_beta))
# sigma2
