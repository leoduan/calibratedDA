
sd(qnorm(rbeta(1000,1E5-1,1)))

sd(qnorm(rbeta(1000,1E5/2,1E5/2)))



sd((rbeta(1000,1E5-1,1)))

sd((rbeta(1000,1E5/2,1E5/2)))



logit(-5)

ilogit<- function(x){
  log(x/(1-x))
}


i<-0
N<- 1E3

var<- numeric(N)
for(i in (1:N)){
  var[i]<- sd(ilogit(rbeta(10000, i+1,N-i+1)))
}

plot(var)



var<- numeric(N)
for(i in (1:N)){
  var[i]<- sd(ilogit(rbeta(50000, i+1,N-i+1)))^2
}

plot(  log( -log(var)))



rtruncate_beta<- function(n,a,b, lb=0, ub =1){
  prob_lb <- pbeta(lb, a, b)
  prob_ub <- pbeta(ub, a, b)
  r<- runif(n)
  prob<- prob_lb + (prob_ub-prob_lb)*r
  qbeta(prob,a,b)
}


hist(rbeta(1E3,1,2))

sd( ilogit(rtruncate_beta(1000, 1,2, 0, 0.00001)))
sd( ilogit(rtruncate_beta(1000, 1,2, 0, 0.01)))

hist( ilogit(rtruncate_beta(1000, 2,1, 0, 1)))
hist( (rtruncate_beta(1000, 1,2, 0, 1)))


sd( ilogit(rtruncate_beta(1000, 1,2, 0.3, 0.4)))
sd( ilogit(rtruncate_beta(1000, 1,2, 0, 1)))



####

require("BayesLogit")

theta<- 10
r<- 0.1
m<- 1
mean(1/rpg(10000, m*r ,z = theta))
r<-1
mean(1/rpg(10000, m*r ,z = theta))

w<- rpg(10000, m*r ,z = theta -log(r))
mean(w)
hist(exp(- (theta -log(r))^2/2* w)/w)


var_pg<- function(b,c){
  b/4/c^3*(sinh(c)-c)/cosh(c/2)^2
}

var_pg(1,3)

var_pg(0.1,3-log(0.1))


plot(pgamma(seq(0,0.1,length.out = 1000), 0.8, 1))


(dgamma(1E-9,1.1,1)^100)/(1E-9*100)


r1<- 1

r2<- 1E-3
# hist(rpg(10000, r ,z = theta-log(r)))


mean(rpg(10000, r1 ,z = theta-log(r1)))/mean(rpg(10000, r2,z = theta-log(r2)))

theta<-3
r<-1
w<- rpg(1E4, r, theta-log(r))
sd(((1+w*log(r)))/w)



biasCorrection<- function(theta, r){
  log((1+exp(theta))^(1/r)-1)-theta
}


theta<- 10
log(0.8)

biasCorrection(theta,0.8)


pg_var<- function(b, c){
  c<-abs(c)
  b/4/c^3*(sinh(c)-c)/cosh(c/2)^2
}


require("BayesLogit")
sd(rpg(1000,5,3))^2

r<- 10

pg_var(5*r,3-log(r))
pg_var(5,3)

pg_var(1E-3,1)



r<- 1
sd(1/rpg(1E4,r, -10-log(r) ))
mean(1/rpg(1E4,r, -10-log(r) ))
hist(1/rpg(1E4,r, -10-log(r) ))


r<- 0.01
sd(1/rpg(1E4,r, -10-log(r) ))
mean(1/rpg(1E4,r, -10-log(r) ))


hist(1/rpg(1E4,r, 10-log(r) ))



lambda<- 10
r<- exp(-3*lambda)*3

-exp(lambda)

r*exp(3*lambda)*log(1/(1+exp(-lambda*2)/r))


approx<- function(t,r1){
  exp_t<- exp(t)
  r2<- 1/(exp(exp_t/r1-t)-exp(-t))
  r1* log(1+exp(t)/r2)
}

approx(3,2)
exp(3)


Xbeta <-  -10


r<- exp(Xbeta)*10
logr <- Xbeta - log( exp(log(1 +exp(Xbeta))/r)-1);

1/(1+exp(Xbeta -logr))^r - 1/(1+exp(Xbeta))
