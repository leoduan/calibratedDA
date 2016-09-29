
pbeta(0.5,3,2)

r<-1.1
y<-2

hb<-function(p){
  gamma(r+2)/gamma(r-y+1)/gamma(y+1) *p^y*(1-p)^(r-y)
}

hb(0.3)
dbeta(0.3,  y+1, r-y+1)


approx_d<-function(psy,y,r){
  p<- exp(psy)/r / (1+exp(psy)/r)
  dbeta(p,  y+1, r-y+1)/pbeta(0.5,  y+1, r-y+1)/p/(1-p)
}

exact_d<-function(psy,y,r){
  dgamma(exp(psy),y+1,rate=1)/pgamma(r,y+1,rate = 1)/abs(psy)
}




psy0<- -3
y0<- 0
r<- 100

approx_d(psy0,y0,r)

exact_d(psy0,y0,r)
