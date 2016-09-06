

rLogit<- function(y,m,xtheta,r){
  a<- xtheta-log(r)
  gamma(r*m+1)/gamma(y+1)/gamma(r*m-y+1) * exp(y*a)/(1+exp(a))^(r*m)
}


rLogit(c(0:3),3,-5,0.1)

