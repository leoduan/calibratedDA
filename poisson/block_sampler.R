# X<- matrix(rnorm(9),3)
# 
# X<- t(X)%*%X



n<-100
p<-3



nu<- 1
sigma2<- 1000
# w<- abs(rnorm(n))
w<- rep(10,n)

X<- matrix(rnorm(n*p),n)

L1 = sqrt(w+ 1/nu)
# L1inv= 1/L1
L2t = w/L1*X 
L3L3= t(X)%*% ((w-w*w/(w+ 1/nu))*X) + diag(1/sigma2,p)
L3t= chol(L3L3)

# L3t= chol(t(X)%*% (w*X) + diag(1/sigma2,p) - t(L2t)%*%L2t)

Lt= rbind(cbind(diag(L1),L2t), cbind(matrix(0,p,n),L3t))

V= t(Lt)%*%Lt
max(abs(diag(V[1:n,1:n])- (w+1/nu)))
max(abs((V[(n+1):(n+p),1:n]) - t(w*X)))
V[(n+1):(n+p),(n+1):(n+p)] - (t(X)%*%(w*X)+diag(1/sigma2,p))

z<- rnorm(n+p)

LtInvMultiply <-function(z){
  z1<- z[1:n]
  z2<- z[(n+1):(n+p)]
  b2<- solve(L3t,z2)
  b1<- 1/L1 * z1 - 1/L1*  L2t%*% b2
  b<- c(b1,b2)
  b
}

LInvMultiply <-function(z){
  z1<- z[1:n]
  z2<- z[(n+1):(n+p)]
  b1<- 1/L1 * z1
  b2<-solve(t(L3t), - t(L2t)%*% b1 + z2)
  b<- c(b1,b2)
  b
}

max(solve(Lt,z)-    LtInvMultiply(z))
max(solve(t(Lt),z)-    LInvMultiply(z))
max(abs(LtInvMultiply(LInvMultiply(z))- solve(V,z)))

#################

beta0= rnorm(p,0)
eta0= rnorm(n,0,sqrt(nu))
y<- c(X%*%beta0 +eta0 +rnorm(n,0,1/sqrt(w)))


m<- LtInvMultiply(LInvMultiply( c( w*y,t(X)%*% (w*y))))
betaEta<- LtInvMultiply(rnorm(n+p))+m

# betaEta
# beta0

plot(m, c(eta0,beta0))
