N<-100000
p<-9

x<- matrix(rnorm(N*p),N,p)
beta<- rnorm(p,sd = 0.01)

s<-svd(x)

d<- s$d[s$d>1E-10]
u<- s$u[,s$d>1E-10]
v<- s$v[,s$d>1E-10]

xinv<- v %*% diag(1/d) %*% t(u)

hist(abs(xinv%*%x))

xbeta_max_norm <-max(abs(x%*%beta))
xinv_1_norm<- sum(abs(xinv))

beta_1_norm <- sum(abs(beta))

beta_1_norm
xinv_1_norm*xbeta_max_norm


beta_2_norm <- sum((beta)^2)
xbeta_2_norm <-sum((x%*%beta)^2)
xinv_2_norm<- sum((xinv)^2)

beta_2_norm
xbeta_2_norm*xinv_2_norm

X1<-apply(X, MARGIN = 2, function(x){
  (x-mean(x))/sd(x)
})
X1[,1]<-1

a<- svd(X1)
xinv<- a$v %*% diag(1/a$d) %*% t(a$u)
sum(abs(xinv))
max(abs(xinv))
sqrt(sum(abs(xinv)^2))
1/a$d
