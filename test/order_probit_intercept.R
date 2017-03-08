#.rs.restartR()
require("ImbalancedPG")
setwd("~/git/ImbalancedPG/test/")
require(truncnorm)

n <- 1E4

X0 <- 1#rnorm(N, 1, 1)
# X1 <- rnorm(n, 1, 1)
# X <- cbind(X0, X1)
# beta <- c(-1, 1)
X<- matrix(rep(1,n))

beta<- 0.5
p<- length(beta)
# hist(X %*% beta)
Z <- rnorm(n,X %*% beta,1)

g <- c(0, 1, Inf) 
J= length(g)

y<- (Z<= g[1])*1 

for(j in 2:J){
  y<- y+  (Z<= g[j] & Z>g[j-1]) *j
}

table(y)


# beta<- c(1,1)


probit_r<- function(r,steps=1000){

  trace_beta<- numeric()
  trace_g<- numeric()
  trace_beta_ub<- numeric()
  trace_beta_lb<- numeric()
  
  for(i in 1:steps){
    
    lb<- rep(-Inf,n)
    ub<- rep(Inf,n)
    
    for(j in c(1:J)){
      if(j>1)
        lb[y==j]<- g[j-1] * sqrt(r) + (1-sqrt(r))*beta
      if(j< J)
        ub[y==j]<- g[j] * sqrt(r) + (1-sqrt(r))*beta
    }
    
    z<- rtruncnorm(n, lb,ub, X%*%beta, sd= sqrt(r))
    
    theta_lb<- -Inf
    theta_ub<- Inf
    for(j in c(1:(J-1))){
      if(r>1){
        theta_ub<- min(theta_ub, min( ((z - g[j]*sqrt(r))/(1- sqrt(r)))[y==j]))
        theta_lb<- max(theta_lb ,max( ((z - g[j]*sqrt(r))/(1- sqrt(r)))[y==j+1]))
        # theta_ub<- min(theta_ub, min((z/(1- sqrt(r) ))[y==j]))
        # theta_lb<- max(theta_lb, max(z/(1- sqrt(r) )[y==j+1]))
        # theta_ub<- Inf
        # theta_lb<- -Inf
      }
      if(r<1){
        theta_ub<- min(theta_ub, min( ((z - g[j]*sqrt(r))/(1- sqrt(r)))[y==j+1]))
        theta_lb<- max(theta_lb ,max( ((z - g[j]*sqrt(r))/(1- sqrt(r)))[y==j]))
      }
      if(r==1){
        theta_ub<- Inf
        theta_lb<- -Inf
      }
    }
    
    beta <- rtruncnorm(1, a=  theta_lb, b= theta_ub, mean = sum(z)/n, sd=sqrt(r/n))
    
    # if(i>(steps/2))
      {
      trace_beta<- rbind(trace_beta,t(beta))
      trace_g<- rbind(trace_g,t(g))
      trace_beta_ub<- rbind(trace_beta_ub, theta_ub)
      trace_beta_lb<- rbind(trace_beta_lb, theta_lb)
      
    }
    print(i)
  }
  
  list('beta'= trace_beta, 'g'=trace_g, "trace_beta_ub"=trace_beta_ub, "trace_beta_lb"=trace_beta_lb)
}

fit<- probit_r(r=1.01,100)
acf(fit$beta, lag.max = 40)
ts.plot(fit$beta)
ts.plot(cbind(fit$beta,fit$trace_beta_ub, fit$trace_beta_lb))


acf(fit$g[,2], lag.max = 40)
ts.plot(fit$g[,2])

ts.plot(fit$beta[,1])



s<- 1
b<- 0.64
trace_b<- numeric()
trace_b_mm<- numeric()
for(i in 1:1000){
  z<- rtruncnorm(200,0,2,0,s)
  s_z<- sort(z)
  b<- runif(1, s_z[100],s_z[101])
  trace_b<- c(trace_b,b)
  trace_b_mm<- rbind(trace_b_mm, c(s_z[100],s_z[101]))
}


ts.plot(trace_b)

trace_b2<- numeric()
trace_b2_mm<- numeric()

for(i in 1:1000){
  
  z1<- rtruncnorm(100,0,b,mean = 0, sd = s )
  z2<- rtruncnorm(100,b,2,mean = 0, sd = s )
  
  b<- runif(1, max(z1),min(z2))
  trace_b2<- c(trace_b2,b)
  trace_b2_mm<- rbind(trace_b2_mm, c(max(z1),min(z2)))
  
}

ts.plot(trace_b2, main="Albert-Chib")



mean(trace_b_mm[,2]-trace_b_mm[,1])
mean(trace_b2_mm[,2]-trace_b2_mm[,1])




png("thresholding.png",width = 8,height = 6,res = 300, units = "in")
par(mfrow=c(1,2))
ts.plot(trace_b2, main="Albert-Chib")
ts.plot(trace_b, main="Leo-James-David")
dev.off()

acf(trace_b2)
acf(trace_b)

ts.plot(trace_b2)



k<- 100
a1<- sapply(c(1:1000), function(x)max(rtruncnorm(k,0,0.64)))
a2<- sapply(c(1:1000), function(x)min(rtruncnorm(200-k,0.64,2)))
ts.plot(cbind(a1,a2))
ts.plot(runif(1000,a1,a2))
acf(runif(1000,a1,a2))
mean(a2-a1)

b1<- sapply(c(1:1000), function(x) (sort(rtruncnorm(200,0,2))[c(k,k+1)]))
ts.plot(runif(1000,min = b1[1,],max= b1[2,]))
acf(runif(1000,min = b1[1,],max= b1[2,]))


mean(b1[2,]-b1[1,])
ts.plot(t(b1))

mean(b1[1,])
mean(b1[2,])

