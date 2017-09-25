#.rs.restartR()
require("ImbalancedPG")
setwd("~/Documents/Projects/calibratedDA/test/")
require(truncnorm)

n <- 1E3

X0 <- 1#rnorm(N, 1, 1)
X1 <- rnorm(n, 1, 1)
X <- cbind(X0, X1)
beta <- c(1, 1)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


# X<- matrix(rep(1,n))
# beta<- 1

p<- length(beta)
# hist(X %*% beta)

Z <- rnorm(n,X %*% beta,1)

g <- c(0, 1,3, Inf) 
J= length(g)

y<- (Z<= g[1])*1 

for(j in 2:J){
  y<- y+  (Z<= g[j] & Z>g[j-1]) *j
}

table(y)

# beta<- c(1,1)
g<-  c(0, 0.1, 0.5, Inf) 

probit_r<- function(r,steps=1000, ver =1){
  
  trace_beta<- numeric()
  trace_g<- numeric()
  AUGY <- matrix(0,n,steps)
  
  for(i in 1:steps){
    
    lb<- rep(-Inf,n)
    ub<- rep(Inf,n)
    
    Xbeta <- X%*%beta
    
    
    count_y1<- sum(y==1)
    if(count_y1>0)
      Z[y==1] <- rtruncnorm(count_y1, mean= Xbeta[ y==1], b=g[1])
    
    
    if(ver==1){
      
      for(j in c(2:(J-1))){
        Z[y == J]  <- rtruncnorm(sum(y == J ),mean= Xbeta[y == J],a= max(Z[y== J-1]), b =Inf)
      }
      
      # s_z<- sort(Z[y>1],decreasing = F)
      
      
      # for(j in c(2:(J-1))){
      #   g[j]<- runif(1,  s_z[sum(y<=j)-count_y1],s_z[sum(y<=j)-count_y1+1])
      # }
      
      for(j in c(2:(J-1))){
        
        min_g<- max(Z[y==j])
        max_g<- min(Z[y==j+1])
        g[j]<- runif(1, min_g,max_g)
      }
      
      
    } 
    
    
    # if(ver==2){
    #   for(j in c(2:(J-1))){
    #     
    #     count_y1<- sum(y==j)
    #     count_y2<- sum(y==j+1)
    #     
    #     z<- rtruncnorm(count_y1+count_y2 ,mean= Xbeta[ (y==j) | (y== j+1)],a= g[j-1], b = g[j+1])
    #     
    #     s_z<- sort(z,decreasing = F)
    #     g[j]<- runif(1,s_z[count_y1],s_z[count_y1+1])
    #     
    #     
    #     Z[y==j]<- z[z < g[j]]
    #     Z[y==j+1]<- z [z > g[j]]
    #     
    #   }}
    
    if(ver==3){
      for(j in c(1:J)){
        if(j>1)
          lb[y==j]<- g[j-1] 
        if(j< J)
          ub[y==j]<- g[j] 
      }
      
      Z<- rtruncnorm(n, lb,ub, X%*%beta, sd= 1)
      
      for(j in c(2:(J-1))){
        if(j>1)
          min_g<- max( max(Z[y==j]), g[j-1])
        else
          min_g<- max(Z[y==1])
        if(j<J)
          max_g<- min(min(Z[y==j+1]), g[j+1])
        g[j]<- runif(1, min_g,max_g)
      }
      
    }
    
    
    XRX <- t(X)%*% X
    XRZ <- t(X) %*% Z
    m<- solve(XRX, XRZ)
    cholvar<- t(chol(solve(XRX)))
    beta<- cholvar%*% rnorm(p)+m
    

    # obtain the value of y from the value of z
    # this should always be right
    repg <- t(matrix(rep(g,n),length(g),n))
    repz <- matrix(rep(Z,length(g)),n,length(g))
    augy <- apply(repz>repg,1,sum)+1
    AUGY[,i] <- augy

    if(i>(steps/2)){
      trace_beta<- rbind(trace_beta,t(beta))
      trace_g<- rbind(trace_g,t(g))
    }
    
    #print(i)
  }
  
  list('beta'= trace_beta, 'g'=trace_g, 'AUGY'=AUGY)
}
# fit<- probit_r(1,10)

fit<- probit_r(1,1000, ver=1)
# fit2<- probit_r(1,1000, ver=2)
fit3<- probit_r(1,1000, ver=3)

nerr <- apply(fit$AUGY,2,function(x){sum(x!=y)})
nerr2 <- apply(fit2$AUGY,2,function(x){sum(x!=y)})
nerr3 <- apply(fit3$AUGY,2,function(x){sum(x!=y)}) 

# these should all be zero, avg # of times per iteration that the thresholded value of z \ne y
sum(nerr/1000)
sum(nerr2/1000)
sum(nerr3/1000)


acf(fit$beta, lag.max = 40)

acf(fit$g[,2], lag.max = 40)
acf(fit3$g[,2], lag.max = 40)

fit$g


da<- c(acf(fit3$g[, 2], lag.max = 40, main="DA",plot = F)$acf)
cda<- c(acf(fit$g[, 2], lag.max = 40, main="CDA",plot = F)$acf)

df<- data.frame("ACF"=c(da,cda),"Method"=rep(c("DA (Albert-Chib)","CDA"),each=41),"Lag"=rep(c(0:40),2))
df$Method<- ordered(df$Method, levels = c("DA (Albert-Chib)","CDA"))


require("ggplot2")
pdf("./ordered_probit_acf.pdf",5,3)
ggplot(data=df, aes(x=Lag, y=ACF,linetype=Method))+ geom_line(size=.75)
dev.off()


df<- data.frame("Value"=c(fit3$g[, 2],fit$g[, 2]),"Method"=rep(c("DA (Albert-Chib)","CDA"),each=500),"Iteration"=rep(c(1:500),2))
df$Method<- ordered(df$Method, levels = c("DA (Albert-Chib)","CDA"))

pdf("./ordered_probit_trace_plot.pdf",5,3)
ggplot(data=df, aes(x=Iteration, y=Value))+ geom_line(size=.75) + facet_grid(Method~.)
dev.off()

# par(mfrow=c(2,1))
# ts.plot(fit3$g[,2],ylab="Value",xlab="Iteration")
# ts.plot(fit$g[,2],ylab="Value",xlab="Iteration")




table(y)

g2<- numeric()
for(i in 1:1000){
  z<- rtruncnorm(1+8372 ,mean=beta,a=0)
  s_z<- sort(z)
  g2<- c(g2,runif(1,s_z[1],s_z[2]))
}

ts.plot(g2)
acf(g2)

sd(g2)
sd(fit$g[5000:10000,2])
