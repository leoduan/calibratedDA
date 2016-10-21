require("truncnorm")

probitCDA<- function(y,X, r_ini=200,burnin=500, run=500){
  
  n<- length(y)
  p<- ncol(X)
  lb<- rep(-Inf,n)
  ub<- rep(Inf,n)
  lb[y==1]<-0
  ub[y==0]<-0
  
  #initialize beta
  beta<- rep(1,p)
  Xbeta<- X%*%beta
  prob<- pnorm(Xbeta)
  
  #initialize r and b
  
  r<- rep(r_ini,n)
  b<- rep(0,n)
  
  #individual log-likelihood
  loglik <- log(1-prob) * (y==0) + log(prob) * (y==1)
  
  
  #objects to store trace
  trace_beta<- numeric()
  
  
  for(i in 1: (burnin+run)){
    
    # generate latent variable
    Z<- rtruncnorm(n,a= lb,b= ub,mean = Xbeta + b, sd = sqrt(r))
    # generate proposal:
    
    V<- solve(t(X) %*% ((1/r)* X))
    m<- V%*%( t(X) %*% ((1/r)*(Z-b)))
    cholV<- t(chol(V))
    new_beta <- cholV%*% rnorm(p) + m
    new_Xbeta <- X%*%new_beta
    new_prob <- pnorm(new_Xbeta)
    new_loglik<-  log(1-new_prob) * (y==0) + log(new_prob) * (y==1)
    
    # inidividual likelihood ratio
    alpha<- exp(new_loglik- loglik)
    
    # fixing NA and INF issue
    alpha[is.na(alpha)]<- 0.001
    alpha[is.infinite(alpha)]<- 1
    
    #use the first half of burnin steps to adaptively tune r and b
    if(i < burnin /2){
      
      # set bias correction to the value that L(y|xbeta) = L_{r,b}(y|xbeta)
      b<- (sqrt(r)-1) * new_Xbeta
      
      # reduce r for those Xbeta in the region that don't have slow mixing problems (> -4) & having likelihood ratio<1.
      r_adapt_set<-  ((alpha<1) & (Xbeta> -4))
      r[r_adapt_set] <- r[r_adapt_set]* alpha[r_adapt_set]
      
      # if any r falls under 1 (which mixes slower that the original Albert-Chib), put it back to 1
      r[r<1]<- 1
    }
    
    #metropolis-hastings
    if(runif(1)< prod(alpha)){   # the transition kernel P(beta|new_beta) = P(new_beta|beta), so the MH ratio is just likelihood ratio
      beta<- new_beta
      Xbeta<- new_Xbeta
      prob<- new_prob
      loglik<- new_loglik
    }
    
    
    if(i>= burnin){
      trace_beta<- rbind(trace_beta,t(beta))
    }
    
    print(i)
  }
  
  
  return(list("beta"=trace_beta , "r"=r, "b"=b))
}