require("truncnorm")

probitCDA<- function(y,X, burnin=500, run=500,fixR = FALSE,r_ini=1, MH= TRUE){
  
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
  trace_proposal<- numeric()

  accept<- 0
  tune_accept<- 0
  tune_step<- 0
  
  tuning<- TRUE
  
  Z<- rtruncnorm(n,a= lb,b= ub,mean = Xbeta + b, sd = sqrt(r))
  V<- solve(t(X) %*% ((1/r)* X))
  m<- V%*%( t(X) %*% ((1/r)*(Z-b)))
  cholV<- t(chol(V))
  ind_normal<- rnorm(p)
  beta <- cholV%*% ind_normal + m
  Xbeta <- X%*%beta
  
  for(i in 1: (burnin+run)){
    
    # generate proposal:
    
    Z<- rtruncnorm(n,a= lb,b= ub,mean = Xbeta + b, sd = sqrt(r))

    V<- solve(t(X) %*% ((1/r)* X))
    m<- V%*%( t(X) %*% ((1/r)*(Z-b)))
    cholV<- t(chol(V))
    ind_normal<- rnorm(p)
    new_beta <- cholV%*% ind_normal + m
    new_Xbeta <- X%*%new_beta

    new_prob <- pnorm(new_Xbeta)
    new_loglik<-  log(1-new_prob) * (y==0) + log(new_prob) * (y==1)
    
    q_prob <- pnorm((Xbeta+b)/sqrt(r))
    q_loglik<-  log(1-q_prob) * (y==0) + log(q_prob) * (y==1)
    
    new_q_prob <- pnorm((new_Xbeta+b)/sqrt(r))
    new_q_loglik<-  log(1-new_q_prob) * (y==0) + log(new_q_prob) * (y==1)
    
    
    # inidividual likelihood ratio
    alpha<- exp(new_loglik + q_loglik  - loglik - new_q_loglik  )

    # fixing NA and INF issue
    alpha[is.na(alpha)]<- 1E-8
    alpha[is.infinite(alpha)]<- 1
    
    #use the first half of burnin steps to adaptively tune r and b
    if(tuning){
      
      # set bias correction to the value that L(y;xbeta) = L_{r,b}(y;xbeta)
      # b<- (sqrt(r)-1) * new_Xbeta
      
      if(!fixR){
        dprob<- dnorm(Xbeta)
        r<-    c(prob*(1 - prob)/dprob^2)
        r[r<1]<- 1
      }
      b<- (sqrt(r)-1) * new_Xbeta
    }
    
    
    #metropolis-hastings
    if(runif(1)< prod(alpha)|| !MH){
      beta<- new_beta
      Xbeta<- new_Xbeta
      prob<- new_prob
      loglik<- new_loglik
      if(i>= burnin)      accept<- accept +1
      if(tuning) {
        tune_accept<- tune_accept+1
      }
    }
    
    if(tuning){
      tune_step<- tune_step+1
      if(tune_accept / tune_step > 0.3 & tune_step>200) tuning = FALSE
    }
    
    
    if(i> burnin & !tuning){
      trace_proposal<- rbind(trace_proposal,t(new_beta))
      trace_beta<- rbind(trace_beta,t(beta))
    }
    
    print(i)
  }
  
  
  return(list("beta"=trace_beta, "proposal"=trace_proposal , "r"=r, "b"=b, "accept_rate"= accept/run))
}