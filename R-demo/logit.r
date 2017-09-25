require("scalableDA")

logitCDA<- function(y,X,N=1, r_ini= 1,tune =200,burnin=500, run=500 ,fixR = FALSE,c=1, MH= FALSE){
  
  n<- length(y)
  p<- ncol(X)
  
  #initialize beta
  beta<- rep(1,p)
  Xbeta<- X%*%beta
  prob<- plogis(Xbeta)
  
  #initialize r and b
  
  r<- rep(r_ini,n)
  b<- rep(0,n)
  
  # generate latent variable
  Z<- rpg(N, abs(Xbeta+b))
  # generate proposal:
  
  V<- solve(t(X) %*% ((Z)* X))
  m<- V%*%( t(X) %*% (y- N*r/2 -Z*b))
  cholV<- t(chol(V))
  beta <- cholV%*% rnorm(p) + m
  Xbeta <- X%*%beta
  new_Xbeta<- Xbeta
  #individual log-likelihood
  loglik <- - N*log(1+exp(Xbeta)) #/ max(N)
  max_loglik<- -Inf
  # loglik<- -Inf
  #objects to store trace
  trace_beta<- numeric()
  trace_proposal<- numeric()
  trace_importance<- numeric()
  
  
  accept<- 0
  tune_accept<- 0
  tune_step<- 0
  
  tuning<- TRUE
  
  for(i in 1: (burnin+run)){
    
    
    
    # #use the first half of burnin steps to adaptively tune r and b
    if(tuning){
      
      if(!fixR){
        
        if(  max_loglik< sum(loglik)){
          dprob<- exp(Xbeta)/(1+exp(Xbeta))^2
          r0<-       dprob/(tanh(abs(Xbeta+b)/2)/2/abs(Xbeta+b)) 
          r<- r0
          r[r>1]<-1
          # r[r*N<y]<- (y/N)[r*N<y]
          b = log( (1+ exp(Xbeta))^{1/r}-1) -Xbeta
          max_loglik = sum(loglik)
          # print(max_loglik)
        }
        if(MH)
        {
          r<- r0*c
          r[r>1]<-1
          r[r*N<y]<- (y/N)[r*N<y]
          
          b = log( (1+ exp(Xbeta))^{1/r}-1) -Xbeta
        }
      }
    }
    # generate latent variable
    Z<- rpg(r*N, abs(Xbeta+b))

    # generate proposal:

    Vinv<- (t(X) %*% ((Z)* X))
    V<- solve(Vinv)
    m<- V%*%( t(X) %*% (y- N*r/2 -Z*b))
    
    
    
    cholV<- solve(chol(Vinv))
    
    new_beta <- cholV%*% rnorm(p) + m
    new_Xbeta <- X%*%new_beta
    # new_prob <- plogis(new_Xbeta)
    new_loglik<-  - N*log(1+exp(new_Xbeta))   #/ max(N)
    
    ###
    q_loglik<-    -N*r*log(1+exp(Xbeta+b))  #/ max(N)
    new_q_loglik<-  -N*r*log(1+exp(new_Xbeta+b))  #/ max(N)
    
    # inidividual likelihood ratio
    alpha<- (new_loglik + q_loglik  - loglik - new_q_loglik  )
    importance<-  sum(new_loglik  - new_q_loglik)
    
    alpha[is.na(alpha)]<- -6
    alpha[is.infinite(alpha)]<- 6
    
    # alpha<- Inf
    # alpha<- exp(new_loglik +   - loglik   )
    # 
    # fixing NA and INF issue
    
    
    
    #metropolis-hastings
    if(log(runif(1))< sum(alpha)|| !MH ){  
      beta<- new_beta
      Xbeta<- new_Xbeta
      # prob<- new_prob
      loglik<- new_loglik
      if(i>= burnin)      accept<- accept +1
      if(tuning) {
        tune_accept<- tune_accept+1
      }
    }
    
    if(tuning){
      tune_step<- tune_step+1
      if( tune_step>=tune) tuning = FALSE
      # c = c * exp((0.4 - tune_accept / tune_step)/10)
      # if(c<1)c=1
      
      print(tune_accept / tune_step)
      
    }
    # 
    if(i> burnin){
      trace_proposal<- rbind(trace_proposal,t(new_beta))
      trace_beta<- rbind(trace_beta,t(beta))
      trace_importance<- c(trace_importance, importance)
      
    }
    
    # print(sum(alpha))
    print(i)
  }
  return(list("beta"=trace_beta, "proposal"=trace_proposal ,
              "r"=r, "b"=b, "accept_rate"= accept/run, "importance"=trace_importance,
              "c"=c))
  
  
}