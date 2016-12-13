require("BayesLogit")

logitCDA<- function(y,X, r_ini= 1,burnin=500, run=500 ,fixR = FALSE){
  
  n<- length(y)
  p<- ncol(X)
 
  #initialize beta
  beta<- rep(1,p)
  Xbeta<- X%*%beta
  prob<- plogis(Xbeta)
  
  #initialize r and b
  
  r<- rep(r_ini,n)
  b<- rep(0,n)
  
  #individual log-likelihood
  loglik <- - log(1+exp(Xbeta)) 
  
  #objects to store trace
  trace_beta<- numeric()
  trace_proposal<- numeric()
  
  accept<- 0
  tune_accept<- 0
  tune_step<- 0
  
  tuning<- TRUE
  
  for(i in 1: (burnin+run)){
    
    # generate latent variable
    Z<- rpg(n, r, abs(Xbeta+b))
    # generate proposal:
    
    V<- solve(t(X) %*% ((Z)* X))
    m<- V%*%( t(X) %*% (y- r/2 -Z*b))
    cholV<- t(chol(V))
    
    new_beta <- cholV%*% rnorm(p) + m
    new_Xbeta <- X%*%new_beta
    # new_prob <- plogis(new_Xbeta)
    new_loglik<-  - log(1+exp(new_Xbeta)) 
    
    ###
    q_loglik<-    -r*log(1+exp(Xbeta+b)) 
    new_q_loglik<-  -r*log(1+exp(new_Xbeta+b)) 
    
    # inidividual likelihood ratio
    alpha<- exp(new_loglik + q_loglik  - loglik - new_q_loglik  )
    # alpha<- exp(new_loglik +   - loglik   )
    # 
    # # fixing NA and INF issue
    alpha[is.na(alpha)]<- 0.00001
    alpha[is.infinite(alpha)]<- 1
    # 
    # #use the first half of burnin steps to adaptively tune r and b
    if(tuning){
    #   
    #   # set bias correction to the value that L(y|xbeta) = L_{r,b}(y|xbeta)
    #   # b<- (sqrt(r)-1) * new_Xbeta
    #   
      if(!fixR){
        
        #
        dprob<- exp(Xbeta)/(1+exp(Xbeta))^2
        # r<-      ( dprob^2 /prob/(1 - prob)/(tanh(abs(Xbeta)/2)/2/(abs(Xbeta))) )
        r<-       dprob/(tanh(abs(Xbeta+b)/2)/2/abs(Xbeta+b)) 
        
        r[Xbeta> -2]<- 1
      
      }
      
      b = log( (1+ exp(Xbeta))^{1/r}-1) -Xbeta
      
    }
    
    #metropolis-hastings
    if(runif(1)< prod(alpha)){  
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
      if(tune_accept / tune_step > 0.3 & tune_step>=100) tuning = FALSE
      # if(tune_accept / tune_step > 0.3 ) tuning = FALSE
      print(tune_accept / tune_step)
    }
    # 
    # 
    if(i> burnin){
      trace_proposal<- rbind(trace_proposal,t(new_beta))
      trace_beta<- rbind(trace_beta,t(beta))
    }
    
    print(i)
  }
  return(list("beta"=trace_beta, "proposal"=trace_proposal , "r"=r, "b"=b, "accept_rate"= accept/run))
  
  
  # return(list("beta"=trace_beta))
}