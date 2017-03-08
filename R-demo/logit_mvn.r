require("BayesLogit")

logitMVN<- function(y,X, r_ini= 1,burnin=500, run=500 ,fixR = FALSE){
  
  n<- length(y)
  p<- ncol(X)
 
  #initialize beta
  beta<- rep(1,p)
  Xbeta<- X%*%beta
  prob<- plogis(Xbeta)
  
  #individual log-likelihood
  loglik <- -log(1+exp(Xbeta)) + y*exp(Xbeta)
  
  #objects to store trace
  trace_beta<- numeric()
  trace_proposal<- numeric()
  
  accept<- 0
  tune_accept<- 0
  tune_step<- 0
  r<-r_ini
  b<-0
  tuning<- TRUE
  
  FI<- 1
  
  for(i in 1: (burnin+run)){
    
    # generate proposal:
    
    V<-  solve(t(X) %*% (c(FI)* X))
    m<- beta
    cholV<- t(chol(V))
    new_beta <- r * cholV%*% rnorm(p) + m
    new_Xbeta <- X%*%new_beta
    # new_prob <- plogis(new_Xbeta)
    new_loglik<-  - log(1+exp(new_Xbeta)) + y*exp(new_Xbeta)
    
    
    
    # inidividual likelihood ratio
    alpha<- exp(new_loglik   - loglik   )
    # # alpha<- exp(new_loglik +   - loglik   )
    # # 
    # # # fixing NA and INF issue
    # alpha[is.na(alpha)]<- 0.00001
    # alpha[is.infinite(alpha)]<- 1
   
    if(runif(1)< prod(alpha)){  
      beta<- new_beta
      Xbeta<- new_Xbeta
      # prob<- new_prob
      loglik<- new_loglik
      if(i>= burnin)      accept<- accept +1
      # if(tuning) {
      #   tune_accept<- tune_accept+1
      # }
    }
    
    
    if(i<1000){
      FI<-  exp(Xbeta)/(1+exp(Xbeta))^2
      FI<- FI
      FI[FI<1E-1]<-1E-1
      # FI<- FI/alpha
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