require("BayesLogit")

logitCDA<- function(y,X, r_ini=200,burnin=500, run=500 ,fixR = FALSE){
  
  n<- length(y)
  p<- ncol(X)
 
  #initialize beta
  beta<- rep(1,p)
  Xbeta<- X%*%beta
  # prob<- pnorm(Xbeta)
  
  #initialize r and b
  
  r<- rep(r_ini,n)
  b<- rep(0,n)
  
  #individual log-likelihood
  # loglik <- log(1-prob) * (y==0) + log(prob) * (y==1)
  
  
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
    
    beta<- cholV%*% rnorm(p) + m
    Xbeta <- X%*%beta
    
    b<- -log(r)
    
    # new_beta <- cholV%*% rnorm(p) + m
    # new_Xbeta <- X%*%new_beta
    # new_prob <- pnorm(new_Xbeta)
    # new_loglik<-  log(1-new_prob) * (y==0) + log(new_prob) * (y==1)
    
    # inidividual likelihood ratio
    # alpha<- exp(new_loglik- loglik)
    # 
    # # fixing NA and INF issue
    # alpha[is.na(alpha)]<- 0.00001
    # alpha[is.infinite(alpha)]<- 1
    # 
    # #use the first half of burnin steps to adaptively tune r and b
    # if(tuning){
    #   
    #   # set bias correction to the value that L(y|xbeta) = L_{r,b}(y|xbeta)
    #   # b<- (sqrt(r)-1) * new_Xbeta
    #   
    #   if(!fixR){
    #     # reduce r for those Xbeta in the region that don't have slow mixing problems (> -4) & having likelihood ratio<1.
    #     # r_adapt_set<-  ((alpha< 1) & (Xbeta> -4))
    #     r_adapt_set<- alpha<1
    #     r[r_adapt_set] <- r[r_adapt_set]* (alpha[r_adapt_set]^0.5)
    #     
    #     r_adapt_set<- (alpha>1)
    #     r[r_adapt_set] <- r[r_adapt_set]* (alpha[r_adapt_set]^0.5)
    # 
    #     # if any r falls under 1 (which mixes slower that the original Albert-Chib), put it back to 1
    #     r[r<1]<- 1
    #     r[r>5000]<- 5000
    #   }
    # }
    # b<- (sqrt(r)-1) * new_Xbeta
    # 
    # 
    # #metropolis-hastings
    # if(runif(1)< prod(alpha)){   # the transition kernel P(beta|new_beta) = P(new_beta|beta), so the MH ratio is just likelihood ratio
    #   beta<- new_beta
    #   Xbeta<- new_Xbeta
    #   prob<- new_prob
    #   loglik<- new_loglik
    #   if(i>= burnin)      accept<- accept +1
    #   if(tuning) {
    #     tune_accept<- tune_accept+1
    #   }
    # }
    # 
    # if(tuning){
    #   tune_step<- tune_step+1
    #   if(tune_accept / tune_step > 0.3 & tune_step>100) tuning = FALSE
    #   
    #   print(tune_accept / tune_step)
    # }
    # 
    # 
    if(i>= burnin){
      # trace_proposal<- rbind(trace_proposal,t(new_beta))
      trace_beta<- rbind(trace_beta,t(beta))
    }
    
    print(i)
  }
  # return(list("beta"=trace_beta, "proposal"=trace_proposal , "r"=r, "b"=b, "accept_rate"= accept/run))
  
  
  return(list("beta"=trace_beta))
}