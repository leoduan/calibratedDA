require("BayesLogit")

poissonCDA<- function(y,X, tune= 100,burnin=500, run=500 ,fixR = FALSE,r_ini= 1, MH= T , lambda=1E9,c0= 0.1, ini_beta= NULL){
  
  n<- length(y)
  p<- ncol(X)
  
  #initialize beta
  
  if(is.null(ini_beta)){
    beta<- rep(1,p)
  }else{
    beta = ini_beta
  }
  Xbeta<- X%*%beta
  
  
  #initialize r and b
  
  r<- rep(r_ini,n)
  b<- rep(0,n)
  
  #individual log-likelihood
  
  #big constant used to approximate (1+exp(Xbeta)/lambda)^lambda approx exp(exp(Xbeta))
  loglambda<- log(lambda)
  
  loglik <-   - exp(Xbeta )
  max_loglik<- -Inf
  
  #objects to store trace
  trace_beta<- numeric()
  trace_proposal<- numeric()
  trace_sigma2<- numeric()
  trace_eta<- numeric()
  
  accept<- 0
  tune_accept<- 0
  tune_step<- 0
  
  tuning<- TRUE
  
  for(i in 1: (burnin+run)){
    
    if(tuning){
    
      if(!fixR){

        dprob<- exp(Xbeta)
        curPGmean<- lambda/2/abs(Xbeta + b - loglambda)*tanh(abs(Xbeta + b - loglambda)/2)
        r0<-       dprob/  curPGmean
        
        r<- r0/c0
        # r[r>1]<-1
        r[r*lambda<y]<- (y/lambda)[r*lambda<y]
        b = log(exp(exp(Xbeta  - loglambda -log(r)))-1) -Xbeta  + loglambda
        b[is.infinite(b)]=-log(r[is.infinite(b)])
        
      }
    }
    
    
    # generate latent variable
    Z<- rpg(n, r*lambda, abs(Xbeta + b - loglambda))
    # generate proposal:
    
    V<- solve(t(X) %*% ((Z)* X))
    m<- V%*%( t(X) %*% (y- (r*lambda)/2 -Z*( b - loglambda)))
    
    cholV<- t(chol(V))
    
    new_beta <- cholV%*% rnorm(p) + m
    new_Xbeta <- X%*%new_beta
    
    new_loglik <- -exp(new_Xbeta )
    
    ###
    q_loglik<-    -lambda*r *log(1+exp(Xbeta+b  -loglambda)) 
    new_q_loglik<-   -lambda*r *log(1+exp(new_Xbeta + b  -loglambda)) 
    
    # inidividual likelihood ratio
    alpha<- (new_loglik + q_loglik  - loglik - new_q_loglik  )
    
    alpha[is.na(alpha)]<- -6
    alpha[is.infinite(alpha)]<- 6
    
    #metropolis-hastings
    if(log(runif(1))< sum(alpha) ||!MH){  
      beta<- new_beta
      Xbeta<- new_Xbeta
      loglik<- new_loglik
      if(i>= burnin)      accept<- accept +1
      if(tuning) {
        tune_accept<- tune_accept+1
      }
    }
    
    if(tuning){
      tune_step<- tune_step+1
      if(tune_step>= tune) tuning = FALSE
      print(c("Acceptance Rate: ", tune_accept / tune_step))
    }
    
    
    ##
    if(i> burnin){
      trace_proposal<- rbind(trace_proposal,t(new_beta))
      trace_beta<- rbind(trace_beta,t(beta))
    }
    
    print(i)
  }
  return(list("beta"=trace_beta, "proposal"=trace_proposal , "r"=r, "b"=b, "accept_rate"= accept/run))
  
}