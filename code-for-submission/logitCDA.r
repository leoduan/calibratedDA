require("BayesLogit")

logitCDA<- function(y,X,N=1,tune= 200, burnin=500, run=500,fixR = FALSE,r_ini= 1,MH= TRUE,c0= 1){
  
  n<- length(y)
  p<- ncol(X)
  
  #initialize beta
  #use MLE as the starting value
  model=glm(y ~ X-1,family=binomial(link=logit))
  beta=model$coefficients
  
  Xbeta<- X%*%beta
  prob<- plogis(Xbeta)
  
  #initialize r and b
  
  r<- rep(r_ini,n)
  b<- rep(0,n)
  
  # generate latent variable
  Z<- rpg(n,r*N, abs(Xbeta+b))
  # generate proposal:
  
  V<- solve(t(X) %*% ((Z)* X))
  m<- V%*%( t(X) %*% (y- N*r/2 -Z*b))
  cholV<- t(chol(V))
  beta <- cholV%*% rnorm(p) + m
  Xbeta <- X%*%beta
  new_Xbeta<- Xbeta
  #individual log-likelihood
  loglik <- - N*log(1+exp(Xbeta))

  #objects to store trace
  trace_beta<- numeric()
  trace_proposal<- numeric()

  accept<- 0
  tune_accept<- 0
  tune_step<- 0
  
  tuning<- TRUE
  
  for(i in 1: (burnin+run)){

    if(tuning){
      
      if(!fixR){
        
          dprob<- exp(Xbeta)/(1+exp(Xbeta))^2
          r0<-  dprob/(tanh(abs(Xbeta+b)/2)/2/abs(Xbeta+b)) 
          r<- r0/c0
          r[r>1]<-1
          r[r*N<y]<- (y/N)[r*N<y]
          
          b = log( (1+ exp(Xbeta))^{1/r}-1) -Xbeta
        }
    }
    
    # generate latent variable
    Z<- rpg(n,r*N, abs(Xbeta+b))

    # generate proposal:

    Vinv<- (t(X) %*% ((Z)* X))
    V<- solve(Vinv)
    m<- V%*%( t(X) %*% (y- N*r/2 -Z*b))
    
    cholV<- solve(chol(Vinv))
    
    new_beta <- cholV%*% rnorm(p) + m
    new_Xbeta <- X%*%new_beta
    new_loglik<-  - N*log(1+exp(new_Xbeta))  
    
    ###
    q_loglik<-    -N*r*log(1+exp(Xbeta+b)) 
    new_q_loglik<-  -N*r*log(1+exp(new_Xbeta+b)) 
    
    # inidividual likelihood ratio
    alpha<- new_loglik + q_loglik  - loglik - new_q_loglik  

    alpha[is.na(alpha)]<- -6
    alpha[is.infinite(alpha)]<- 6
    

    #metropolis-hastings
    if(log(runif(1))< sum(alpha)|| !MH ){  
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
      if(tune_step>tune) {
        if(tune_accept / tune_step<0.2)
          print("acceptance rate is too low, consider reducing c0")
        tuning = FALSE
      }
      # print(c("Acceptance Rate: ", tune_accept / tune_step))
    }
    
    if(i> burnin  & !tuning){
      trace_proposal<- rbind(trace_proposal,t(new_beta))
      trace_beta<- rbind(trace_beta,t(beta))
    }
    
    print(i)
  }
  return(list("beta"=trace_beta, "proposal"=trace_proposal ,
              "r"=r, "b"=b, "accept_rate"= accept/run))
  
  
}