require("scalableDA")

poissonCDA<- function(y,X, r_ini= 1,tune= 100,burnin=500, run=500 ,fixR = FALSE, randomEffect= FALSE,sigma2= 0.1, MH= T , C=1E4,c=1){
  
  n<- length(y)
  p<- ncol(X)
  
  #initialize beta
  beta<- rep(1,p)
  Xbeta<- X%*%beta
  
  # initialize random effect
  # sigma2<- 1
  eta<- rnorm(n,0,sd=sqrt(sigma2))
  eta[1]<- 0
  if(!randomEffect)
    eta<- rep(0,n)
  
  #initialize r and b
  
  r<- rep(r_ini,n)
  b<- rep(0,n)
  
  #individual log-likelihood
  
  #big constant used to approximate (1+exp(Xbeta)/C)^C approx exp(exp(Xbeta))
  logC<- log(C)
  
  loglik <-   - exp(Xbeta + eta)
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
    
    # #use the first half of burnin steps to adaptively tune r and b
    if(tuning){
      
      if(!fixR){
        if(  max_loglik< sum(loglik)){
          
          dprob<- exp(Xbeta+eta)
          curPGmean<- C/2/abs(Xbeta+eta + b - logC)*tanh(abs(Xbeta+eta + b - logC)/2)
          r0<-       dprob/  curPGmean
          }
        r<- r0*c
        r[r>1]<-1
        r[r*C<y]<- (y/C)[r*C<y]
        b = log(exp(exp(Xbeta + eta - logC -log(r)))-1) -Xbeta -eta + logC
        b[is.infinite(b)]=-log(r[is.infinite(b)])
        
        max_loglik = sum(loglik)
        
      }
      
      
    }
    
    
    # generate latent variable
    Z<- rpg(r*C, abs(Xbeta+eta + b - logC))
    # generate proposal:
    
    V<- solve(t(X) %*% ((Z)* X))
    m<- V%*%( t(X) %*% (y- (r*C)/2 -Z*(eta + b - logC)))
    
    cholV<- t(chol(V))
    
    new_beta <- cholV%*% rnorm(p) + m
    new_Xbeta <- X%*%new_beta
    
    new_loglik <- -exp(new_Xbeta + eta)
    
    ###
    q_loglik<-    -C*r *log(1+exp(Xbeta+b + eta -logC)) 
    new_q_loglik<-   -C*r *log(1+exp(new_Xbeta + b + eta -logC)) 
    
    # inidividual likelihood ratio
    alpha<- (new_loglik + q_loglik  - loglik - new_q_loglik  )
    
    alpha[is.na(alpha)]<- -10
    alpha[is.infinite(alpha)]<- 10
    
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
      print(tune_accept / tune_step)
      if(i>50)
        c = c* exp(  (0.6- (tune_accept / tune_step))/10)
      if(c<1)c=1
    }
    print(c)
    
    # 
    # #random effect
    # if(randomEffect){
    #   bigC<- 1E6
    #   # generate latent variable
    #   Z<- rpg(bigC+y, abs(Xbeta -log(bigC) + eta))
    #   # generate proposal:
    #   V<- 1/(Z + 1/sigma2)
    #   m<- V* (y- (y+bigC)/2 -Z*(-log(bigC)+ Xbeta))
    #   cholV<- sqrt(V)
    #   eta <- cholV * rnorm(n) + m
    #   eta[1]<- 0
    #   sigma2<- 1/rgamma(1, (n-1)/2 + 2, rate= sum(eta^2)/2+1)
    # }
    
    
    ##
    if(i> burnin){
      trace_proposal<- rbind(trace_proposal,t(new_beta))
      trace_beta<- rbind(trace_beta,t(beta))
      trace_sigma2<- c(trace_sigma2,sigma2)
      # trace_eta<- c(trace_eta,eta[1])
      
    }
    
    
    
    
    
    # print(exp(sum(alpha)))
    print(i)
  }
  return(list("beta"=trace_beta, "proposal"=trace_proposal , "r"=r, "b"=b, "accept_rate"= accept/run, "eta"=eta, "sigma2"=trace_sigma2))
  
  
  # return(list("beta"=trace_beta))
}