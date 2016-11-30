require("truncnorm")

probitCDA<- function(y,X, r_ini=200,burnin=500, run=500 ,fixR = FALSE){
  
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
  Xbeta <- X%*%new_beta
  
  for(i in 1: (burnin+run)){
    
    # generate proposal:
    
    # 1. generate latent variable
    newZ<- rtruncnorm(n,a= lb,b= ub,mean = Xbeta + b, sd = sqrt(r))
    densityNewZ<- log(dtruncnorm(newZ, a=lb, b=ub, mean = Xbeta + b, sd= sqrt(r)))
    
    V<- solve(t(X) %*% ((1/r)* X))
    m<- V%*%( t(X) %*% ((1/r)*(newZ-b)))
    cholV<- t(chol(V))
    ind_normal<- rnorm(p)
    new_beta <- cholV%*% ind_normal + m
    new_Xbeta <- X%*%new_beta
    densityNewBeta<- sum(dnorm(ind_normal,log = T)) -0.5* log(det(V))
    
    # 2. compute reverse density
    
    m<- V%*%( t(X) %*% ((1/r)*(newZ-b)))
    ind_normal<- solve(cholV,(beta - m))
    densityOldBeta<- sum(dnorm(ind_normal,log = T)) -0.5* log(det(V))
    
    densityOldZ<- log(dtruncnorm(Z, a=lb, b=ub, mean = Xbeta + b, sd= sqrt(r)))
    

    new_prob <- pnorm(new_Xbeta)
    new_loglik<-  log(1-new_prob) * (y==0) + log(new_prob) * (y==1)
    
    
    q_prob <- pnorm((Xbeta+b)/sqrt(r))
    q_loglik<-  log(1-q_prob) * (y==0) + log(q_prob) * (y==1)
    
    new_q_prob <- pnorm((new_Xbeta+b)/sqrt(r))
    new_q_loglik<-  log(1-new_q_prob) * (y==0) + log(new_q_prob) * (y==1)
    
    
    # inidividual likelihood ratio
    alpha<- exp(new_loglik + q_loglik  - loglik - new_q_loglik  )
    # alpha<- exp(new_loglik  - loglik  )

    # if(i>burnin)
    # alpha<- exp(new_loglik + densityOldZ    - loglik - densityNewZ +  (densityOldBeta -  densityNewBeta )/n )
    
    
    # fixing NA and INF issue
    alpha[is.na(alpha)]<- 0.00001
    alpha[is.infinite(alpha)]<- 1
    
    #use the first half of burnin steps to adaptively tune r and b
    if(tuning){
      
      # set bias correction to the value that L(y|xbeta) = L_{r,b}(y|xbeta)
      # b<- (sqrt(r)-1) * new_Xbeta
      
      if(!fixR){
        # reduce r for those Xbeta in the region that don't have slow mixing problems (> -4) & having likelihood ratio<1.
        # r_adapt_set<-  ((alpha< 1) & (Xbeta> -4))
        r_adapt_set<- alpha<1
        r[r_adapt_set] <- r[r_adapt_set]* (alpha[r_adapt_set]^0.5)
        
        r_adapt_set<- (alpha>1)
        r[r_adapt_set] <- r[r_adapt_set]* (alpha[r_adapt_set]^0.5)

        # if any r falls under 1 (which mixes slower that the original Albert-Chib), put it back to 1
        r[r<1]<- 1
        r[r>5000]<- 5000
      }
    }
    b<- (sqrt(r)-1) * new_Xbeta
    
    
    #metropolis-hastings
    if(runif(1)< prod(alpha)){   # the transition kernel P(beta|new_beta) = P(new_beta|beta), so the MH ratio is just likelihood ratio
      beta<- new_beta
      Xbeta<- new_Xbeta
      prob<- new_prob
      loglik<- new_loglik
      Z<- newZ
      if(i>= burnin)      accept<- accept +1
      if(tuning) {
        tune_accept<- tune_accept+1
      }
    }
    
    if(tuning){
      tune_step<- tune_step+1
      if(tune_accept / tune_step > 0.3 & tune_step>100) tuning = FALSE
      
      print(tune_accept / tune_step)
      print(beta)
    }
    
    
    if(i>= burnin & !tuning){
      trace_proposal<- rbind(trace_proposal,t(new_beta))
      trace_beta<- rbind(trace_beta,t(beta))
    }
    
    print(i)
  }
  
  
  return(list("beta"=trace_beta, "proposal"=trace_proposal , "r"=r, "b"=b, "accept_rate"= accept/run))
}