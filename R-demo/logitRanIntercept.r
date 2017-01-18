require("scalableDA")

logitRanIntercept<- function(y,N=1, r_ini= 1,tune =200,burnin=500, run=500,fixR = FALSE,c=1, MH= FALSE, priorMean=0, priorVar=1E5){

  n<- length(y)
  p<- n

  #initialize beta
  beta<- rep(1,n)

  beta0<- priorMean
  sigma0<- priorVar

  #initialize r and b
  r<- rep(r_ini,n)
  b<- rep(0,n)

  # generate latent variable
  Z<- rpg(N, abs(beta+b))
  # generate proposal:

  V<- 1/(Z + 1/sigma0)
  m<- V* (y- N*r/2 -Z*b + beta0/sigma0)
  cholV<- sqrt(V)
  beta <- cholV  * rnorm(p) + m
  new_beta<- beta
  max_beta= rep(-Inf,n)

  #individual log-likelihood
  # loglik <- - N*log(1+exp(beta)) # Original
  loglik <- - N*log(1+exp(beta)) - (beta-beta0)^2/2/sigma0 # Modified
  max_loglik<- -Inf
  # loglik<- -Inf
  #objects to store trace
  trace_beta<- numeric()
  trace_proposal<- numeric()
  trace_importance<- numeric()
  trace_beta0<- numeric()
  trace_sigma0<- numeric()

  accept<- 0
  tune_accept<- 0
  tune_step<- 0

  tuning<- TRUE

  for(i in 1: (burnin+run)){

    # #use the first half of burnin steps to adaptively tune r and b
    if(tuning)
    {
      if(!fixR){

        if(  max_loglik< sum(loglik)){
          dprob<- exp(beta)/(1+exp(beta))^2
          r0<-       dprob/(tanh(abs(beta+b)/2)/2/abs(beta+b))
          # r0<- exp(beta)*c
          r<- r0 *c
          r[r>1]<-1
          r[r*N<y]<- (y/N)[r*N<y]
          # b= -log(r)
          b = log( (1+ exp(beta))^{1/r}-1) -beta
          b[is.infinite(b)]=-log(r[is.infinite(b)])
          max_loglik = sum(loglik)
          # print(max_loglik)
        }
      }
    }
    # generate latent variable
    Z<- rpg(r*N, abs(beta+b))

    # generate proposal:

    V<- 1/(Z + 1/sigma0)
    m<- V* (y- N*r/2 -Z*b + beta0/sigma0)
    cholV<- sqrt(V)
    new_beta <- cholV  * rnorm(p) + m
    # new_loglik<-  - N*log(1+exp(new_beta)) # Original
    new_loglik<-  - N*log(1+exp(new_beta)) - (new_beta-beta0)^2/2/sigma0 # Modified

    q_loglik<-    -N*r*log(1+exp(beta+b))  
    new_q_loglik<-  -N*r*log(1+exp(new_beta+b)) 

    # inidividual likelihood ratio
    alpha<- (new_loglik + q_loglik  - loglik - new_q_loglik  )
    importance<-  sum(new_loglik  - new_q_loglik)

    alpha[is.na(alpha)]<- -10
    alpha[is.infinite(alpha)]<- 10

    #metropolis-hastings
    # if(log(runif(1))< sum(alpha)|| !MH ){
    #   beta<- new_beta
    #   loglik<- new_loglik
    #   if(i>= burnin)      accept<- accept +1
    #   if(tuning) {
    #     tune_accept<- tune_accept+1
    #   }
    # }

    #point-wise metropolis-hastings
    if(MH){
      accept_set <-  log(runif(n))< alpha
    }else{
      accept_set<- rep(TRUE, n)
    }

    if( sum(accept_set)>0){
      beta[accept_set]<- new_beta[accept_set]
      loglik[accept_set]<- new_loglik[accept_set]
      if(i>= burnin)      accept<- accept +1
      if(tuning) {
        tune_accept<- tune_accept+1
      }
    }


    #update the hierarchical mean and variance
    v <- 1/(n/sigma0 + 1/ priorVar)
    m <- v* ( sum(beta)/sigma0 + priorMean/ priorVar)
    beta0 =  rnorm(1)* sqrt(v) + m

    a1 = n/2 + 2
    a2 = sum((beta-beta0)^2)/2+ 1
    sigma0 = 1/ rgamma(1, a1, rate=a2)

    if(tuning){
      tune_step<- tune_step+1
      if( tune_step>=tune) tuning = FALSE
      print(tune_accept / tune_step)

    }

    if(i> burnin){
      trace_proposal<- rbind(trace_proposal,t(new_beta))
      trace_beta<- rbind(trace_beta,t(beta))
      trace_importance<- c(trace_importance, importance)
      trace_beta0<- c(trace_beta0, beta0)
      trace_sigma0<- c(trace_sigma0, sigma0)

    }

    print(i)
  }
  return(list("beta"=trace_beta,
            "beta0"=trace_beta0,
            "sigma0"=trace_sigma0,
           "proposal"=trace_proposal ,
              "r"=r, "b"=b, "accept_rate"= accept/run, "importance"=trace_importance,
              "c"=c))

}
