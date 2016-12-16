data {
  int<lower=0> N; //
  int<lower=0> p; //
  int<lower=0> y[N]; 
  matrix[N,p] X;
}
parameters {
  vector[p] beta1;
  vector[p] beta2;
}
transformed parameters {
  vector[N] mu1;
  vector[N] mu2;

  mu1= X * beta1;
  mu2= X * beta2;

}
model {
	
	for (j in 1:N) {
       if (y[j] == 0)
	      	 target +=( log_sum_exp( bernoulli_lpmf( 1|inv_logit(mu1[j] ) ),
						bernoulli_lpmf(0|inv_logit(mu1[j] ))
						+ poisson_log_lpmf( y[j]|mu2[j])));
       else
		target +=( bernoulli_lpmf(0|inv_logit(mu1[j] ))
						+ poisson_log_lpmf( y[j]|mu2[j])) ;
  }


	for(i in 1:p){
		beta1[i] ~ normal(0, 1000);
		beta2[i] ~ normal(0, 1000);
	}
}

