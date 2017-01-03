data {
  int<lower=0> N; //
  int<lower=0> p; //
  int<lower=0> y[N]; 
  matrix[N,p] X;
}
parameters {
  vector[p] beta;
  real<lower=0> nu;
  vector<lower=0>[N] eta;
  // real eta0;

}
model {
  vector[N] mu;

	for (i in 1:N){
	  eta[i] ~ normal(0, sqrt(nu));
	}

	for(i in 1:p){
		beta[i] ~ normal(0, 100);
	}
	
	mu= X * beta + eta;
	
	for (i in 1:N) {
		y[i] ~ poisson_log(mu[i]);
	}

  nu ~ uniform(0, 1000);
  // eta0 ~ normal(0,100);
  
}

