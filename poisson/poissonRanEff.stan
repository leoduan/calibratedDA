data {
  int<lower=0> N; //
  int<lower=0> p; //
  int<lower=0> y[N]; 
  matrix[N,p] X;
}
parameters {
  vector[p] beta;
  real<lower=0> sigma2;
}
transformed parameters {
  
  real sigma;
  sigma= sqrt(sigma2); 
}
model {
  vector[N] mu;
  vector[N] eta;

	for (i in 1:N){
	  eta[i] ~ normal(0, sigma);
	}
	for(i in 1:p){
		beta[i] ~ normal(0, 1000);
	}
	
	mu= X * beta + eta;
	
	sigma2 ~ inv_gamma(2, 1);
	
	for (i in 1:N) {
		y[i] ~ poisson_log(mu[i]);
	}
}

