data {
  int<lower=0> N; //
  int<lower=0> p; //
  int<lower=0> y[N]; 
  matrix[N,p] X;
}
parameters {
  vector[p] beta;
  real<lower=0> sigma2;
  vector[N] eta;
}
transformed parameters {
  vector[N] mu;
  real sigma;
  mu= exp(X * beta + eta);
  sigma= sqrt(sigma2); 
}
model {

	for (i in 1:N) {
		y[i] ~ poisson(mu[i]);
	}
	for (i in 1:N){
	  eta[i] ~ normal(0, sigma);
	}
	for(i in 1:p){
		beta[i] ~ normal(0, 1000);
	}
	sigma2 ~ inv_gamma(2, 1);
}

