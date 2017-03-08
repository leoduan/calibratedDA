data {
  int<lower=0> N; //
  int<lower=0> p; //
  int<lower=0> y[N]; 
  matrix[N,p] X;
}
parameters {
  vector[p] beta; 
}
transformed parameters {
  vector[N] mu;
  mu= exp(X * beta);
}
model {

	for (i in 1:N) {
		y[i] ~ poisson(mu[i]);
	}
	for(i in 1:p){
		beta[i] ~ normal(0, 1000);
	}
}

