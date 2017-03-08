data {
  int<lower=0> N; //
  int<lower=0> p; //
  int<lower=0> y[N]; // number of schools
  matrix[N,p] X;
}
parameters {
  vector[p] beta; 
  real<lower=0, upper=1> prob;
}
model {
    vector[N] mu;
    mu= exp(X * beta);
  
    for (i in 1:N) {
    if (y[i] == 0) {
      // not present
      // Bernoulli(0|θ) + Bernoulli(1|θ) * Poisson(0|λ)
      target += log(    (1 - prob) + prob * exp(-mu[i])   );
    } else {
      // present
      // Bernoulli(1|θ) * Poisson(y|λ)
      target +=  (bernoulli_lpmf(1 |  prob)
                    + poisson_lpmf(y[i] |  mu[i]));
    }
  }
	for(i in 1:p){
		beta[i] ~ normal(0, 1000);
	}
	prob ~ beta(1,1);
}

