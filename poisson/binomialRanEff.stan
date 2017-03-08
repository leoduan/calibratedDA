data {
  int<lower=0> n;
  int<lower=0> y[n]; 
  int<lower=0> N[n]; 
}
parameters {
  real beta0;
  vector[n] beta;
  real<lower=0> sigma2;
}
model {
    vector[n] prob;
    real sigma;
    sigma=sqrt(sigma2);

    for(i in 1:n){
        prob[i] = 1/(1+exp(-beta[i]));
    }

	for(i in 1:n){
		beta[i] ~ normal(beta0, sigma);
	}
	
	y ~ binomial(N, prob);

    beta0 ~ normal(-12,7);
    sigma2 ~ uniform(0, 20);
}

