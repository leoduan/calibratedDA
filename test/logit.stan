data {
  int<lower=0> y; // number of schools 
  int<lower=0> n; // 
  real<lower=0> r; // 

}
parameters {
  real theta; 
}
transformed parameters {
  real<lower=0> prob;
  prob <- 1- 1.0/(1.0+ exp(theta)/r)^r;
}
model {
  theta ~ normal(0, 10000);
  // r~ normal(0, 1)T[0,];;
  y ~ binomial(n ,prob);
}