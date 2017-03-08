using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
SEXP binomial_reg(SEXP y_r, SEXP n_r, int burnin = 500, int run = 500,
                  double r0 = 20, int mc_draws = 1E4) {
  Rcpp::NumericVector yr(y_r);
  Rcpp::NumericVector nr(n_r);

  const colvec y(yr.begin(), yr.size(), false);
  const colvec n(nr.begin(), nr.size(), false);

  vec theta = zeros(y.n_elem);
  vec w = ones(y.n_elem);

  mat trace_theta(run, y.n_elem);

  for (int i = 0; i < (burnin + run); ++i) {
    vec r = exp(theta) * r0;
    r(find(r > 1)).fill(1);

    vec alpha1 = r % n;
    vec alpha2 = theta - log(r);
    w = rpg(alpha1, alpha2);

    vec k = y - alpha1 / 2.0 + w % log(r);

    vec V = 1.0 / w;
    vec m = V % k;

    R_CheckUserInterrupt();

    theta = randn(y.n_elem) % sqrt(V) + m;

    if (i >= burnin) {
      trace_theta.row(i - burnin) = theta.t();
    }
  }

  return Rcpp::List::create(Rcpp::Named("theta") = trace_theta);
}
