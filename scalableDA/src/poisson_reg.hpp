using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
SEXP poisson_reg(SEXP y, SEXP X, int tune = 100, int burnin = 500,
                 int run = 500, int fixR = false, double bigR = 20) {
  C11RNG c11r;

  Rcpp::NumericVector yr(y);
  Rcpp::NumericMatrix Xr(X);
  const colvec y_vec(yr.begin(), yr.size(), false);
  const mat X_mat(Xr.begin(), Xr.nrow(), Xr.ncol(), false);

  int n = X_mat.n_rows;
  int p = X_mat.n_cols;

  // vec beta = ones(p);
  vec beta = solve(X_mat.t() * X_mat, X_mat.t() * log(y_vec + 1));
  vec Xbeta = X_mat * beta;

  vec r = ones(n);
  vec b = zeros(n);

  vec loglik = -exp(Xbeta);

  mat trace_beta(run, p);

  int accept = 0;
  int tune_accept = 0;

  for (int i = 0; i < (burnin + run); ++i) {
    if (i < tune) {
      vec dprob = exp(Xbeta);
      r = dprob / (tanh(abs(Xbeta) / 2.0) / 2.0 / abs(Xbeta));
      uvec bigXbeta = find(Xbeta > -2);
      r(bigXbeta) = exp(Xbeta(bigXbeta)) * bigR;
      b = log(r) + log(exp(exp(Xbeta - log(y_vec + r))) - 1) - Xbeta;
    }

    // compute  Xbeta

    vec alpha1 = r + y_vec;
    vec alpha2 = abs(Xbeta + b - log(r));
    vec Z = rpg(alpha1, alpha2);

    mat ZX_tilde = X_mat.each_col() % Z;
    mat invV = (X_mat.t() * ZX_tilde + eye(p, p) * 1E-5);
    mat cholInvV = chol(invV);
    mat V = inv(invV);

    vec k = y_vec - (y_vec + r) / 2.0 - Z % (b - log(r));

    vec m = V * (X_mat.t() * k);

    // vec new_beta = trans(chol(V)) * randn(p) + m;
    vec new_beta = inv(cholInvV) * randn(p) + m;
    vec new_Xbeta = X_mat * new_beta;

    // compute likelihood
    vec new_loglik = -exp(new_Xbeta);
    vec q_loglik = -(y_vec + r) % log(1 + exp(Xbeta + b) / r);
    vec new_q_loglik = -(y_vec + r) % log(1 + exp(new_Xbeta + b) / r);

    // metropolis-hastings

    vec alpha = new_loglik + q_loglik - loglik - new_q_loglik;
    alpha(find_nonfinite(alpha)).ones();

    if (log(c11r.runif()) < accu(alpha)) {
      beta = new_beta;
      Xbeta = new_Xbeta;
      loglik = new_loglik;
      if (i >= burnin) accept += 1;
      if (i < tune) tune_accept += 1;
    }

    if (i < tune) {
      cout << exp(accu(alpha)) << endl;
      cout << (double)tune_accept / (double)(i + 1) << endl;
    }

    cout << i << endl;

    if (i >= burnin) {
      trace_beta.row(i - burnin) = beta.t();
      // trace_r.row(i - burnin) = r.t();
      // trace_w.row(i - burnin) = w.t();
    }
  }

  return Rcpp::List::create(
      Rcpp::Named("beta") = trace_beta,
      Rcpp::Named("acceptance_rate") = (double)accept / run);
  // Rcpp::Named("r") = trace_r,
  // Rcpp::Named("w") = trace_w
}
