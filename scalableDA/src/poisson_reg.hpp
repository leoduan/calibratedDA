using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
SEXP poisson_reg(SEXP y, SEXP X, double r_ini = 1, int tune = 100,
                 int burnin = 500, int run = 500, bool fixR = false,
                 double C = 1E4, double c_ini = 1, bool MH = true) {
  C11RNG c11r;

  Rcpp::NumericVector yr(y);
  Rcpp::NumericMatrix Xr(X);
  const colvec y_vec(yr.begin(), yr.size(), false);
  const mat X_mat(Xr.begin(), Xr.nrow(), Xr.ncol(), false);

  int n = X_mat.n_rows;
  int p = X_mat.n_cols;

  vec beta = solve(X_mat.t() * X_mat, X_mat.t() * log(y_vec + 1));
  vec Xbeta = X_mat * beta;

  vec r = ones(n) * r_ini;
  vec b = zeros(n);

  double c = c_ini;

  double logC = log(C);

  vec loglik = -exp(Xbeta);
  double max_loglik = -INFINITY;

  mat trace_beta(run, p);
  mat trace_proposal(run, p);
  vec trace_accept_alpha(run);

  int accept = 0;
  int tune_accept = 0;
  vec r0 = ones(n);

  for (int i = 0; i < (burnin + run); ++i) {
    if (i < tune) {
      if (!fixR) {
        if (max_loglik < accu(loglik)) {
          vec dprob = exp(Xbeta);
          r0 = dprob / (tanh(abs(Xbeta + b - logC) / 2.0) / 2.0 /
                        abs(Xbeta + b - logC)) /
               C;
        }
        r = r0 * c;
        r(find(r > 1)).fill(1);
        r(find(r * C < y_vec)) = y_vec(find(r * C < y_vec)) / C;

        b = log(exp(exp(Xbeta - logC - log(r))) - 1.0) - Xbeta + logC;
        b(find_nonfinite(b)) = -log(r(find_nonfinite(b)));
        max_loglik = accu(loglik);
      }
    }

    vec alpha1 = r * C;
    vec alpha2 = abs(Xbeta + b - logC);
    vec Z = rpg(alpha1, alpha2);

    mat ZX_tilde = X_mat.each_col() % Z;
    mat invV = (X_mat.t() * ZX_tilde + eye(p, p) * 1E-5);
    mat cholInvV = chol(invV);
    mat V = inv(invV);

    vec k = y_vec - (r * C) / 2.0 - Z % (b - logC);

    vec m = V * (X_mat.t() * k);

    // vec new_beta = trans(chol(V)) * randn(p) + m;
    vec new_beta = inv(cholInvV) * randn(p) + m;
    vec new_Xbeta = X_mat * new_beta;

    // compute likelihood
    vec new_loglik = -exp(new_Xbeta);
    vec q_loglik = -C * r % log(1 + exp(Xbeta + b - logC));
    vec new_q_loglik = -C * r % log(1 + exp(new_Xbeta + b - logC));

    // metropolis-hastings

    vec alpha = new_loglik + q_loglik - loglik - new_q_loglik;
    alpha(find_nonfinite(alpha)).fill(10);

    if (log(c11r.runif()) < accu(alpha) || !MH) {
      beta = new_beta;
      Xbeta = new_Xbeta;
      loglik = new_loglik;
      if (i >= burnin) accept += 1;
      if (i < tune) tune_accept += 1;
    }

    if (i < tune) {
      cout << exp(accu(alpha)) << endl;
      cout << (double)tune_accept / (double)(i + 1) << endl;
      if (i > 50)
        c *= exp((0.6 - (double)tune_accept / (double)(i + 1)) / 10.0);
      if (c < 1) c = 1;
      if (!MH) c = c_ini;
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
