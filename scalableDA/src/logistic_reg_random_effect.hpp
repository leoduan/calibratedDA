using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
SEXP logistic_reg_random_effect(SEXP y, SEXP X, int burnin = 500, int run = 500,
                                double tau = 10, double c = 1,
                                int mc_draws = 1E4, int da_ver = 1,
                                double sigma2_ini = 0.01, bool track_r = false,
                                bool update_sigma2 = false) {
  C11RNG c11r;

  Rcpp::NumericVector yr(y);

  Rcpp::NumericMatrix Xr(X);

  const colvec y_vec(yr.begin(), yr.size(), false);

  const mat X_mat(Xr.begin(), Xr.nrow(), Xr.ncol(), false);

  int n = X_mat.n_rows;
  int p = X_mat.n_cols;
  vec beta = solve(X_mat.t() * X_mat, X_mat.t() * log(y_vec + 1));
  double sigma2 = sigma2_ini;

  // PG latent variable
  colvec w = ones(n);

  mat trace_beta(run, p);
  mat trace_theta(run, n);
  mat trace_r;
  if (track_r) {
    trace_r = zeros<mat>(run, n);
  }
  mat trace_w(run, n);
  vec trace_sigma2(run);

  vec r, log_r;

  vec max_theta = -ones(n) * INFINITY;

  vec Xbeta = X_mat * beta;
  vec theta =
      Xbeta + randn(n) * sqrt(sigma2);  // ones with random effects xbeta+eta
  vec eta = theta - Xbeta;              // random effects eta

  r = ones(n) * tau;
  log_r = log(r);

  mat X_mat2 = X_mat.t() * X_mat;
  mat X_mat2inv = inv(X_mat2 + eye<mat>(p, p) * 1E-5);
  mat cholX_mat2inv = trans(chol(X_mat2inv));

  for (int i = 0; i < (burnin + run); ++i) {
    vec Xbeta = X_mat * beta;

    if (i > burnin / 10) {
      max_theta = max(theta, max_theta);
      r = max(exp(max_theta) * tau, exp(2 * max_theta) / 1000000);
    } else {
      r = max(exp(theta) * tau, exp(2 * theta) / 1000000);
    }

    r(find(r > 1)).fill(1);
    log_r = log(r);

    {
      vec alpha1 = r;
      vec alpha2 = abs(theta - log_r);
      w = rpg(alpha1, alpha2);

      uvec fail = find_nonfinite(w);
      w(fail) = alpha1(fail) / 2. / alpha2(fail) % tanh(alpha2(fail) / 2.0);
    }

    {
      vec var_theta = 1.0 / (w + 1.0 / sigma2);

      vec m_theta = var_theta % (w % log_r + Xbeta / sigma2 + y_vec - r / 2.0);

      theta = randn(n) % sqrt(var_theta) + m_theta;

      eta = theta - Xbeta;

      // vec var_eta = 1.0 / (w + 1.0 / sigma2);

      // vec m_eta = var_eta % (y_vec - r / 2.0 + w % (log_r - Xbeta));

      // eta = randn(n) % sqrt(var_eta) + m_eta;

      // theta = Xbeta + eta;
    }

    if (update_sigma2) {
      double a = (double)n / 2.0;
      vec diff = eta;
      double b = dot(diff, diff) / 2.0;
      sigma2 = 1.0 / c11r.draw_gamma(a, b);
    }

    if (da_ver == 1) {
      mat wX_tilde = X_mat;
      vec v = w - w % w / (w + 1 / sigma2);
      wX_tilde.each_col() %= v;
      mat V = (X_mat.t() * wX_tilde + eye<mat>(p, p) * 1E-8);

      vec k = y_vec - r / 2.0 + w % log_r;

      vec m = solve(V, (X_mat.t() * (k % (1 - w / (w + 1 / sigma2)))));

      mat cholV = chol(V);
      beta = solve(cholV, randn(p)) + m;
      Xbeta = X_mat * beta;
    }

    if (da_ver == 2) {
      mat wX_tilde = X_mat;
      wX_tilde.each_col() %= w;

      vec k = y_vec - r / 2.0 + w % (log_r - eta);

      mat V = (X_mat.t() * wX_tilde + eye<mat>(p, p) * 1E-8);
      vec m = solve(V, (X_mat.t() * k));

      mat cholV = chol(V);
      beta = solve(cholV, randn(p)) + m;
      Xbeta = X_mat * beta;
    }

    if (da_ver == 3) {
      beta = cholX_mat2inv * randn(p) + X_mat2inv * (X_mat.t() * theta);

      Xbeta = X_mat * beta;
    }

    // sigma2 = b / a;

    if (w.has_nan()) {
      throw std::runtime_error("w nan, consider reduce c & increase tau");
      cout << tau << endl;
      // beta = ones(p);
      // tau *= 1.1;
      // i = 0;
      continue;
    }

    if (beta.has_nan()) {
      throw std::runtime_error("beta nan, consider reduce c & increase tau");
      cout << tau << endl;
      // beta = ones(p);
      // tau *= 1.1;
      // i = 0;
      continue;
    }

    cout << i << endl;

    if (i >= burnin) {
      trace_beta.row(i - burnin) = beta.t();
      // trace_theta.row(i - burnin) = theta.t();
      if (track_r) {
        trace_r.row(i - burnin) = r.t();
      }
      // trace_w.row(i - burnin) = w.t();
      trace_sigma2(i - burnin) = sigma2;
    }
  }

  return Rcpp::List::create(
      Rcpp::Named("beta") = trace_beta, Rcpp::Named("sigma2") = trace_sigma2,
      Rcpp::Named("r") = r, Rcpp::Named("w") = w, Rcpp::Named("theta") = theta,
      Rcpp::Named("trace_r") = trace_r);
}
