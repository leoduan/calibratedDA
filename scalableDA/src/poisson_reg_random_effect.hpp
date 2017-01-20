using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
SEXP poisson_reg_block_random(SEXP y, SEXP X, double r_ini = 1, int tune = 100,
                              int burnin = 500, int run = 500,
                              bool fixR = false, double C = 1E4,
                              double c_ini = 1, double c_ini2 = 1,
                              double theta_ts = 1, bool MH = true,
                              double nu_ini = 1, double sigma2 = 100,
                              bool adaptC = true, bool fixNu = false,
                              bool centeredRanEff = false) {
  C11RNG c11r;

  Rcpp::NumericVector yr(y);
  Rcpp::NumericMatrix Xr(X);
  const colvec y_vec(yr.begin(), yr.size(), false);
  const mat X_mat(Xr.begin(), Xr.nrow(), Xr.ncol(), false);

  int n = X_mat.n_rows;
  int p = X_mat.n_cols;

  vec Z = ones(n);

  vec L1;
  mat L2t, L3t;
  double nu = nu_ini;
  double eta0 = 0;

  vec beta, logtau;

  auto LtInvMultiply = [&](vec z) -> vec {
    vec z1 = z(span(0, n - 1));
    vec z2 = z(span(n, n + p - 1));
    vec b2 = solve(L3t, z2);
    vec b1 = 1.0 / L1 % z1 - 1.0 / L1 % (L2t * b2);
    return join_cols(b1, b2);
  };

  auto LInvMultiply = [&](vec z) -> vec {
    vec z1 = z(span(0, n - 1));
    vec z2 = z(span(n, n + p - 1));
    vec b1 = 1.0 / L1 % z1;
    vec b2 = solve(L3t.t(), -L2t.t() * b1 + z2);
    return join_cols(b1, b2);
  };

  // initialization
  {
    L1 = sqrt(Z + 1.0 / nu);
    L2t = X_mat.each_col() % (Z / L1);
    {
      vec temp = Z - Z % Z / (Z + 1.0 / nu);
      mat tempX = X_mat.each_col() % temp;
      mat L3L3 = trans(X_mat) * tempX + eye(p, p) / sigma2;
      L3t = chol(L3L3);
    }

    vec k = log(y_vec + 0.01);
    vec XtildeK = join_cols(k, X_mat.t() * k);

    L1 = sqrt(Z + 1.0 / nu);
    L2t = X_mat.each_col() % (Z / L1);
    {
      vec temp = Z - Z % Z / (Z + 1.0 / nu);
      mat tempX = X_mat.each_col() % temp;
      mat L3L3 = trans(X_mat) * tempX + eye(p, p) / sigma2;
      L3t = chol(L3L3);
    }
    vec m = LtInvMultiply(LInvMultiply(XtildeK));
    logtau = m(span(0, n - 1));
    beta = m(span(n, n + p - 1));
  }

  vec Xbeta = X_mat * beta;
  // vec logtau = randn(n) * sqrt(nu);

  vec r = ones(n) * r_ini;
  vec b = zeros(n);

  double c = c_ini;

  double logC = log(C);

  vec loglik = -exp(Xbeta + logtau);
  double max_loglik = -INFINITY;

  mat trace_beta(run, p);
  mat trace_tau(run, n);

  mat trace_proposal(run, p);
  vec trace_accept_alpha(run);
  vec trace_nu(run);
  vec trace_eta0(run);

  int accept = 0;
  int tune_accept = 0;
  vec r0 = ones(n);

  for (int i = 0; i < (burnin + run); ++i) {
    vec loglik = -exp(Xbeta + logtau);

    vec alpha1 = r * C;
    vec alpha2 = abs(Xbeta + logtau + b - logC);
    Z = rpg(alpha1, alpha2);

    // mat ZX_tilde = X_mat.each_col() % Z;
    // mat invV = (X_mat.t() * ZX_tilde + eye(p, p) / sigma2);
    // mat cholInvV = chol(invV);
    // mat V = inv(invV);

    vec k = y_vec - alpha1 / 2.0 + Z % (logC - b);
    vec k_ranEff = k;

    if (centeredRanEff) k_ranEff += eta0 / nu;

    vec XtildeK = join_cols(k_ranEff, X_mat.t() * k);

    L1 = sqrt(Z + 1.0 / nu);
    L2t = X_mat.each_col() % (Z / L1);
    {
      vec temp = Z - Z % Z / (Z + 1.0 / nu);
      mat tempX = X_mat.each_col() % temp;
      mat L3L3 = trans(X_mat) * tempX + eye(p, p) / sigma2;
      L3t = chol(L3L3);
    }
    vec m = LtInvMultiply(LInvMultiply(XtildeK));

    vec stdNormal = randn(n + p);

    vec new_theta = LtInvMultiply(stdNormal) + m;
    vec new_logtau = new_theta(span(0, n - 1));
    vec new_beta = new_theta(span(n, n + p - 1));

    vec new_Xbeta = X_mat * new_beta;

    // compute likelihood
    vec new_loglik = -exp(new_Xbeta + new_logtau);
    vec q_loglik = -C * r % trunc_log(1 + exp(Xbeta + logtau + b - logC));
    vec new_q_loglik =
        -C * r % trunc_log(1 + exp(new_Xbeta + new_logtau + b - logC));

    // metropolis-hastings

    vec alpha = (new_loglik + q_loglik - loglik - new_q_loglik);
    alpha(find_nonfinite(alpha)).fill(0);

    if (log(c11r.runif()) < accu(alpha) || !MH) {
      beta = new_beta;
      Xbeta = new_Xbeta;
      logtau = new_logtau;
      loglik = new_loglik;
      if (i >= tune) accept += 1;
      if (i < tune) tune_accept += 1;
    }

    if (!fixNu)
      nu = 1 /
           c11r.draw_gamma(n / 2.0 - 1,
                           accu((logtau - eta0) % (logtau - eta0)) / 2.0);
    if (centeredRanEff) eta0 = c11r.rnorm(mean(logtau), sqrt(nu / n));

    cout << exp(accu(alpha)) << endl;
    if (i < tune) {
      if (i % 30 == 0) {
        cout << "tuning accept rate " << (double)tune_accept / 30 << endl;
        tune_accept = 0;
      }

      // if (!MH) c = c_ini;
    } else {
      cout << (double)accept / (double)(i - tune + 1) << endl;
    }

    // cout << nu << endl;

    R_CheckUserInterrupt();
    cout << i << endl;

    if (adaptC) {
      if (i < burnin - 1)
        c = c_ini2 * c_ini;
      else
        c = c_ini;
    }

    if (i < tune) {
      if (!fixR) {
        vec dprob = exp(Xbeta + logtau);
        r0 = dprob / (C / 2.0 / abs(Xbeta + logtau - logC) %
                      tanh(abs(Xbeta + logtau - logC)));

        r = r0 * c;
        r(find(r > 1)).fill(1);
        uvec problem_set = find(r * C < y_vec);
        r(problem_set) = y_vec(problem_set) / C;

        // uvec big_theta = find(Xbeta + logtau > theta_ts);
        // r(big_theta) = r0(big_theta) * c_ini2;

        b = trunc_log(exp(exp(Xbeta + logtau - logC - trunc_log(r))) - 1.0) -
            ((Xbeta + logtau)) + logC;
        b(find_nonfinite(b)) = -trunc_log(r(find_nonfinite(b)));

        max_loglik = accu(loglik);
      }
    }

    if (i >= burnin) {
      trace_beta.row(i - burnin) = beta.t();
      trace_tau.row(i - burnin) = logtau.t();
      trace_nu(i - burnin) = nu;
      trace_eta0(i - burnin) = eta0;
      // trace_r.row(i - burnin) = r.t();
      // trace_w.row(i - burnin) = w.t();
    }
  }

  return Rcpp::List::create(
      Rcpp::Named("beta") = trace_beta, Rcpp::Named("tau") = trace_tau,
      Rcpp::Named("nu") = trace_nu, Rcpp::Named("r") = r,
      Rcpp::Named("eta0") = trace_eta0,
      Rcpp::Named("acceptance_rate") = (double)accept / run,
      Rcpp::Named("c") = c);
  // Rcpp::Named("r") = trace_r,
  // Rcpp::Named("w") = trace_w
}
