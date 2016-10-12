using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
SEXP probit_reg_simple(SEXP y, SEXP X, SEXP b, SEXP B, int burnin = 500,
                       int run = 500, double r0 = 20, int mc_draws = 1E4) {
  C11RNG c11r;

  Rcpp::NumericVector yr(y);
  Rcpp::NumericVector br(b);

  Rcpp::NumericMatrix Xr(X);
  Rcpp::NumericMatrix Br(B);

  const colvec y_vec(yr.begin(), yr.size(), false);
  const colvec b_vec(br.begin(), br.size(), false);

  const mat X_mat(Xr.begin(), Xr.nrow(), Xr.ncol(), false);
  const mat B_mat(Br.begin(), Br.nrow(), Br.ncol(), false);

  int n = X_mat.n_rows;
  int p = X_mat.n_cols;
  vec beta = ones(p);

  mat B_inv = inv(B_mat);
  vec B_invb = solve(B_mat, b_vec);

  // PG latent variable
  colvec w = ones(yr.size());

  // X corresponds to y=1
  // mat X1 = X_mat.rows(pos_group);
  // X corresponds to y=0
  // mat X0 = X_mat.rows(zero_group);

  mat trace_beta(run, p);
  mat trace_beta_prop(run, p);

  vec lb = -ones(n) * INFINITY;
  vec ub = ones(n) * INFINITY;
  uvec pos_set = find(y_vec == 1);
  uvec zero_set = find(y_vec == 0);
  vec r = ones(n) * r0;
  vec C = ones(n);
  vec Z = y_vec;

  r.fill(r0);
  // r(r < 1).fill(1);

  auto genBeta = [&]() {

    mat V = inv(X_mat.t() * diagmat(1.0 / r / r) * X_mat + B_inv);

    vec m = V * (X_mat.t() * diagmat(1.0 / r / r) * Z + B_invb);

    mat cholV = trans(chol(V));

    return cholV * randn(p) + m;
  };

  // compute  Xbeta
  vec Xbeta = X_mat * beta;

  vec prob = pnorm_std(Xbeta);

  vec loglik = (1 - prob) % (y_vec == 0) + prob % (y_vec == 1);

  for (int i = 0; i < (burnin + run); ++i) {
    C = (1.0 - r) % Xbeta;
    // C = zeros(n);

    lb(pos_set) = C(pos_set);
    ub(zero_set) = C(zero_set);

    Z = c11r.rtruncnorm(Xbeta, r, lb, ub);

    // mat V = inv(X_mat.t() * diagmat(1.0 / r / r) * X_mat + B_inv);

    // vec m = V * (X_mat.t() * diagmat(1.0 / r / r) * Z + B_invb);

    // mat cholV = trans(chol(V));

    // beta = cholV * randn(p) + m;

    vec new_beta = genBeta();

    vec new_Xbeta = X_mat * new_beta;

    vec mu1 = (sqrt(2) / r) % ((0.5 - r) % Xbeta - 0.5 * new_Xbeta);
    vec pmu1 = pnorm_std(mu1);
    vec Qf = pmu1 % (y_vec == 0) + (1 - pmu1) % (y_vec == 1);

    vec mu2 = (sqrt(2) / r) % ((0.5 - r) % new_Xbeta - 0.5 * Xbeta);
    vec pmu2 = pnorm_std(mu2);
    vec Qb = pmu2 % (y_vec == 0) + (1 - pmu2) % (y_vec == 1);

    vec new_prob = pnorm_std(new_Xbeta);
    vec new_loglik = (1 - new_prob) % (y_vec == 0) + new_prob % (y_vec == 1);

    vec mh = new_loglik % Qb / loglik / Qf;

    if (mh.has_nan() | mh.has_inf()) {
      uvec problem = find_nonfinite(mh);
      // cout << trans(loglik(problem)) << endl
      // << trans(new_loglik(problem)) << endl
      // << trans(Qf(problem)) << endl
      // << trans(Qb(problem)) << endl;
      mh(problem) = new_loglik(problem) / loglik(problem);
    }

    if (i < burnin) {
      uvec mh_less1 = find(mh < 1);
      r(mh_less1) %= pow(mh(mh_less1), 1);
      // r %= exp(0.8 - mh);
      r(find(r < 1)).fill(1);
    }

    if (c11r.runif() < prod(mh) & mh.is_finite()) {
      beta = new_beta;
      Xbeta = new_Xbeta;
      prob = new_prob;
      loglik = new_loglik;
    }

    R_CheckUserInterrupt();

    if (i >= burnin) {
      trace_beta.row(i - burnin) = beta.t();
      trace_beta_prop.row(i - burnin) = new_beta.t();
    }
    cout << i << " " << prod(mh) << endl;
  }

  return Rcpp::List::create(Rcpp::Named("beta") = trace_beta,
                            Rcpp::Named("beta_prop") = trace_beta_prop,
                            Rcpp::Named("r") = r);
}
