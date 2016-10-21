using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
SEXP probit_reg(SEXP y, SEXP X, SEXP b, SEXP B, int burnin = 500, int run = 500,
                double r0 = 20, int mc_draws = 1E4) {
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

  vec alpha = zeros(n);

  vec C = ones(n);
  vec Z = y_vec;

  lb(pos_set).fill(0);
  ub(zero_set).fill(0);

  r.fill(r0);
  // r(r < 1).fill(1);

  auto genBeta = [&]() {

    mat V = inv(X_mat.t() * diagmat(1.0 / r) * X_mat);

    vec m = V * (X_mat.t() * diagmat(1.0 / r) * (Z - alpha));

    mat cholV = trans(chol(V));

    return cholV * randn(p) + m;
  };

  // compute  Xbeta
  vec Xbeta = X_mat * beta;

  vec prob = pnorm_std(Xbeta);

  vec loglik = log(1 - prob) % (y_vec == 0) + log(prob) % (y_vec == 1);

  for (int i = 0; i < (burnin + run); ++i) {
    Z = c11r.rtruncnorm(Xbeta + alpha, sqrt(r), lb, ub);

    vec new_beta = genBeta();

    vec new_Xbeta = X_mat * new_beta;

    vec new_prob = pnorm_std(new_Xbeta);

    vec new_loglik =
        log(1 - new_prob) % (y_vec == 0) + log(new_prob) % (y_vec == 1);

    vec mh = exp(new_loglik - loglik);

    if (mh.has_nan() | mh.has_inf()) {
      uvec problem = find_nonfinite(mh);
      mh(problem).fill(1);
    }

    if (!std::isfinite(prod(mh))) {
      mh.fill(1);
    }

    // adaptively changes  r and b, in the first half of burnin

    if (i < burnin / 2) {
      uvec mh_less1 = find((mh < 1) && (Xbeta > -4));
      alpha = (sqrt(r) - 1) % new_Xbeta;
    }

    if (i < burnin / 2) {
      uvec mh_less1 = find((mh < 1) && (Xbeta > -4));
      r(mh_less1) %= pow(mh(mh_less1), 0.5);
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
                            Rcpp::Named("r") = r, Rcpp::Named("alpha") = alpha);
}
