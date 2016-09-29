using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
SEXP poisson_reg(SEXP y, SEXP X, SEXP b, SEXP B, int burnin = 500,
                 int run = 500, double r0ini = 10, double c = 1,
                 int fixed_R = false) {
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
  vec beta = solve(X_mat.t() * X_mat, X_mat.t() * log(y_vec + 1));

  mat B_inv = inv(B_mat);
  vec B_invb = solve(B_mat, b_vec);

  // PG latent variable
  colvec w = ones(yr.size());

  // X corresponds to y=1
  // mat X1 = X_mat.rows(pos_group);
  // X corresponds to y=0
  // mat X0 = X_mat.rows(zero_group);

  mat trace_beta(run, p);
  mat trace_r(run, n);
  mat trace_w(run, n);

  vec r = ones(n);
  vec log_r = log(r);

  vec max_xbeta = -ones(n) * INFINITY;

  double r0 = r0ini;

  for (int i = 0; i < (burnin + run); ++i) {
    // compute  Xbeta
    vec Xbeta = X_mat * beta;

    // double r0 = c;
    // double r0 = sqrt(accu(expBeta)) / c;
    // if (r0 < 10) r0 = 10;

    if (fixed_R) {
      r.fill(1000);
      log_r = log(r);
      r(find(y_vec > 1000)) = y_vec(find(y_vec > 1000)) * 10;
    } else {
      if (i < burnin) {
        vec expBeta = exp(max(Xbeta, c * Xbeta));
        r = expBeta * r0;
      } else {
        max_xbeta = max(max_xbeta, max(Xbeta, Xbeta * c));
        r = exp(max_xbeta) * r0;
      }
      log_r = log(r);
    }

    // r(find(Xbeta > 20)) = pow(expBeta(find(Xbeta > 20)), 1.5) * r0;

    // r(find(r>1)) = exp( Xbeta(fid(r>1)) *2 ) *c;

    /*
    uvec setExpBetaGt1 = find_finite(Xbeta);
    // uvec setExpBetaGt1 = find(Xbeta>0);
    vec log_r2 = -log(
        exp(expBeta(setExpBetaGt1) / r(setExpBetaGt1) - Xbeta(setExpBetaGt1)) -
        1.0 / expBeta(setExpBetaGt1));
    // r(find(Xbeta > 0)) = exp(Xbeta(find(Xbeta > 0)) * 2) * r0;
        log_r(setExpBetaGt1) = log_r2;
    */

    vec alpha1 = ones(n) % r;
    vec alpha2 = abs(Xbeta - log_r);
    w = rpg(alpha1, alpha2);

    if (w.has_nan()) {
      throw std::runtime_error("w nan, consider reduce c & increase r0ini");
      cout << r0 << endl;
      // beta = ones(p);
      // r0 *= 1.1;
      // i = 0;
      continue;
    }

    mat wX_tilde = X_mat;
    wX_tilde.each_col() %= w;

    vec k = y_vec - r / 2.0 + w % log_r;

    mat V = inv(X_mat.t() * wX_tilde + B_inv);
    vec m = V * (X_mat.t() * k + B_invb);

    beta = trans(chol(V)) * randn(p) + m;

    if (beta.has_nan()) {
      throw std::runtime_error("beta nan, consider reduce c & increase r0ini");
      cout << r0 << endl;
      // beta = ones(p);
      // r0 *= 1.1;
      // i = 0;
      continue;
    }

    cout << i << endl;

    if (i >= burnin) {
      trace_beta.row(i - burnin) = beta.t();
      trace_r.row(i - burnin) = r.t();
      trace_w.row(i - burnin) = w.t();
    }
  }

  return Rcpp::List::create(Rcpp::Named("beta") = trace_beta,
                            Rcpp::Named("r") = trace_r,
                            Rcpp::Named("w") = trace_w);
}
