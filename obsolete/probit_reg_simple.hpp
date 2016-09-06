using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
SEXP probit_reg_simple(SEXP y, SEXP X, SEXP b, SEXP B, int burnin = 500,
                       int run = 500, double r0 = 20, int mc_draws = 1E4,
                       double c = 1) {
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

  vec trace_r0(run);

  vec lb = -ones(n) * INFINITY;
  vec ub = ones(n) * INFINITY;
  uvec pos_set = find(y_vec == 1);
  uvec zero_set = find(y_vec == 0);
  vec r = ones(n) * r0;
  vec C = ones(n);
  vec Z = y_vec;

  r.fill(r0);

  lb(pos_set).fill(0);
  ub(zero_set).fill(0);

  for (int i = 0; i < (burnin + run); ++i) {
    // compute  Xbeta

    vec Xbeta = X_mat * beta;

    // r.fill(r0);
    // r(find(Xbeta > -2)).fill(1);

    // C = (1.0 - r) * c;

    Z = c11r.rtruncnorm(Xbeta, ones(n), lb, ub) % r;

    double a = (double)n / 2.0;
    vec diff = (Z - r % Xbeta);
    double b = dot(diff, diff) / 2.0;
    double r1 = sqrt(1.0 / c11r.draw_gamma(a, b));
    // if (r1 > r0 & c > 1) r0 = r1;
    // r.fill(r0);

    // vec a = ones(n) / 2 + 0.1;
    // vec b = diff % diff / 2.0 + 0.1;
    // r = sqrt(1.0 / c11r.draw_gamma(a, b));

    mat V = inv(X_mat.t() * X_mat + B_inv);

    vec m = V * (X_mat.t() * diagmat(1.0 / r) * Z + B_invb);

    mat cholV = trans(chol(V));

    beta = cholV * randn(p) + m;

    // Xbeta = X_mat * beta;
    // C = (1.0 - r) % Xbeta;
    R_CheckUserInterrupt();

    // uvec violate1 = find((Z <= C) && (y_vec == 1));
    // uvec violate0 = find((Z > C) && (y_vec == 0));

    // r(violate1).fill(1);  // = 1- Z(violate1)/Xbeta(violate1);
    // r(violate0).fill(1);  // = 1 - Z(violate0) / Xbeta(violate0);
    // r(find(r < 1)).fill(1);

    if (i >= burnin) {
      trace_beta.row(i - burnin) = beta.t();
      trace_r0(i - burnin) = r0;
    }
    cout << i << endl;
  }

  return Rcpp::List::create(Rcpp::Named("beta") = trace_beta,
                            Rcpp::Named("r") = r,
                            Rcpp::Named("trace_r0") = trace_r0);
}
