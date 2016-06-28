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

  vec lb = -ones(n) * INFINITY;
  vec ub = ones(n) * INFINITY;
  uvec pos_set = find(y_vec == 1);
  uvec zero_set = find(y_vec == 0);
  vec r = ones(n) * r0;
  vec C = ones(n);
  
  for (int i = 0; i < (burnin + run); ++i) {
    // compute  Xbeta
    vec Xbeta = X_mat * beta;

    // vec tail_prob = 0.5 - erf(abs(Xbeta)/sqrt(2))/2;

    // r = 1.0 / tail_prob * r0;
    if (i < burnin) {
      r.fill(r0);
      r(find(abs(Xbeta) < 2)).fill(1);
      C = (1.0 - sqrt(r)) % Xbeta;
    }

    

    lb(pos_set) = C(pos_set);
    ub(zero_set) = C(zero_set);

    vec Z = c11r.rtruncnorm(Xbeta, sqrt(r), lb, ub);

    mat V = inv(X_mat.t() * diagmat(1.0 / r) * X_mat + B_inv);

    vec m = V * (X_mat.t() * diagmat(1.0 / r) * Z + B_invb);

    mat cholV = trans(chol(V));

    beta = cholV * randn(p) + m;

    if (i >= burnin) {
      trace_beta.row(i - burnin) = beta.t();
    }
    cout << i << endl;
  }

  return Rcpp::List::create(Rcpp::Named("beta") = trace_beta);
}
