using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
SEXP logit_reg(SEXP y, SEXP X, SEXP b, SEXP B, int burnin = 500, int run = 500, double r0_ratio = 20, int mc_draws = 1E4, double r1 = 2.0, bool downsampling= true) {
  C11RNG c11r;

  Rcpp::NumericVector yr(y);
  Rcpp::NumericVector br(b);

  Rcpp::NumericMatrix Xr(X);
  Rcpp::NumericMatrix Br(B);

  colvec y_vec(yr.begin(), yr.size(), false);
  colvec b_vec(br.begin(), br.size(), false);

  mat X_mat(Xr.begin(), Xr.nrow(), Xr.ncol(), false);
  mat B_mat(Br.begin(), Br.nrow(), Br.ncol(), false);

  int n = X_mat.n_rows;
  int p = X_mat.n_cols;
  vec beta = ones(p);

  mat B_inv = inv(B_mat);
  vec B_invb = solve(B_mat, b_vec);

  // colvec w2 = ones(nr.size());

  // colvec theta = ones(nr.size());
  // colvec var = ones(nr.size());
  // colvec m = ones(nr.size());

  // set y=1
  uvec pos_group = find(y_vec == 1);

  // set y=0
  uvec zero_group = find(y_vec == 0);
  int n_y0 = zero_group.n_elem;

  // thresholding r_tilde
  double r_tilde = accu(y_vec) / (double)y_vec.size() * r0_ratio;
  double log_r = log(r_tilde);

  if (r_tilde > 1) {
    stop("r_tilde>1");
  }

  // PG latent variable
  colvec w = ones(yr.size());

  mat X_tilde1 = X_mat.rows(pos_group);

  mat trace_beta(run, p);
  vec filter_count(run);

  vec eta;
  for (int i = 0; i < (burnin + run); ++i) {
    // compute  Xbeta
    vec Xbeta = X_mat * beta;

    // sample Monte Carlo weight, eta: keep zero_group_subsample, throw away
    // other y=0.
    if(downsampling){
    eta = c11r.draw_binomial(mc_draws, r_tilde, n_y0);
    }else{
      eta = ones(n_y0) * mc_draws * r_tilde;
    }
    uvec zero_group_subsample = zero_group(find(eta > 0));

    // sample PG: w

    w.fill(0);

    vec alpha1 = zeros(n);
    vec alpha2 = zeros(n);

    // w for y=1
    alpha1(pos_group).fill(r1);
    alpha2(pos_group) = Xbeta(pos_group) - log(r1);

    w(pos_group) = rpg(alpha1(pos_group), alpha2(pos_group));

    // w for y=0 and eta>0

    alpha1(zero_group_subsample) = eta(find(eta > 0)) / (double)mc_draws;
    alpha2(zero_group_subsample) = Xbeta(zero_group_subsample) - log_r;

    w(zero_group_subsample) =
        rpg(alpha1(zero_group_subsample), alpha2(zero_group_subsample));

    // sample beta

    mat X_tilde2 = X_mat.rows(zero_group_subsample);
    mat X_tilde = join_vert(X_tilde1, X_tilde2);

    vec w1 = w(pos_group);
    vec w2 = w(zero_group_subsample);
    vec w12 = join_vert(w1, w2);

    vec k1 = y_vec(pos_group) - r1 / 2 + w1 * log(r1);
    vec k2 = -eta(find(eta > 0)) / (double)mc_draws / 2.0 + w2 * log_r;
    vec k12 = join_vert(k1, k2);

    mat wX_tilde = X_tilde;
    wX_tilde.each_col() %= w12;

    mat V = inv(X_tilde.t() * wX_tilde + B_inv);
    vec m = V * (X_tilde.t() * k12 + B_invb);

    beta = trans(chol(V)) * randn(p) + m;

    if (i >= burnin) {
      trace_beta.row(i - burnin) = beta.t();
      uvec filtered = find(eta == 0);
      filter_count(i - burnin) = (double)filtered.n_elem / (double)n;
    }
  }

  return Rcpp::List::create(Rcpp::Named("beta") = trace_beta,
                            Rcpp::Named("eta") = eta,
                            Rcpp::Named("filter_count") = filter_count);
}
