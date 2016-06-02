using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
SEXP logit_reg(SEXP y, SEXP X, SEXP b, SEXP B, int burnin = 500, int run = 500,
               double r0_ratio = 20, int mc_draws = 1E4, double r1 = 2.0,
               bool downsampling = true) {
  C11RNG c11r;

  Rcpp::NumericVector yr(y);
  Rcpp::NumericVector br(b);

  Rcpp::NumericMatrix Xr(X);
  Rcpp::NumericMatrix Br(B);

  colvec y_vec(yr.begin(), yr.size(), true);
  const colvec b_vec(br.begin(), br.size(), false);

  const mat X_mat(Xr.begin(), Xr.nrow(), Xr.ncol(), false);
  const mat B_mat(Br.begin(), Br.nrow(), Br.ncol(), false);

  int n = X_mat.n_rows;
  int p = X_mat.n_cols;
  vec beta = ones(p);

  mat B_inv = inv(B_mat);
  vec B_invb = solve(B_mat, b_vec);

  // set y=1
  uvec pos_group = find(y_vec == 1);
  int n_y1 = pos_group.n_elem;

  // set y=0
  uvec zero_group = find(y_vec == 0);

  // threshold: r_tilde
  double r_tilde = n_y1 / (double)y_vec.size() * r0_ratio;
  double log_r0 = log(r_tilde);

  if (r_tilde > 1) {
    stop("r_tilde>1");
  }

  // PG latent variable
  colvec w = ones(yr.size());

  // X corresponds to y=1
  // mat X1 = X_mat.rows(pos_group);
  // X corresponds to y=0
  // mat X0 = X_mat.rows(zero_group);

  mat trace_beta(run, p);
  vec filter_count(run);

  vec eta;

  for (int i = 0; i < (burnin + run); ++i) {
    // compute  Xbeta
    vec Xbeta = X_mat * beta;
    vec exp_Xbeta = exp(Xbeta);

    // keep those with exp(Xbeta) > r_tilde

    const double r_star = 1.1;//little wiggle room for the threshold
    uvec group_r_fixed_1 = find(exp_Xbeta > r_star * r_tilde);
    uvec pos_group_r1 =  find((exp_Xbeta <= r_star * r_tilde) && (y_vec == 1));
    uvec zero_group_r0 = find((exp_Xbeta <= r_star * r_tilde) && (y_vec == 0));

    // sample Monte Carlo weight, eta: keep zero_group_r0_subsample, throw away
    // other y=0.
    if (downsampling) {
      eta = c11r.draw_binomial(mc_draws, r_tilde, zero_group_r0.n_elem);
    } else {
      eta = ones(zero_group_r0.n_elem) * mc_draws * r_tilde; //eta/k set to equally r_tilde
    }

    uvec zero_group_r0_subsample = zero_group_r0(find(eta>0));
    vec eta_subsample = eta(find(eta > 0));

    // sample PG: w

    w.fill(0);
    vec alpha1 = zeros(n);
    vec alpha2 = zeros(n);

    // w for r = 1
    alpha1(group_r_fixed_1).fill(1);
    alpha2(group_r_fixed_1) = Xbeta(group_r_fixed_1);
    w(group_r_fixed_1) = rpg(alpha1(group_r_fixed_1), alpha2(group_r_fixed_1));

    // w for r = r1
    alpha1(pos_group_r1).fill(r1);
    alpha2(pos_group_r1) = Xbeta(pos_group_r1) - log(r1);

    w(pos_group_r1) = rpg(alpha1(pos_group_r1), alpha2(pos_group_r1));

    // w for r = r0 under downsampling

    alpha1(zero_group_r0_subsample) = eta_subsample / (double)mc_draws;
    alpha2(zero_group_r0_subsample) = Xbeta(zero_group_r0_subsample) - log_r0;

    w(zero_group_r0_subsample) =
        rpg(alpha1(zero_group_r0_subsample), alpha2(zero_group_r0_subsample));

    // sample beta

    uvec keep_idx = join_vert( join_vert(group_r_fixed_1, pos_group_r1), zero_group_r0_subsample);
    
    mat X123 = X_mat.rows(keep_idx);
    vec w123 = w(keep_idx);

    vec k1 = y_vec(group_r_fixed_1) - 0.5;
    vec k2 = y_vec(pos_group_r1) - r1 / 2.0 + w(pos_group_r1) * log(r1);
    vec k3 = - eta_subsample / (double)mc_draws / 2.0 + w(zero_group_r0_subsample) * log_r0;
    vec k123 = join_vert(join_vert(k1, k2), k3);

    mat wX_tilde = X123;
    wX_tilde.each_col() %= w123;

    mat V = inv(X123.t() * wX_tilde + B_inv);
    vec m = V * (X123.t() * k123 + B_invb);

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
