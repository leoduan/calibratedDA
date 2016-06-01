using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
SEXP binomial(SEXP n, SEXP y, double b = 0, double B = 1000, int burnin = 500,
              int run = 500, double r_ratio = 20) {
  RNG r;
#ifdef USE_R
  GetRNGstate();
#endif

  Rcpp::NumericVector nr(n);
  Rcpp::NumericVector yr(y);

  colvec n_vec(nr.begin(), nr.size(), false);
  colvec y_vec(yr.begin(), yr.size(), false);

  // n_vec = n_vec/y_vec * 1;
  // y_vec.fill(1);
  // colvec power_reducing = y_vec / 1;
  // colvec B_vec = B * power_reducing;
  double B_vec = B;

  colvec r_vec = y_vec / n_vec * r_ratio;

  // PG latent variable
  colvec w = ones(nr.size());
  colvec theta = ones(nr.size());
  colvec var = ones(nr.size());
  colvec m = ones(nr.size());

  colvec log_r = log(r_vec);
  colvec rn = r_vec % n_vec;

  colvec k = y_vec - rn / 2.0;

  // prepare trace_matrix
  mat trace_theta(run, nr.size());
  mat trace_w(run, nr.size());

  for (int i = 0; i < (burnin + run); ++i) {
    var = 1.0 / (w + 1.0 / B_vec);
    m = var % (k + w % log_r + b / B_vec);
    r.norm(theta, 0.0, 1.0);
    // cout<<psi.t()<<endl;
    theta = theta % sqrt(var) + m;
    w = rpg(rn, theta - log_r);
    if (i >= burnin) {
      trace_theta.row(i - burnin) = theta.t();
      trace_w.row(i - burnin) = w.t();
    }
  }

#ifdef USE_R
  PutRNGstate();
#endif

  // Rcpp::NumericVector(out.begin(),out.end());

  return Rcpp::List::create(Rcpp::Named("theta") = trace_theta,
                            Rcpp::Named("w") = trace_w);
}
