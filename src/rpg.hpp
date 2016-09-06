using namespace Rcpp;
using namespace arma;

colvec rpg(colvec shape, colvec scale) {
  // C++-only interface to PolyaGamma class
  // draws random PG variates from arma::vectors of n's and psi's
  RNG r;
  PolyaGamma dv;
  PolyaGammaSmallB pg_smallb;
  PolyaGammaAlt alt;
  PolyaGammaSP sp;

#ifdef USE_R
  GetRNGstate();
#endif
  int d = shape.n_elem;
  colvec result(d);
  for (int i = 0; i < d; i++) {
    double b = shape(i);
    double z = scale(i);

    if (shape(i) > 170) {
      double m = dv.pg_m1(b, z);
      double v = dv.pg_m2(b, z) - m * m;
      result[i] = r.norm(m, sqrt(v));
    } else if (b > 13) {
      sp.draw(result[i], b, z, r);
    } else if (b == 1 || b == 2) {
      result[i] = dv.draw((int)b, z, r);
    } else if (b > 1) {
      result[i] = alt.draw(b, z, r);
    } else if (b > 0) {
      result[i] = pg_smallb.draw(shape(i), scale(i), r);
    }
  }
#ifdef USE_R
  PutRNGstate();
#endif
  return result;
}

// [[Rcpp::export]]
SEXP rpg(SEXP b, SEXP c) {
  Rcpp::NumericVector br(b);
  Rcpp::NumericVector cr(c);

  colvec pgscale(br.begin(), br.size(), false);
  colvec pgshape(cr.begin(), cr.size(), false);

  colvec out = rpg(pgscale, pgshape);

  return Rcpp::NumericVector(out.begin(), out.end());
}