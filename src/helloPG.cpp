#include "RcppArmadillo.h"
#include <R_ext/Utils.h>
#include <iostream>
#include <exception>
#include "RNG.h"

#include "PolyaGamma.h"
#include "PolyaGammaAlt.h"
#include "PolyaGammaSP.h"
#include "PolyaGammaSmallB.h"
// Rcpp::depends(RcppArmadillo)

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
  for(int i=0; i<d; i++) {
    double b = shape(i);
    double z = scale(i);

    if(shape(i)>170)
    {
      double m = dv.pg_m1(b,z);
      double v = dv.pg_m2(b,z) - m*m;
      result[i] =  r.norm(m, sqrt(v));
    }
    else if (b > 13) {
      sp.draw(result[i], b, z, r);
    }
    else if (b==1 || b==2) {
      result[i] = dv.draw((int)b, z, r);
    }
    else if (b > 1) {
      result[i] = alt.draw(b, z, r);
    }
    else if (b > 0) {
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

  return Rcpp::NumericVector(out.begin(),out.end());
}


// [[Rcpp::export]]
SEXP multinomial(SEXP n, SEXP y, double b=0, double B=1000,
  int burnin = 500, int run = 500, double r_ratio = 20) {


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
  double B_vec= B;

  colvec r_vec = y_vec / n_vec * r_ratio;



  //PG latent variable
  colvec w = ones(nr.size());
  colvec psi = ones(nr.size());
  colvec theta = ones(nr.size());
  colvec var = ones(nr.size());
  colvec m = ones(nr.size());

  colvec log_r = log(r_vec);
  colvec rn = r_vec % n_vec;

  colvec k = y_vec - rn /2.0;

  //prepare trace_matrix
  mat trace_theta(run, nr.size());
  mat trace_w(run, nr.size());


  for (int i = 0; i < (burnin + run); ++i)
  {
    var= 1.0 / ( w +1.0 / B_vec);
    m = var % (k + b/ B_vec);
    r.norm(psi, 0.0, 1.0);
      // cout<<psi.t()<<endl;
    psi = psi  % sqrt(var) + m;
    w= rpg(rn, psi);
    theta = psi + log_r;
    if(i >= burnin){
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

