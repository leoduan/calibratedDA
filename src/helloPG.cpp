#include "RcppArmadillo.h"
#include <R_ext/Utils.h>
#include <iostream>
#include <exception>
#include "RNG.h"
#include "C11RNG.hpp"
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
SEXP binomial(SEXP n, SEXP y, double b=0, double B=1000,
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
    m = var % (k + w % log_r + b/ B_vec ) ;
    r.norm(theta, 0.0, 1.0);
      // cout<<psi.t()<<endl;
    theta = theta  % sqrt(var) + m;
    w= rpg(rn, theta - log_r);
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

/*
// [[Rcpp::export]]
SEXP binomial2(SEXP n, SEXP y, double b=0, double B=1000,
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


  //PG latent variable
  colvec w1 = ones(nr.size());
  colvec w2 = ones(nr.size());

  colvec theta = ones(nr.size());
  colvec var = ones(nr.size());
  colvec m = ones(nr.size());


  //prepare trace_matrix
  mat trace_theta(run, nr.size());
  mat trace_w1(run, nr.size());
  mat trace_w2(run, nr.size());


  for (int i = 0; i < (burnin + run); ++i)
  { 

    colvec r_vec = max( exp(theta)*1.2, y_vec / (n_vec - y_vec) * r_ratio);
    
    r_vec(find(r_vec>1)).fill(1);
    
    colvec log_r = log(r_vec);

    colvec rny = r_vec % (n_vec - y_vec);
    colvec y_minus_rny = y_vec - rny;

  // colvec k = y_vec - rny /2.0;

    var= 1.0 / ( w1 + w2 +1.0 / B_vec);
    m = var % ( y_minus_rny/2.0 + w2% log_r + b/ B_vec);
    r.norm(theta, 0.0, 1.0);
    theta = theta  % sqrt(var) + m;
    w1 = rpg(y_vec, theta);
    w2 = rpg(rny, theta - log_r);

    if(i >= burnin){
      trace_theta.row(i - burnin) = theta.t();
      trace_w1.row(i - burnin) = w1.t();
      trace_w2.row(i - burnin) = w2.t();
    }
  }
  
  #ifdef USE_R
  PutRNGstate();
  #endif

   // Rcpp::NumericVector(out.begin(),out.end());

  return Rcpp::List::create(
    Rcpp::Named("theta") = trace_theta,
    Rcpp::Named("w1") = trace_w1,
    Rcpp::Named("w2") = trace_w2
    );

}*/



// [[Rcpp::export]]
    SEXP logit_reg(SEXP y, SEXP X, SEXP b, SEXP B,
      int burnin = 500, int run = 500, double r_ratio = 20,
      int mc_draws= 1E4) {


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
      vec B_invb= solve(B_mat, b_vec);

  // colvec w2 = ones(nr.size());

  // colvec theta = ones(nr.size());
  // colvec var = ones(nr.size());
  // colvec m = ones(nr.size());

  //set y=1
      uvec pos_group = find(y_vec == 1);

  //set y=0
      uvec zero_group = find(y_vec ==0);
      int n_y0 = zero_group.n_elem;

  //thresholding r_tilde
      double r_tilde = accu(y_vec) / (double) y_vec.size() * r_ratio;
      double log_r = log(r_tilde);

      if(r_tilde >1){
        stop("r_tilde>1");
      }

    //PG latent variable
      colvec w = ones(yr.size());

      mat X_tilde1 = X_mat.rows(pos_group);
      vec k1 = y_vec(pos_group) - 0.5;

        mat trace_beta(run, p);
        vec filter_count(run);

        vec eta;
      for (int i = 0; i < (burnin + run); ++i)
      {


        //compute  Xbeta
        vec Xbeta = X_mat * beta;

        //sample Monte Carlo weight, eta: keep zero_group_subsample, throw away other y=0.
        eta = c11r.draw_binomial( mc_draws, r_tilde, n_y0 );
        uvec zero_group_subsample = zero_group(find(eta>0));

        //sample PG: w

        w.fill(0);

        vec alpha1 = zeros(n);
        vec alpha2 = zeros(n);

        // w for y=1
        alpha1(pos_group).fill(1);
        alpha2(pos_group) = Xbeta(pos_group);

        w(pos_group) = rpg(alpha1(pos_group), alpha2(pos_group));

        // w for y=0 and eta>0

        alpha1(zero_group_subsample) = eta(find(eta>0))/(double)mc_draws;
        alpha2(zero_group_subsample) = Xbeta(zero_group_subsample) - log_r;

        w(zero_group_subsample) = rpg(alpha1(zero_group_subsample), alpha2(zero_group_subsample));

        //sample beta

        mat X_tilde2 = X_mat.rows(zero_group_subsample);
        mat X_tilde = join_vert(X_tilde1, X_tilde2);

        vec w1 = w(pos_group);
        vec w2 = w(zero_group_subsample);
        vec w12 = join_vert(w1, w2);

        vec k2 = - eta(find(eta>0))/(double)mc_draws/2.0 + w2 * log_r;
        vec k12 = join_vert(k1, k2);

        mat wX_tilde = X_tilde;
        wX_tilde.each_col() %= w12;

        mat V = inv(X_tilde.t() * wX_tilde + B_inv);
        vec m = V*(X_tilde.t()*k12 + B_invb);

        beta = trans(chol(V)) * randn(p)  + m;

       
        if(i >= burnin){
            trace_beta.row(i - burnin) = beta.t();
            uvec filtered= find(eta==0);
            filter_count(i-burnin) = (double)filtered.n_elem / (double)n; 
          }
      }



      return Rcpp::List::create(
       Rcpp::Named("beta") = trace_beta,
       Rcpp::Named("eta") = eta,
       Rcpp::Named("filter_count") = filter_count
  //   Rcpp::Named("theta") = trace_theta,
  //   Rcpp::Named("w1") = trace_w1,
  //   Rcpp::Named("w2") = trace_w2
       );

    }

