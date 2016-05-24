// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// rpg
SEXP rpg(SEXP b, SEXP c);
RcppExport SEXP ImbalancedPG_rpg(SEXP bSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type b(bSEXP);
    Rcpp::traits::input_parameter< SEXP >::type c(cSEXP);
    __result = Rcpp::wrap(rpg(b, c));
    return __result;
END_RCPP
}
// multinomial
SEXP multinomial(SEXP n, SEXP y, double b, double B, int burnin, int run, double r_ratio);
RcppExport SEXP ImbalancedPG_multinomial(SEXP nSEXP, SEXP ySEXP, SEXP bSEXP, SEXP BSEXP, SEXP burninSEXP, SEXP runSEXP, SEXP r_ratioSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type n(nSEXP);
    Rcpp::traits::input_parameter< SEXP >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type run(runSEXP);
    Rcpp::traits::input_parameter< double >::type r_ratio(r_ratioSEXP);
    __result = Rcpp::wrap(multinomial(n, y, b, B, burnin, run, r_ratio));
    return __result;
END_RCPP
}