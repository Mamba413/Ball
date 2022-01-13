// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// sbd_cpp
double sbd_cpp(Rcpp::NumericMatrix& x, int n1, int n1_total, int n2, int n2_total);
RcppExport SEXP _Ball_sbd_cpp(SEXP xSEXP, SEXP n1SEXP, SEXP n1_totalSEXP, SEXP n2SEXP, SEXP n2_totalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n1_total(n1_totalSEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< int >::type n2_total(n2_totalSEXP);
    rcpp_result_gen = Rcpp::wrap(sbd_cpp(x, n1, n1_total, n2, n2_total));
    return rcpp_result_gen;
END_RCPP
}