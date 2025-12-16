#include <Rcpp.h>
using namespace Rcpp;

// Function declaration for CJS_loglik
double CJS_loglik(NumericVector beta, NumericVector cap_X, NumericVector surv_X, NumericMatrix ch);

// [[Rcpp::export]]
double CJS_loglik_wrapper(NumericVector beta, NumericVector cap_X, NumericVector surv_X, NumericMatrix ch) {
  return -CJS_loglik(beta, cap_X, surv_X, ch);  // Return the negative log-likelihood
}