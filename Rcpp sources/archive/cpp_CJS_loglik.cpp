#include <Rcpp.h>
using namespace Rcpp;

// Function declaration for location and pro_capsur
List cpp_location(int nan, int ns, NumericMatrix ch);
List cpp_pro_capsur(int i, int j, NumericMatrix ch, NumericVector beta, NumericVector cap_X, NumericVector surv_X);

// [[Rcpp::export]]
double cpp_CJSloglik(NumericVector beta, NumericVector cap_X, NumericVector surv_X, NumericMatrix ch) {
  double xlnlik = 0;
  int nan = ch.nrow();
  int ns = ch.ncol();
  int p = beta.size();
  
  List loc = cpp_location(nan, ns, ch);
  IntegerVector first = loc["first"];
  IntegerVector last = loc["last"];
  
  for (int i = 0; i < nan; ++i) {
    double sum1 = 0;
    double sum2 = 0;
    NumericVector vp_ij(ns, NA_REAL);
    NumericVector vs_ij(ns, NA_REAL);
    
    int init_cap = (first[i] == 0) ? ns + 1 : first[i];
    int init_surv = (first[i] > 0) ? first[i] - 1 : ns + 1;
    
    if (init_cap <= ns) {
      for (int j = init_cap; j <= ns; ++j) {
        vp_ij[j - 1] = as<List>(cpp_pro_capsur(i, j, ch, beta, cap_X, surv_X))["p.hat"];
      }
    }
    if (init_surv < ns) {
      for (int j = init_surv; j < ns; ++j) {
        vs_ij[j] = as<List>(cpp_pro_capsur(i, j + 1, ch, beta, cap_X, surv_X))["s.hat"];
      }
    }
    
    if (first[i] > 0 && first[i] <= last[i]) {
      for (int j = first[i]; j <= last[i]; ++j) {
        int hij = (ch(i, j - 1) >= 1) ? 1 : 0;
        sum1 += hij * log(vp_ij[j - 1]) + (1 - hij) * log(1 - vp_ij[j - 1]) + log(vs_ij[j - 2]);
      }
    }
    
    if (ch(i, last[i] - 1) >= 2) {
      sum2 = 0;
    } else if (last[i] > 0 && last[i] < ns) {
      sum2 = 1 - vs_ij[last[i] - 1];
      for (int ii = last[i] + 1; ii <= ns; ++ii) {
        double prod = 1;
        for (int jj = last[i]; jj < ii; ++jj) {
          prod *= vs_ij[jj - 1] * (1 - vp_ij[jj]);
        }
        if (ii < ns) {
          prod *= (1 - vs_ij[ii - 1]);
        }
        sum2 += prod;
      }
      sum2 = log(sum2);
    }
    
    xlnlik += sum1 + sum2;
  }
  
  return xlnlik;
}