#include <RcppArmadillo.h>
#include <iostream>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
List cpp_pro_capsur(int i, int j, NumericVector beta, NumericMatrix ch, NumericMatrix cap_X, NumericMatrix surv_X) {
  int p = beta.size();
  NumericVector cap_beta = beta[Range(0, p / 2 - 1)];
  NumericVector surv_beta = beta[Range(p / 2, p - 1)];
  
 
  
  double zp = exp(cap_beta[0] + cap_beta[1] * cap_X(i-1, j-1));
  double zs = exp(surv_beta[0] + surv_beta[1] * surv_X(i-1, j-1));
  
  double p_hat = zp / (1 + zp);
  double s_hat = zs / (1 + zs);
  
  return List::create(Named("p.hat") = p_hat, Named("s.hat") = s_hat);
}

// [[Rcpp::export]]
List cpp_location(int nan, int ns, NumericMatrix ch) {
  NumericVector first(nan, 0.0);
  NumericVector last(nan, 0.0);

  for (int i = 0; i < nan; ++i) {
    bool findic = true;
    for (int j = 0; j < ns; ++j) {
      if (ch(i, j) >= 1) {
        if (findic && j < (ns - 1)) {
          first(i) = static_cast<double>(j + 2);
          findic = false;
        }
        last(i) = static_cast<double>(j + 1);
      }
    }
  }
  return List::create(Named("first") = first, Named("last") = last);
}

// [[Rcpp::export]]
double safe_log(double x) {
  return log(std::max(x, 1e-10));
}

// [[Rcpp::export]]
double cpp_CJSloglik_testing(NumericVector beta, NumericMatrix ch, NumericMatrix cap_X, NumericMatrix surv_X) {
  double xlnlik = 0;
  int nan = ch.nrow();
  int ns = ch.ncol();
  
  List loc = cpp_location(nan, ns, ch);
  NumericVector first = loc["first"];
  NumericVector last = loc["last"];

  for (int i = 0; i < nan; ++i) {
    double sum1 = 0.0;
    double sum2 = 0.0;
    NumericVector vp_ij(ns, 0.0);
    NumericVector vs_ij(ns, 0.0);
    
    int init1 = (first[i-1] > 0) ? first[i-1] : ns + 1;
    int init2 = (first[i-1] > 0) ? first[i-1] - 1 : ns + 1;

    if (init1 <= ns) {
      for (int j = init1; j <= ns; ++j) {
        vp_ij[j - 1] = cpp_pro_capsur(i, j, beta, ch, cap_X, surv_X)["p.hat"];
      }
    }

    if (init2 < ns) {
      for (int j = init2; j < ns; ++j) {
        vs_ij[j-1] = cpp_pro_capsur(i, j, beta, ch, cap_X, surv_X)["s.hat"];
      }
    }
    //compute log-likelihood contribution for carcass i
    if (first[i-1] > 0 && first[i-1] <= last[i-1]) {
            int start = first[i-1];
            int end = last[i-1];
        for (int j = start; j <= end; ++j) {
            double hij = (ch(i-1, j-1) >= 1) ? 1 : 0;
            sum1 += hij * safe_log(vp_ij[j-1]) + (1 - hij) * safe_log(1 - vp_ij[j-1]) + safe_log(vs_ij[j - 2]);
            }

        }

    if (last[i-1] > 0 && last[i-1] <= ns) {
      if (ch(i-1, last[i-1] - 1) >= 2) {
        sum2 = 0.0;
      } else {
        sum2 = 1 - vs_ij[last[i-1] - 1];
        for (int ii = last[i-1]; ii < ns; ++ii) {
          double prod = 1.0;
          for (int jj = last[i-1] - 1; jj < ii; ++jj) {
            prod *= vs_ij[jj] * (1 - vp_ij[jj + 1]);
          }
          if (ii < ns - 1) {
            prod *= (1 - vs_ij[ii]);
          }
          sum2 += prod;
        }
        sum2 = safe_log(sum2);
      }
    }
    xlnlik += sum1 + sum2;
        // Debugging output
    //Rcpp::Rcout << "i: " << i << ", sum1: " << sum1 << ", sum2: " << sum2 << ", xlnlik: " << xlnlik << std::endl;
  }
  return xlnlik;
}

// [[Rcpp::export]]
double cpp_CJSwrapper_testing(NumericVector beta, NumericMatrix ch, NumericMatrix cap_X, NumericMatrix surv_X) {
  double lnlik = -1 * cpp_CJSloglik_testing(beta, ch, cap_X, surv_X);
  return lnlik;
}

// [[Rcpp::export]]
List cpp_optim_testing(NumericVector beta, NumericMatrix ch, NumericMatrix cap_X, NumericMatrix surv_X) {
  Rcpp::Environment stats("package:stats");
  Rcpp::Function optim = stats["optim"];
  
  Rcpp::List result = optim(Rcpp::_["par"] = beta,
                            Rcpp::_["fn"] = Rcpp::InternalFunction(&cpp_CJSwrapper_testing),
                            Rcpp::_["method"] = "BFGS",
                            Rcpp::_["ch"] = ch,
                            Rcpp::_["cap_X"] = cap_X,
                            Rcpp::_["surv_X"] = surv_X,
                            Rcpp::_["control"] = List::create(Named("maxit") = 100),
                            Rcpp::_["hessian"] = true);
  
  return result;
}