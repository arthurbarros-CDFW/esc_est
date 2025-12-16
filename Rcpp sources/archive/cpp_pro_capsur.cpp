#include <Rcpp.h>

//[[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List pro_capsur(int i, int j, Rcpp::NumericMatrix ch,
                      Rcpp::NumericVector beta, Rcpp::NumericVector cap_X,
                      Rcpp::NumericVector surv_X) {
    int nan = ch.nrow();
    int ns = ch.ncol();
    
    Rcpp::NumericMatrix cap_matrix(nan, ns);
    Rcpp::NumericMatrix surv_matrix(nan, ns);
    
    for (int row = 0; row < nan; ++row) {
        for (int col = 0; col < ns; ++col) {
            cap_matrix(row, col) = cap_X[row];
            surv_matrix(row, col) = surv_X[row];
        }
    }
    
    int p = beta.size();
    Rcpp::NumericVector cap_beta = beta[Rcpp::Range(0, p / 2 - 1)];
    Rcpp::NumericVector surv_beta = beta[Rcpp::Range(p / 2, p - 1)];
    
    double zp = std::exp(cap_beta[0] * 1 + cap_beta[1] * cap_matrix(i, j));
    double zs = std::exp(surv_beta[0] * 1 + surv_beta[1] * surv_matrix(i, j));
    
    double p_hat = zp / (1 + zp);
    double s_hat = zs / (1 + zs);
    
    return Rcpp::List::create(Rcpp::Named("p.hat") = p_hat,
                              Rcpp::Named("s.hat") = s_hat);
}