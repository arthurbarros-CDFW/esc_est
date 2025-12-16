#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List location(int nan, int ns, NumericMatrix ch) {
  NumericVector first(nan, 0);
  NumericVector last(nan, 0);
  
  for (int i = 0; i < nan; ++i) {
    bool findch = true;
    for (int j = 0; j < ns; ++j) {
      if (ch(i, j) >= 1) {
        if (findch && (j < ns)) {
          first[i] = j + 1;
          findch = false;
        }
        last[i] = j;
      }
    }
  }
  
  return List::create(Named("first") = first, Named("last") = last);
}