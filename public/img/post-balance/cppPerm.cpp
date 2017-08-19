#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]

NumericVector stl_sort(NumericVector x) {
  
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y;
  
}

// [[Rcpp::export]]

NumericMatrix cppPerm(int ncol, int nchoose, int N) {
  
  NumericVector X(ncol);
  NumericVector Xtemp(ncol);
  NumericMatrix Y(N, nchoose);
  
  for (int j = 0; j < ncol; j++) {
    X(j) = j + 1;
  }
  
  for (int i = 0; i < N; i++) {
    Xtemp = RcppArmadillo::sample(X, nchoose, FALSE);
    Y(i, _) = stl_sort(Xtemp);
  }
  
  return(Y);
}
