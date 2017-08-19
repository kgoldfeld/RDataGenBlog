#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]

NumericMatrix groupX(IntegerVector Group1, IntegerVector Group2) {
  
  int n = Group1.size();
  
  ivec Group(Group1.begin(), n, false);
  ivec Group_unique = unique(Group);
  
  ivec GroupX(Group2.begin(), n, false);
  ivec Group_uniqueX = unique(GroupX);
  
  int k1 = Group_unique.size();
  int k2 = Group_uniqueX.size();
  
  NumericMatrix result(k1, k2);
  
  for(int i = 0; i < k1; i++) {
    for(int j = 0; j < k2; j++) {
      
      uvec test = find(Group == Group_unique(i) && GroupX == Group_uniqueX(j));
      result(i, j) = test.size(); 
      
    }
  }
  
  return result;
  
}

// [[Rcpp::export]]

bool cppChk( IntegerVector dtrow, 
             IntegerMatrix dorig) {
  
  int balanced = TRUE;
  int varnum = 0;
  
  while (balanced==TRUE && varnum <= dorig.ncol()-1 ) {
    
    IntegerVector colX = dorig(_, varnum);
    NumericMatrix X = groupX(colX, dtrow );
    LogicalVector Y = ( abs(X(_, 1) - (X(_, 0))) > 1);
    if (sum(Y) != 0) balanced = FALSE;
    varnum += 1;
    
  }
  
  return(balanced);
}