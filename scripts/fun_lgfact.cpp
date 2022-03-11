#include "RcppArmadillo.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::vec lgfactorial(int n, double sigma){
  
  arma::vec LogCk1(n+1); LogCk1.zeros();
  arma::vec LogCk(n+1);  LogCk.zeros();
  
  LogCk1(1) = std::log(sigma);
  LogCk1(0) = -arma::datum::inf;
  LogCk(0)  = -arma::datum::inf;
  
  for(int i = 2; i <= n; i++){
    LogCk(1)   = std::log(i - 1 - sigma) + LogCk1(1);
    for(int j = 2; j < i; j++){
      LogCk(j) = LogCk1(j-1)+ std::log(sigma + (i - 1 - sigma*j)*std::exp(LogCk1(j) - LogCk1(j-1)));
    }
    LogCk(i)  = i*std::log(sigma);
    LogCk1    = LogCk;
  }
  return(LogCk.rows(1,n));
}

// This function was created by Tommaso Rigon and was taken from 
// https://github.com/danieledurante/ESBM/blob/master/Source/stirling.cpp