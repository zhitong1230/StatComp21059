#include <Rcpp.h>
using namespace Rcpp;

//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp
//' @param N the number of samples
//' @param n the number of between-sample random numbers
//' @param a the first parameter of Beta distribution
//' @param b the second parameter of Beta distribution
//' @return a random sample of size \code{n}
//' @examples
//' \dontrun{
//' set.seed(1230)
//' NumericMatrix gibbsC(1000,1000,1,1)
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix gibbsC(int N, int n,int a,int b) {
  NumericMatrix mat(N, 2);
  double x = 0, y = 0;
  for(int i = 1; i < N; i++) {
    y=mat(i-1,1);
    mat(i,0)=rbinom(1, n, y)[0];
    x=mat(i,0);
    mat(i,1)=rbeta(1, x+a, n-x+b)[0];
  }
  return(mat);
}