#include <Rcpp.h>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
int FirstOne(NumericVector x)
  {
  int index = 0;
  int length = x.size();
  for (int i = 0; i < length; i++) {
    if(x[i]==1) {
      index = i + 1;
      break;
    }}
  return index;
}
