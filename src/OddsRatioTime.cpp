#include <Rcpp.h>
using namespace Rcpp;


//// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List OddsRatioTime(NumericMatrix YNT, NumericMatrix YT, NumericMatrix riskNT, NumericMatrix riskT)
     {
  int n = YT.nrow();
  int J = YT.ncol();
  NumericVector N(J);
  NumericVector OR(J);
  NumericVector n11(J);
  NumericVector n10(J);
  NumericVector n01(J);
  NumericVector n00(J);
  // Rcpp::Rcout << (n01) << std::endl;
  // Rcpp::Rcout << (ExpXBetaNT) << std::endl;

  for (int j = 0; j < J; ++j)
   {
    for (int i = 0; i < n; ++i)
     {
  //    Rcpp::Rcout << (YNT(i,j)) << std::endl;
       if (riskNT(i,j)==1 && riskT(i,j)==1) {
         // Rcpp::Rcout << "YNT(i,j)  "<< YNT(i,j) << std::endl;
         // Rcpp::Rcout << "YT(i,j)  "<< YT(i,j) << std::endl;
         N[j] += 1;
           if (YNT(i,j)==1 && YT(i,j)==1) {
              n11[j] += 1;
             }
           if (YNT(i,j)==0 && YT(i,j)==1) {
             n01[j] += 1;
             }
           if (YNT(i,j)==1 && YT(i,j)==0) {
             n10[j] += 1;
           }
           if (YNT(i,j)==0 && YT(i,j)==0) {
             n00[j] += 1;
           }}
     }
//    Rcpp::Rcout << (n01[j]) << std::endl;
    OR[j] = n11[j]*n00[j]/(n01[j]*n10[j]);
   }
  List out = List::create(Named("OddsRatios") = OR , Named("N")  = N, Named("n00") = n00, Named("n01") = n01, Named("n10") = n10, Named("n11") = n11);
  return(out);
     }
