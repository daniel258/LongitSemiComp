#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
//' @export
// [[Rcpp::export]]
double PenalLogLikORconstNoX(arma::vec param, arma::mat YNT, arma::mat riskNT, arma::mat riskT, arma::mat YT, arma::mat TimeBase, arma::mat TimePen, arma::vec lambda)
     {
  int n = YT.n_rows;
  int J = YT.n_cols;
  int Q = TimeBase.n_cols; // Q is the the number of B-splines (number of rows in TimeBase should be J)
  // Rcpp::Rcout << "p= "<< p << std::endl;
// Decomposing "param" to individual paramters as in LogLik
  double gamma = param[0];
  // Rcpp::Rcout << "gamma "<< gamma << std::endl;
  // double alpha0NT = param[1];
  // double alpha0T = param[2];
  // double alpha0OR = param[3];
 // arma::vec alphaNT = param.subvec(1,J);
//  arma::vec alphaT = param.subvec(J+1,2*J);
//  arma::vec alphaOR = param.subvec(2*J+1,3*J);
  arma::vec alphaNT = param.subvec(1,Q);
  arma::vec alphaT = param.subvec(Q+1,2*Q);
  double alphaOR = param[2*Q+1];
  arma::mat penaltermNT = lambda[0] * alphaNT.t() * TimePen * alphaNT;
  arma::mat penaltermT = lambda[1] * alphaT.t() * TimePen * alphaT;
  // Rcpp::Rcout << "penaltermNT "<< penaltermNT << std::endl;
  // Rcpp::Rcout << "penaltermT "<< penaltermT << std::endl;
  // Rcpp::Rcout << "penaltermOR "<< penaltermOR << std::endl;
  // Rcpp::Rcout << "alphaNT "<< alphaNT << std::endl;
  // Rcpp::Rcout << "alphaT "<< alphaT << std::endl;
  // Rcpp::Rcout << "alphaOR "<< alphaOR << std::endl;
  double loglik = 0;
  double iContrib = 0;
  double iProbTafterNT = 0;
  double iProb1 = 0;
  double iProb2 = 0;
  double iOR = 0;
  double iProb12 = 0;
  double ExpAlphaNTnow = 0;
  double ExpAlphaTnow = 0;
  arma::vec ExpAlphaNT = exp(TimeBase * alphaNT);
  arma::vec ExpAlphaT = exp(TimeBase * alphaT);
  double ExpAlphaOR = exp(alphaOR);
  for (int j = 0; j < J; ++j)
  {
    ExpAlphaNTnow = ExpAlphaNT[j];
    ExpAlphaTnow = ExpAlphaT[j];
    for (int i = 0; i < n; ++i)
    {
        if (riskT(i,j)==0) {
          iContrib=0;
          //   nocontrib
        } else {
          if (riskNT(i,j)==0) {
            iProbTafterNT = (ExpAlphaTnow*exp(gamma))/(1 + (ExpAlphaTnow*exp(gamma)));
            if(YT(i,j)==1) {
              iContrib = log(iProbTafterNT);
            }
            else {
              iContrib = log(1-iProbTafterNT);
            }
          } else {
            iProb1 = ExpAlphaNTnow/(1 + ExpAlphaNTnow);
            iProb2 = ExpAlphaTnow/(1 + ExpAlphaTnow);
            iOR = ExpAlphaOR;
            if ((iOR>0.99) & (iOR<1.01))
            {
              iProb12 = iProb1*iProb2;
            } else {
              iProb12 = (1 + (iProb1 + iProb2)*(iOR - 1) - sqrt(pow(1 + (iProb1 + iProb2)*(iOR - 1), 2.0) -
                4*iOR*(iOR - 1)*iProb1*iProb2)) / (2 * (iOR - 1));
            }
            // Rcpp::Rcout << "iProb12  "<< iProb12 << std::endl;
            if (YNT(i,j)==1 && YT(i,j)==1) {
              iContrib = log(iProb12);
            }
            if (YNT(i,j)==1 && YT(i,j)==0) {
              iContrib = log(iProb1 - iProb12);
            }
            if (YNT(i,j)==0 && YT(i,j)==1) {
              iContrib = log(iProb2 - iProb12);
            }
            if (YNT(i,j)==0 && YT(i,j)==0) {
              iContrib = log(1 - iProb1 - iProb2 + iProb12);
            }
          }}
      // Rcpp::Rcout << "iContrib "<< iContrib << std::endl;
      loglik += iContrib;
    }}
double  penalloglik = loglik - as_scalar(penaltermNT) - as_scalar(penaltermT);
 // Rcpp::Rcout << "penalloglik "<< penalloglik << std::endl;
  return(-penalloglik);
     }
