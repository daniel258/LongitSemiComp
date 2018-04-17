#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double PenalLogLik(arma::vec param, arma::mat X, arma::mat YNT, arma::mat riskNT, arma::mat riskT, arma::mat YT, arma::mat TimeBase, arma::mat TimePen, arma::vec lambda)
     {
  int n = YT.n_rows;
  int J = YT.n_cols;
  int p = X.n_cols;
  int Q = TimeBase.n_cols; // Q is the the number of B-splines (number of rows in TimeBase should be J)
  // Rcpp::Rcout << "p= "<< p << std::endl;
// Decomposing "param" to individual paramters as in LogLik
// Rcpp::Rcout << "param: "<< param << std::endl;
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
  arma::vec alphaOR = param.subvec(2*Q+1,3*Q);
  arma::mat penaltermNT = lambda[0] * alphaNT.t() * TimePen * alphaNT;
  arma::mat penaltermT = lambda[1] * alphaT.t() * TimePen * alphaT;
  arma::mat penaltermOR = lambda[2] * alphaOR.t() * TimePen * alphaOR;
  // Rcpp::Rcout << "penaltermNT "<< penaltermNT << std::endl;
  // Rcpp::Rcout << "penaltermT "<< penaltermT << std::endl;
  // Rcpp::Rcout << "penaltermOR "<< penaltermOR << std::endl;
  // Rcpp::Rcout << "alphaNT "<< alphaNT << std::endl;
  // Rcpp::Rcout << "alphaT "<< alphaT << std::endl;
  // Rcpp::Rcout << "alphaOR "<< alphaOR << std::endl;
  arma::vec betaNT = param.subvec(3*Q+1,3*Q+p);
  // Rcpp::Rcout << "betaNT "<< betaNT << std::endl;
  arma::vec betaT = param.subvec(3*Q+p+1,3*Q+2*p);
  // Rcpp::Rcout << "betaT "<< betaT << std::endl;
  arma::vec betaOR = param.subvec(3*Q+2*p+1,3*Q+3*p);
  // Rcpp::Rcout << "betaOR "<< betaOR << std::endl;
  double loglik = 0;
  double iContrib = 0;
  double iProbTafterNT = 0;
  double iProb1 = 0;
  double iProb2 = 0;
  double iOR = 0;
  double iProb12 = 0;
  double ExpAlphaNTnow = 0;
  double ExpAlphaTnow = 0;
  double ExpAlphaORnow = 0;
  arma::vec ExpXBetaNT = exp(X * betaNT);
  arma::vec ExpXBetaT = exp(X * betaT);
  // Rcpp::Rcout << (ExpXBetaT) << std::endl;
  // Rcpp::Rcout << (ExpXBetaNT) << std::endl;
  arma::vec ExpXBetaOR = exp(X * betaOR);
  arma::vec ExpAlphaNT = exp(TimeBase * alphaNT);
  arma::vec ExpAlphaT = exp(TimeBase * alphaT);
  arma::vec ExpAlphaOR = exp(TimeBase * alphaOR);
  for (int j = 0; j < J; ++j)
  {
    //YNTnow = YNT.submat(0, j, n-1, j);
    //  YTnow = YT.submat(0, j, n-1, j);
    ExpAlphaNTnow = ExpAlphaNT[j];
    ExpAlphaTnow = ExpAlphaT[j];
    ExpAlphaORnow = ExpAlphaOR[j];
    for (int i = 0; i < n; ++i)
    {
        if (riskT(i,j)==0) {
          iContrib=0;
          //   nocontrib
        } else {
          if (riskNT(i,j)==0) {
            iProbTafterNT = (ExpAlphaTnow*ExpXBetaT[i]*exp(gamma))/(1 + (ExpAlphaTnow*ExpXBetaT[i]*exp(gamma)));
            // Rcpp::Rcout << "ExpAlphaTnow "<< ExpAlphaTnow << std::endl;
            // Rcpp::Rcout << "ExpXBetaT[i] "<< ExpXBetaT[i] << std::endl;
            // Rcpp::Rcout << "exp(gamma) "<< exp(gamma) << std::endl;
            // Rcpp::Rcout << "iProbTafterNT "<< iProbTafterNT << std::endl;
            if(YT(i,j)==1) {
              iContrib = log(iProbTafterNT);
            }
            else {
              iContrib = log(1-iProbTafterNT);
            }
          } else {
            iProb1 = ExpAlphaNTnow*ExpXBetaNT[i]/(1 + (ExpAlphaNTnow*ExpXBetaNT[i]));
            iProb2 = ExpAlphaTnow*ExpXBetaT[i]/(1 + ExpAlphaTnow*ExpXBetaT[i]);
            iOR = ExpAlphaORnow*ExpXBetaOR[i];
            // Rcpp::Rcout << "iProb1  "<< iProb1 << std::endl;
            // Rcpp::Rcout << "iProb2 "<< iProb2 << std::endl;
            // Rcpp::Rcout << "iOR "<< iOR << std::endl;
            if (iOR>0.9999 & iOR<1.0001)
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
      // if (arma::is_finite(iContrib)) {} else {
      //   Rcpp::Rcout << "iContrib "<< iContrib << std::endl;
      //   Rcpp::Rcout << "j =  "<< j << std::endl;
      //   Rcpp::Rcout << "i =  "<< i << std::endl;
      //   Rcpp::Rcout << "iProb1  "<< iProb1 << std::endl;
      //   Rcpp::Rcout << "iProb2  "<< iProb2 << std::endl;
      //   Rcpp::Rcout << "iProb12  "<< iProb12 << std::endl;
      //   Rcpp::Rcout << "iOR  "<< iOR << std::endl;
      //   Rcpp::Rcout << "Risk-NT"<< riskNT(i,j) << std::endl;
      //   Rcpp::Rcout << "Risk-T  "<< riskT(i,j) << std::endl;
      //   Rcpp::Rcout << "YNT "<< YNT(i,j) << std::endl;
      //   Rcpp::Rcout << "YT  "<< YT(i,j) << std::endl;}
//        Rcpp::Rcout << "iORCOND  "<< iORcond << std::endl;}
      loglik += iContrib;
    }}
double  penalloglik = loglik - as_scalar(penaltermNT) - as_scalar(penaltermT) - as_scalar(penaltermOR);
 // Rcpp::Rcout << "penalloglik "<< penalloglik << std::endl;
  return(-penalloglik);
     }
