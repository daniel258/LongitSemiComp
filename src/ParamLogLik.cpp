#include <RcppArmadillo.h>
using namespace Rcpp;
//using namespace arma;
// arma::vec ExpitModel(double alpha, arma::vec beta, arma::mat X, bool expitInd)
// {
//   arma::vec lin = alpha + X * beta;
//   arma::vec ret = arma::zeros(X.n_rows);
//   if(expitInd == true)  {
//      ret = exp(lin)/(1 + exp(lin));
//    }  else {
//   ret = exp(lin);
//   };
//   return ret;
// }

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double ParamLogLik(arma::vec param, arma::mat X, arma::mat YNT, arma::mat YT)
     {
  int n = YT.n_rows;
  int J = YT.n_cols;
  int p = X.n_cols;
  // Rcpp::Rcout << "p= "<< p << std::endl;
// Decomposing "param" to individual paramters as in LogLik
  double gamma = param[0];
  // Rcpp::Rcout << "gamma "<< gamma << std::endl;
  // double alpha0NT = param[1];
  // double alpha0T = param[2];
  // double alpha0OR = param[3];
  arma::vec alphaNT = param.subvec(1,J);
  arma::vec alphaT = param.subvec(J+1,2*J);
  arma::vec alphaOR = param.subvec(2*J+1,3*J);
  // Rcpp::Rcout << "alphaNT "<< alphaNT << std::endl;
  // Rcpp::Rcout << "alphaT "<< alphaT << std::endl;
  // Rcpp::Rcout << "alphaOR "<< alphaOR << std::endl;
  arma::vec betaNT = param.subvec(3*J+1,3*J+p);
  // Rcpp::Rcout << "betaNT "<< betaNT << std::endl;
  arma::vec betaT = param.subvec(3*J+p+1,3*J+2*p);
  // Rcpp::Rcout << "betaT "<< betaT << std::endl;
  arma::vec betaOR = param.subvec(3*J+2*p+1,3*J+3*p);
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
//  LogicalVector NTstatus(n,true);
//  LogicalVector Tstatus(n,true);
//  arma::vec YNTnow = arma::zeros(n);
//  arma::vec YTnow = arma::zeros(n);
  // double Expalpha0NT = exp();
  // double Expalpha0T = exp();
  // double Expalpha0OR = exp();
  arma::vec ExpXBetaNT = exp(X * betaNT);
  arma::vec ExpXBetaT = exp(X * betaT);
  // Rcpp::Rcout << (ExpXBetaT) << std::endl;
  // Rcpp::Rcout << (ExpXBetaNT) << std::endl;
  arma::vec ExpXBetaOR = exp(X * betaOR);
  arma::vec ExpAlphaNT = exp(alphaNT);
  arma::vec ExpAlphaT = exp(alphaT);
  arma::vec ExpAlphaOR = exp(alphaOR);
  //   arma::mat contrib=1 + ps*(exp(beta)-1);
  for (int j = 0; j < J; ++j)
  {
    //YNTnow = YNT.submat(0, j, n-1, j);
    //  YTnow = YT.submat(0, j, n-1, j);
    ExpAlphaNTnow = ExpAlphaNT[j];
    ExpAlphaTnow = ExpAlphaT[j];
    ExpAlphaORnow = ExpAlphaOR[j];
    for (int i = 0; i < n; ++i)
    {
      if (j==0) {
        iProb1 = ExpAlphaNTnow*ExpXBetaNT[i]/(1 + (ExpAlphaNTnow*ExpXBetaNT[i]));
        iProb2 = ExpAlphaTnow*ExpXBetaT[i]/(1 + ExpAlphaTnow*ExpXBetaT[i]);
        iOR = ExpAlphaORnow*ExpXBetaOR[i];
        // Rcpp::Rcout << "iProb1  "<< iProb1 << std::endl;
        // Rcpp::Rcout << "iProb2 "<< iProb2 << std::endl;
        // Rcpp::Rcout << "iOR "<< iOR << std::endl;
        if (iOR==1)
        {
          iProb12 = iProb1*iProb2;
        } else {
          iProb12 = (1 + (iProb1 + iProb2)*(iOR - 1) - sqrt(pow(1 + (iProb1 + iProb2)*(iOR - 1), 2.0) -
            4*iOR*(iOR - 1)*iProb1*iProb2)) / (2 * (iOR - 1));
        }
        // Rcpp::Rcout << "iProb12  "<< iProb12 << std::endl;
        if (YNT(i,0)==1 && YT(i,0)==1) {
          iContrib = log(iProb12);
        }
        if (YNT(i,0)==1 && YT(i,0)==0) {
          iContrib = log(iProb1 - iProb12);
        }
        if (YNT(i,0)==0 && YT(i,0)==1) {
          iContrib = log(iProb2 - iProb12);
        }
        if (YNT(i,0)==0 && YT(i,0)==0) {
          iContrib = log(1 - iProb1 - iProb2 + iProb12);
        }
      } else {
        if (YT(i,j-1)==1) {
          iContrib=0;
          //   nocontrib
        } else {
          if (YNT(i,j-1)==1) {
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
            if (iOR==1)
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
          }}}
      // Rcpp::Rcout << "iContrib "<< iContrib << std::endl;
      loglik += iContrib;
    }}
  return(-loglik);
     }
