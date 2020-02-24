#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double ParamLogLik(arma::vec param, arma::mat X, arma::mat YNT, arma::mat YT, 
                   arma::mat riskNT, arma::mat riskT, double epsOR)
     {
  int n = YT.n_rows;
  int J = YT.n_cols;
  int p = X.n_cols;
  double gamma = param[0];
  arma::vec alphaNT = param.subvec(1,J);
  arma::vec alphaT = param.subvec(J+1,2*J);
  arma::vec alphaOR = param.subvec(2*J+1,3*J);
  arma::vec betaNT = param.subvec(3*J+1,3*J+p);
  arma::vec betaT = param.subvec(3*J+p+1,3*J+2*p);
  arma::vec betaOR = param.subvec(3*J+2*p+1,3*J+3*p);
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
  arma::vec ExpXBetaOR = exp(X * betaOR);
  arma::vec ExpAlphaNT = exp(alphaNT);
  arma::vec ExpAlphaT = exp(alphaT);
  arma::vec ExpAlphaOR = exp(alphaOR);
  for (int j = 0; j < J; ++j)
  {
    ExpAlphaNTnow = ExpAlphaNT[j];
    ExpAlphaTnow = ExpAlphaT[j];
    ExpAlphaORnow = ExpAlphaOR[j];
    for (int i = 0; i < n; ++i)
    {
        if (riskT(i,j)==0) {
          iContrib=0;          //   nocontrib
        } else {
          if (riskNT(i,j)==0) {
            iProbTafterNT = (ExpAlphaTnow*ExpXBetaT[i]*exp(gamma))/(1 + (ExpAlphaTnow*ExpXBetaT[i]*exp(gamma)));
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
            if ((iOR > 1 - epsOR) & (iOR < 1 + epsOR))
            {
              iProb12 = iProb1*iProb2;
            } else {
              iProb12 = (1 + (iProb1 + iProb2)*(iOR - 1) - sqrt(pow(1 + (iProb1 + iProb2)*(iOR - 1), 2.0) -
                4*iOR*(iOR - 1)*iProb1*iProb2)) / (2 * (iOR - 1));
            }
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
      loglik += iContrib;
    }}
  return(-loglik);
     }
