#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double ParamLogLik(arma::vec param, arma::mat YT, arma::mat YNT, arma::mat riskNT, arma::mat riskT, 
                            arma::mat XNT, arma::mat XT, arma::mat XOR, arma::mat XinterMat, double epsOR)
{
  int n = YT.n_rows;
  int K = YT.n_cols;
  int pNT = XNT.n_cols;
  int pT = XT.n_cols;
  int pOR = XOR.n_cols;
  int pInter = XinterMat.n_cols;
  // Decomposing "param" to individual paramters 
//  double gamma = param[0];
  arma::vec alphaNT = param.subvec(0, K - 1);
  arma::vec alphaT = param.subvec(K, 2*K - 1);
  arma::vec alphaOR = param.subvec(2*K, 3*K - 1);
  arma::vec betaNT = param.subvec(3*K, 3*K + pNT - 1);
  arma::vec betaT = param.subvec(3*K + pNT, 3*K + pNT + pT - 1);
  arma::vec betaOR = param.subvec(3*K + pNT + pT, 3*K + pNT + pT + pOR - 1);
  arma::vec betaInt = param.subvec(3*K + pNT + pT + pOR, 3*K +  pNT + pT + pOR + pInter - 1);
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
  arma::vec ExpXBetaNT = exp(XNT * betaNT);
  arma::vec ExpXBetaT = exp(XT * betaT);
  arma::vec ExpXBetaOR = exp(XOR * betaOR);
  arma::vec ExpXBetaInt = exp(XinterMat * betaInt);
  arma::vec ExpAlphaNT = exp(alphaNT);
  arma::vec ExpAlphaT = exp(alphaT);
  arma::vec ExpAlphaOR = exp(alphaOR);
  for (int k = 0; k < K; ++k)
  {
    ExpAlphaNTnow = ExpAlphaNT[k];
    ExpAlphaTnow = ExpAlphaT[k];
    ExpAlphaORnow = ExpAlphaOR[k];
    for (int i = 0; i < n; ++i)
    {
      if (riskT(i,k)==0) {
        iContrib=0;
        //   nocontrib
      } else {
        if (riskNT(i,k)==0) {
          iProbTafterNT = (ExpAlphaTnow*ExpXBetaT[i]*ExpXBetaInt[i])/
            (1 + (ExpAlphaTnow*ExpXBetaT[i]*ExpXBetaInt[i]));
          if(YT(i,k)==1) {
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
          if (YNT(i,k)==1 && YT(i,k)==1) {
            iContrib = log(iProb12);
          }
          if (YNT(i,k)==1 && YT(i,k)==0) {
            iContrib = log(iProb1 - iProb12);
          }
          if (YNT(i,k)==0 && YT(i,k)==1) {
            iContrib = log(iProb2 - iProb12);
          }
          if (YNT(i,k)==0 && YT(i,k)==0) {
            iContrib = log(1 - iProb1 - iProb2 + iProb12);
          }
        }}
      loglik += iContrib;
    }}
  return(-loglik);
}
