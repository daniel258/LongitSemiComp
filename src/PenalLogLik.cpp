#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double PenalLogLik(arma::vec param, arma::mat YT, arma::mat YNT, arma::mat riskNT, arma::mat riskT, 
                  arma::mat XNT, arma::mat XT, arma::mat XOR, arma::mat XinterMat,
                  arma::mat TimeBase, arma::mat TimePen, arma::vec lambda, double epsOR)
{
  int n = YT.n_rows;
  int K = YT.n_cols;
  int pNT = XNT.n_cols;
  int pT = XT.n_cols;
  int pOR = XOR.n_cols;
  int pInter = XinterMat.n_cols;
  int J = TimeBase.n_cols; // J is the the number of B-splines (number of rows in TimeBase should be K)
  // Decomposing "param" to individual paramters 
  arma::vec alphaNT = param.subvec(0, J - 1);
  arma::vec alphaT = param.subvec(J, 2*J - 1);
  arma::vec alphaOR = param.subvec(2*J, 3*J - 1);
  arma::mat penaltermNT = lambda[0] * alphaNT.t() * TimePen * alphaNT;
  arma::mat penaltermT = lambda[1] * alphaT.t() * TimePen * alphaT;
  arma::mat penaltermOR = lambda[2] * alphaOR.t() * TimePen * alphaOR;
  arma::vec betaNT = param.subvec(3*J, 3*J + pNT - 1);
  arma::vec betaT = param.subvec(3*J + pNT, 3*J + pNT + pT - 1);
  arma::vec betaOR = param.subvec(3*J + pNT + pT,3*J +  pNT + pT + pOR - 1);
  arma::vec betaInt = param.subvec(3*J + pNT + pT + pOR, 3*J +  pNT + pT + pOR + pInter - 1);
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
  arma::vec ExpAlphaNT = exp(TimeBase * alphaNT);
  arma::vec ExpAlphaT = exp(TimeBase * alphaT);
  arma::vec ExpAlphaOR = exp(TimeBase * alphaOR);
  for (int j = 0; j < K; ++j)
  {
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
          iProbTafterNT = (ExpAlphaTnow*ExpXBetaT[i]*ExpXBetaInt[i])/
            (1 + (ExpAlphaTnow*ExpXBetaT[i]*ExpXBetaInt[i]));
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
  double  penalloglik = loglik - as_scalar(penaltermNT) - as_scalar(penaltermT) - as_scalar(penaltermOR);
  return(-penalloglik);
}
