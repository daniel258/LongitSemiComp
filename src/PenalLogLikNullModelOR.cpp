#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double PenalLogLikNullModelOR(arma::vec param, arma::mat YT, arma::mat YNT, arma::mat riskNT, arma::mat riskT, 
                              arma::mat XNT, arma::mat XT, //arma::mat XOR,  
                              arma::mat TimeBase, arma::mat TimePen, 
                              arma::vec lambda, double epsOR)
{
  int n = YT.n_rows;
  int J = YT.n_cols;
  int pNT = XNT.n_cols;
  int pT = XT.n_cols;
  //int pOR = XOR.n_cols;
  //int pOR = 0;
  int Q = TimeBase.n_cols; // Q is the the number of B-splines (number of rows in TimeBase should be J)
  // Decomposing "param" to individual paramters 
  double gamma = param[0];
  arma::vec alphaNT = param.subvec(1,Q);
  arma::vec alphaT = param.subvec(Q+1,2*Q);
  arma::vec alphaOR = param.subvec(2*Q+1,3*Q);
  arma::mat penaltermNT = lambda[0] * alphaNT.t() * TimePen * alphaNT;
  arma::mat penaltermT = lambda[1] * alphaT.t() * TimePen * alphaT;
  arma::mat penaltermOR = lambda[2] * alphaOR.t() * TimePen * alphaOR;
  arma::vec betaNT = param.subvec(3*Q + 1, 3*Q + pNT);
  arma::vec betaT = param.subvec(3*Q + pNT + 1,3*Q + pNT + pT);
  //arma::vec betaOR = param.subvec(3*Q + pNT + pT + 1,3*Q +  pNT + pT + pOR);
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
  //arma::vec ExpXBetaOR = exp(XOR * betaOR);
  arma::vec ExpAlphaNT = exp(TimeBase * alphaNT);
  arma::vec ExpAlphaT = exp(TimeBase * alphaT);
  arma::vec ExpAlphaOR = exp(TimeBase * alphaOR);
  for (int j = 0; j < J; ++j)
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
          iOR = ExpAlphaORnow;
          // iOR = ExpAlphaORnow*ExpXBetaOR[i];
          // Rcpp::Rcout << "iProb1  "<< iProb1 << std::endl;
          // Rcpp::Rcout << "iProb2 "<< iProb2 << std::endl;
          // Rcpp::Rcout << "iOR "<< iOR << std::endl;
          if ((iOR > 1 - epsOR) & (iOR < 1 + epsOR))
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
      loglik += iContrib;
    }}
    double  penalloglik = loglik - as_scalar(penaltermNT) - as_scalar(penaltermT) - as_scalar(penaltermOR);
    return(-penalloglik);
}
