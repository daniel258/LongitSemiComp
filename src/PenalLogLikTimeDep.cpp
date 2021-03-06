#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double PenalLogLikTimeDep(arma::vec param, arma::vec ID, arma::uvec TM,
                      arma::vec YT, arma::vec YNT, 
                   arma::mat XNT, arma::mat XT, arma::mat XOR,
                   arma::mat TimeBase, arma::mat TimePen, arma::vec lambda, double epsOR)
{
  arma::vec IDunq = arma::unique(ID); // IDs once
  int n = IDunq.size(); // sample size
  arma::uvec TMunq = arma::unique(TM); // unique "times" (intervals)
 // int J = TMunq.size(); // number of "times" (intervals)
  int pNT = XNT.n_cols;
  int pT = XT.n_cols;
  int pOR = XOR.n_cols;
  int Q = TimeBase.n_cols; // Q is the the number of B-splines (number of rows in TimeBase should be J)
  // Decomposing "param" to individual paramters
  // Current verison of the code: only one time-dep variable (but throgout the code, I made preperations for generalizing it)
  //int pNTtimeDep = 1;
  //int pTtimeDep = 1;
  //int pORtimeDep = 1;
  double betay = param[0];
  arma::vec alphaNT = param.subvec(1,Q);
  arma::vec alphaT = param.subvec(Q+1,2*Q);
  arma::vec alphaOR = param.subvec(2*Q+1,3*Q);
  arma::mat penaltermNT = lambda[0] * alphaNT.t() * TimePen * alphaNT;
  arma::mat penaltermT = lambda[1] * alphaT.t() * TimePen * alphaT;
  arma::mat penaltermOR = lambda[2] * alphaOR.t() * TimePen * alphaOR;
  arma::vec betaNT = param.subvec(3*Q + 1, 3*Q + pNT);
  arma::vec betaT = param.subvec(3*Q + pNT + 1, 3*Q + pNT + pT);
  arma::vec betaOR = param.subvec(3*Q + pNT + pT + 1, 3*Q + pNT + pT + pOR);
  // Rcpp::Rcout << "betaNT:  "<< betaNT << std::endl;
  // Rcpp::Rcout << "betaT:  "<< betaT << std::endl;
  // Rcpp::Rcout << "betaOR:  "<< betaOR << std::endl;
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
  double ExpXBetaNTnow = 0;
  double ExpXBetaTnow = 0;
  double ExpXBetaORnow = 0;
  arma::vec ExpAlphaNT = exp(TimeBase * alphaNT);
  arma::vec ExpAlphaT = exp(TimeBase * alphaT);
  arma::vec ExpAlphaOR = exp(TimeBase * alphaOR);
  arma::vec ExpXBetaNT = exp(XNT * betaNT);
  arma::vec ExpXBetaT = exp(XT * betaT);
  arma::vec ExpXBetaOR = exp(XOR * betaOR);
  for (int i = 0; i < n; ++i)
  {
    int iID = IDunq[i]; // get observation ID
    // Get all data on that observation
    arma::uvec pos = find(ID==iID); // rows associated with the observation
    arma::vec iYNT = YNT.elem(pos); // YNT for that observation
    arma::vec iYT = YT.elem(pos); // YT for that observation
    arma::uvec iTM = TM.elem(pos); // TM for that observation, they are already unique
    arma::vec iExpXBetaNT = ExpXBetaNT.elem(pos); // exp(X%*%\beta) for NT, for that observation (all times)
    arma::vec iExpXBetaT = ExpXBetaT.elem(pos); // exp(X%*%\beta) for T, for that observation (all times)
    arma::vec iExpXBetaOR = ExpXBetaOR.elem(pos); // exp(X%*%\beta) for OR, for that observation (all times)
    arma::vec iExpAlphaNT = ExpAlphaNT.elem(iTM - 1);
    arma::vec iExpAlphaT = ExpAlphaT.elem(iTM - 1);
    arma::vec iExpAlphaOR = ExpAlphaOR.elem(iTM - 1);
    int iJ = iTM.size(); // number of "times" (intervals)
    arma::vec iRiskNT = arma::ones<arma::vec>(iJ);
    for (int j = 0; j < iJ; ++j)
    {
      int jTM = iTM[j]-1;
      ExpAlphaNTnow = iExpAlphaNT[jTM];
      ExpAlphaTnow = iExpAlphaT[jTM];
      ExpAlphaORnow = iExpAlphaOR[jTM];
      ExpXBetaNTnow = iExpXBetaNT[j];
      ExpXBetaTnow = iExpXBetaT[j];
      ExpXBetaORnow = iExpXBetaOR[j];
      if (iRiskNT[j]==0) {
        iProbTafterNT = (ExpAlphaTnow*ExpXBetaTnow*exp(betay)) / (1 + (ExpAlphaTnow*ExpXBetaTnow*exp(betay)));
        if (j < iJ) { 
          iRiskNT[j+1] = 0;}
        if(iYT[j]==1) {
          iContrib = log(iProbTafterNT);
          }
        else {
          iContrib = log(1-iProbTafterNT);
          }} 
      else {
        iProb1 = ExpAlphaNTnow*ExpXBetaNTnow/(1 + (ExpAlphaNTnow*ExpXBetaNTnow));
        iProb2 = ExpAlphaTnow*ExpXBetaTnow/(1 + ExpAlphaTnow*ExpXBetaTnow);
        iOR = ExpAlphaORnow*ExpXBetaORnow;
        if ((iOR > 1 - epsOR) & (iOR < 1 + epsOR))
        {
          iProb12 = iProb1*iProb2;
        } 
        else {
          iProb12 = (1 + (iProb1 + iProb2)*(iOR - 1) - sqrt(pow(1 + (iProb1 + iProb2)*(iOR - 1), 2.0) -
            4*iOR*(iOR - 1)*iProb1*iProb2)) / (2 * (iOR - 1));
          }
        if (iYNT[j]==1 && iYT[j]==1) {
          iContrib = log(iProb12);
          }
        if (iYNT[j]==1 && iYT[j]==0) {
          iContrib = log(iProb1 - iProb12);
          if (j < iJ) { 
            iRiskNT[j+1] = 0;}
          }
        if (iYNT[j]==0 && iYT[j]==1) {
          iContrib = log(iProb2 - iProb12);
          }
        if (iYNT[j]==0 && iYT[j]==0) {
          iContrib = log(1 - iProb1 - iProb2 + iProb12);
          }}
      loglik += iContrib;
  }}
    double  penalloglik = loglik - as_scalar(penaltermNT) - 
      as_scalar(penaltermT) - as_scalar(penaltermOR);
    return(-penalloglik);
}
