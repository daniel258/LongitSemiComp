#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
//' @export
// [[Rcpp::export]]
double PenalLogLikTimeDepOLD(arma::vec param, arma::mat YT, arma::mat YNT, arma::mat riskNT, arma::mat riskT, 
                   arma::mat XNT, arma::mat XT, arma::mat XOR,
                   arma::mat XNTtimeDep, arma::mat XTtimeDep, arma::mat XORtimeDep,
                   arma::mat TimeBase, arma::mat TimePen, arma::vec lambda, double epsOR)
{
  //int pNTtimeDep, int pTtimeDep, int pORtimeDep
  int n = YT.n_rows;
  int J = YT.n_cols;
  int pNT = XNT.n_cols;
  int pT = XT.n_cols;
  int pOR = XOR.n_cols;
  int Q = TimeBase.n_cols; // Q is the the number of B-splines (number of rows in TimeBase should be J)
  // Decomposing "param" to individual paramters 
  // Current verison of the code: only one time-dep variable (but throgout the code, I made preperations for generalizing it)
  int pNTtimeDep = 1;
  int pTtimeDep = 1;
  int pORtimeDep = 1;
  double gamma = param[0];
  arma::vec alphaNT = param.subvec(1,Q);
  arma::vec alphaT = param.subvec(Q+1,2*Q);
  arma::vec alphaOR = param.subvec(2*Q+1,3*Q);
  arma::mat penaltermNT = lambda[0] * alphaNT.t() * TimePen * alphaNT;
  arma::mat penaltermT = lambda[1] * alphaT.t() * TimePen * alphaT;
  arma::mat penaltermOR = lambda[2] * alphaOR.t() * TimePen * alphaOR;
  // beta order is betaNT,betaNTtimeDep, betaT, betaTtimeDep, betaOR, betaORtimeDep
  arma::vec betaNT = param.subvec(3*Q + 1, 3*Q + pNT);
  arma::vec betaNTtimeDep = param.subvec(3*Q + pNT +  1, 3*Q + pNT + pNTtimeDep);
  arma::vec betaT = param.subvec(3*Q + pNT + pNTtimeDep + 1, 3*Q + pNT + pNTtimeDep + pT);
  arma::vec betaTtimeDep = param.subvec(3*Q + pNT + pNTtimeDep + pT + 1, 3*Q + pNT + pNTtimeDep + pT + pTtimeDep);
  arma::vec betaOR = param.subvec(3*Q + pNT + pNTtimeDep + pT + pTtimeDep + 1, 3*Q + pNT + pNTtimeDep + pT + pTtimeDep + pOR);
  arma::vec betaORtimeDep = param.subvec(3*Q + pNT + pNTtimeDep + pT + pTtimeDep + pOR + 1, 
                                         3*Q + pNT + pNTtimeDep + pT + pTtimeDep + pOR + pORtimeDep);
//  Rcpp::Rcout << "betaTtimeDep:"<< betaTtimeDep << std::endl;
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
  arma::vec ExpXBetaNTfixed = exp(XNT * betaNT);
  arma::vec ExpXBetaTfixed = exp(XT * betaT);
  arma::vec ExpXBetaORfixed = exp(XOR * betaOR);
  arma::vec ExpAlphaNT = exp(TimeBase * alphaNT);
  arma::vec ExpAlphaT = exp(TimeBase * alphaT);
  arma::vec ExpAlphaOR = exp(TimeBase * alphaOR);
  for (int j = 0; j < J; ++j)
  {
    ExpAlphaNTnow = ExpAlphaNT[j];
    ExpAlphaTnow = ExpAlphaT[j];
    ExpAlphaORnow = ExpAlphaOR[j];
    // Take current value of time-dep variables
    arma::vec ExpXBetaNTtimeDep = arma::zeros(n);
    arma::vec ExpXBetaTtimeDep = arma::zeros(n);
    arma::vec ExpXBetaORtimeDep = arma::zeros(n);
    arma::vec XNTtimeDepNow = XNTtimeDep.col(j);
    arma::vec XTtimeDepNow = XTtimeDep.col(j);
    arma::vec XORtimeDepNow = XORtimeDep.col(j);
    for(int k = 0; k < n; ++k)
    {
      arma::uvec indexNT =  as<arma::uvec>(wrap(Rcpp::Range((k-1)*pNTtimeDep + 1, k*pNTtimeDep)));
      arma::uvec indexT =  as<arma::uvec>(wrap(Rcpp::Range((k-1)*pTtimeDep + 1, k*pTtimeDep)));
      arma::uvec indexOR =  as<arma::uvec>(wrap(Rcpp::Range((k-1)*pORtimeDep + 1, k*pORtimeDep)));
      arma::rowvec tempNT = XNTtimeDepNow(indexNT);
      arma::rowvec tempT = XTtimeDepNow(indexT);
      arma::rowvec tempOR = XORtimeDepNow(indexOR);
     // Rcpp::Rcout << "tempT:"<< tempT << std::endl;
      ExpXBetaNTtimeDep[k] = exp(as_scalar(tempNT * betaNTtimeDep));
      ExpXBetaTtimeDep[k] = exp(as_scalar(tempT * betaTtimeDep));
      ExpXBetaORtimeDep[k] = exp(as_scalar(tempOR * betaORtimeDep));
    }
    for (int i = 0; i < n; ++i)
    {
      if (riskT(i,j)==0) {
        iContrib=0;
        //   nocontrib
      } else {
        if (riskNT(i,j)==0) {
          iProbTafterNT = (ExpAlphaTnow*ExpXBetaTfixed[i]*ExpXBetaTtimeDep[i]*exp(gamma))/
            (1 + (ExpAlphaTnow*ExpXBetaTfixed[i]*ExpXBetaTtimeDep[i]*exp(gamma)));
        //  Rcpp::Rcout << "i:"<< i + 1 << std::endl;
        //  Rcpp::Rcout << "ExpXBetaTtimeDep:"<< ExpXBetaTtimeDep[i] << std::endl;
            if(YT(i,j)==1) {
            iContrib = log(iProbTafterNT);
          }
          else {
            iContrib = log(1-iProbTafterNT);
          }
      //    Rcpp::Rcout << "iProbTafterNT  "<< iProbTafterNT << std::endl;
        } else {
          iProb1 = ExpAlphaNTnow*ExpXBetaNTfixed[i]*ExpXBetaNTtimeDep[i]/(1 + (ExpAlphaNTnow*ExpXBetaNTfixed[i]*ExpXBetaNTtimeDep[i]));
          iProb2 = ExpAlphaTnow*ExpXBetaTfixed[i]*ExpXBetaTtimeDep[i]/(1 + ExpAlphaTnow*ExpXBetaTfixed[i]*ExpXBetaTtimeDep[i]);
          iOR = ExpAlphaORnow*ExpXBetaORfixed[i]*ExpXBetaORtimeDep[i];
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
