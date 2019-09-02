#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
//' @export
// [[Rcpp::export]]
arma::mat GradPenalLogLikPersTimeDep(arma::vec param, arma::vec ID, arma::uvec TM,
                      arma::vec YT, arma::vec YNT, 
                   arma::mat XNT, arma::mat XT, arma::mat XOR,
                   arma::mat TimeBase, arma::mat TimePen, arma::vec lambda, double epsOR)
{
  arma::vec IDunq = arma::unique(ID); // IDs once
  int n = IDunq.size(); // sample size
  arma::uvec TMunq = arma::unique(TM); // unique "times" (intervals)
 // int J = TMunq.size(); // number of "times" (intervals)
 // int nrose = XNT.n_rows;
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
  arma::vec penaltermNT = 2 * lambda[0] * TimePen * alphaNT;
  arma::vec penaltermT = 2 * lambda[1] * TimePen * alphaT;
  arma::vec penaltermOR = 2 * lambda[2]  * TimePen * alphaOR;
  // Rcpp::Rcout << "penaltermNT:  "<< penaltermNT << std::endl;
  // Rcpp::Rcout << "penaltermT:  "<< penaltermT << std::endl;
  // Rcpp::Rcout << "penaltermOR:  "<< penaltermOR << std::endl;
  arma::vec betaNT = param.subvec(3*Q + 1, 3*Q + pNT);
  arma::vec betaT = param.subvec(3*Q + pNT + 1, 3*Q + pNT + pT);
  arma::vec betaOR = param.subvec(3*Q + pNT + pT + 1, 3*Q + pNT + pT + pOR);
  arma::rowvec Grad(1 + 3*Q + pNT + pT + pOR);
  arma::rowvec iGrad(1 + 3*Q + pNT + pT + pOR);
  arma::mat GradSquare(1 + 3*Q + pNT + pT + pOR, 1 + 3*Q + pNT + pT + pOR);
  // Rcpp::Rcout << "betaNT:  "<< betaNT << std::endl;
  // Rcpp::Rcout << "betaT:  "<< betaT << std::endl;
  // Rcpp::Rcout << "betaOR:  "<< betaOR << std::endl;
 // double loglik = 0;
 // double iContrib = 0;
  double iProbTafterNT = 0;
  double iProb1 = 0;
  double iProb2 = 0;
  double iOR = 0;
  double iProb12 = 0;
  Grad.fill(0);
  iGrad.fill(0);
  GradSquare.fill(0);
  double cij = 0;
  double nuij = 0;
  double ExpAlphaNTnow = 0;
  double ExpAlphaTnow = 0;
  double ExpAlphaORnow = 0;
  double ExpXBetaNTnow = 0;
  double ExpXBetaTnow = 0;
  double ExpXBetaORnow = 0;
  //double b = 0;
  arma::rowvec iXNTnow(pT);
  arma::rowvec iXTnow(pT);
  arma::rowvec iXORnow(pOR);
  arma::rowvec bnow;
  iXNTnow.fill(0);  
  iXTnow.fill(0);  
  iXORnow.fill(0);
  bnow.fill(0);
  arma::vec ExpAlphaNT = exp(TimeBase * alphaNT);
  arma::vec ExpAlphaT = exp(TimeBase * alphaT);
  arma::vec ExpAlphaOR = exp(TimeBase * alphaOR);
  arma::vec ExpXBetaNT = exp(XNT * betaNT);
  arma::vec ExpXBetaT = exp(XT * betaT);
  arma::vec ExpXBetaOR = exp(XOR * betaOR);
 // Rcpp::Rcout << "ExpAlphaNT:  "<< ExpAlphaNT << std::endl;
  for (int i = 0; i < n; ++i)
  {
   // Rcpp::Rcout << "i =   "<< i << std::endl;
    int iID = IDunq[i]; // get observation ID
    // Get all data on that observation
    arma::uvec pos = find(ID==iID); // rows associated with the observation
    arma::vec iYNT = YNT.elem(pos); // YNT for that observation
    arma::vec iYT = YT.elem(pos); // YT for that observation
    arma::uvec iTM = TM.elem(pos); // TM for that observation, they are already unique
    arma::mat iXNT = XNT.rows(pos);
    arma::mat iXT = XT.rows(pos);
    arma::mat iXOR = XOR.rows(pos);
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
    // Rcpp::Rcout << "j =   "<< j << std::endl;
    int jTM = iTM[j]-1;
    ExpAlphaNTnow = iExpAlphaNT[jTM];
    ExpAlphaTnow = iExpAlphaT[jTM];
    ExpAlphaORnow = iExpAlphaOR[jTM];
    ExpXBetaNTnow = iExpXBetaNT[j];
    ExpXBetaTnow = iExpXBetaT[j];
    ExpXBetaORnow = iExpXBetaOR[j];
    iXNTnow = iXNT.row(j);
    iXTnow = iXT.row(j);
    iXORnow = iXOR.row(j);
    bnow = TimeBase.row(j);
    iGrad.fill(0);
    if (iRiskNT[j]==0) {
          iProbTafterNT = (ExpAlphaTnow*ExpXBetaTnow*exp(betay)) / (1 + (ExpAlphaTnow*ExpXBetaTnow*exp(betay)));
          if(iYT[j]==1) {
            iGrad[0] = 1-iProbTafterNT;
            iGrad(arma::span(1 + Q, 2*Q))= bnow*(1-iProbTafterNT);
            iGrad(arma::span(3*Q + pNT + 1, 3*Q + pNT + pT)) = iXTnow*(1-iProbTafterNT);
            } 
          else {
            iGrad[0] = -iProbTafterNT;
            iGrad(arma::span(1 + Q, 2*Q)) = -bnow*iProbTafterNT;
            iGrad(arma::span(3*Q + pNT + 1, 3*Q + pNT + pT)) = -iXTnow*iProbTafterNT;    
              }} else {
          iGrad[0] = 0;
          iProb1 = ExpAlphaNTnow*ExpXBetaNTnow/(1 + (ExpAlphaNTnow*ExpXBetaNTnow));
          iProb2 = ExpAlphaTnow*ExpXBetaTnow/(1 + ExpAlphaTnow*ExpXBetaTnow);
          iOR = ExpAlphaORnow*ExpXBetaORnow;
          if ((iOR > 1 - epsOR) & (iOR < 1 + epsOR))
          {
            iProb12 = iProb1*iProb2;
          // Rcpp::Rcout << "iProb12  "<< iProb12 << std::endl;
          if (iYNT[j]==1 && iYT[j]==1) {
            iGrad(arma::span(1, Q)) = bnow*(1 - iProb1);
            iGrad(arma::span(Q + 1,2*Q)) = bnow*(1 - iProb2);
            iGrad(arma::span(2*Q + 1,3*Q)).fill(0);
            iGrad(arma::span(3*Q + 1, 3*Q + pNT)) = iXNTnow*(1 - iProb1);
            iGrad(arma::span(3*Q + pNT + 1, 3*Q + pNT + pT)) = iXTnow*(1 - iProb2);
            iGrad(arma::span(3*Q + pNT + pT + 1, 3*Q + pNT + pT + pOR)).fill(0);
            }
          if (iYNT[j]==1 && iYT[j]==0) {
            iGrad(arma::span(1, Q)) =  bnow*(1 - iProb1);
            iGrad(arma::span(Q + 1, 2*Q)) = -bnow*iProb2;
            iGrad(arma::span(2*Q + 1, 3*Q)).fill(0);
            iGrad(arma::span(3*Q + 1, 3*Q + pNT)) = iXNTnow*(1 - iProb1);
            iGrad(arma::span(3*Q + pNT + 1, 3*Q + pNT + pT)) = -iXTnow*iProb2;
            iGrad(arma::span(3*Q + pNT + pT + 1, 3*Q + pNT + pT + pOR)).fill(0);
            if (j < iJ) { 
             iRiskNT[j+1] = 0;}
          }
          if (iYNT[j]==0 && iYT[j]==1) {
            iGrad(arma::span(1, Q)) =  -bnow*iProb1;
            iGrad(arma::span(Q + 1, 2*Q)) = bnow*(1 - iProb2);
            iGrad(arma::span(2*Q + 1, 3*Q)).fill(0);
            iGrad(arma::span(3*Q + 1, 3*Q + pNT)) = -iXNTnow*iProb1;
            iGrad(arma::span(3*Q + pNT + 1, 3*Q + pNT + pT)) = iXTnow*(1 - iProb2);
            iGrad(arma::span(3*Q + pNT + pT + 1, 3*Q + pNT + pT + pOR)).fill(0);
          }
          if (iYNT[j]==0 && iYT[j]==0) {
            iGrad(arma::span(1, Q)) =  -bnow * iProb1;
            iGrad(arma::span(Q + 1, 2*Q)) = -bnow * iProb2;
            iGrad(arma::span(2*Q + 1, 3*Q)).fill(0);
            iGrad(arma::span(3*Q + 1, 3*Q + pNT)) = -iXNTnow*iProb1;
            iGrad(arma::span(3*Q + pNT + 1, 3*Q + pNT + pT)) = -iXTnow*iProb2;
            iGrad(arma::span(3*Q + pNT + pT + 1, 3*Q + pNT + pT + pOR)).fill(0);
            }
        } 
      else {
        cij = (iProb1 + iProb2)*(iOR - 1);
        nuij = sqrt(pow(1 + cij, 2.0) - 4*iOR*(iOR - 1)*iProb1*iProb2);
        iProb12 = (1 + cij - nuij) / (2 * (iOR - 1));
        // Rcpp::Rcout << "iProb12  "<< iProb12 << std::endl;
        if (iYNT[j]==1 && iYT[j]==1) {
          iGrad(arma::span(1, Q)) = 0.5*bnow*iProb1*(1 - iProb1) * (nuij - 1 - cij + 2*iOR*iProb2)/(nuij*iProb12);
          iGrad(arma::span(Q + 1, 2*Q)) = 0.5*bnow*iProb2*(1 - iProb2) * (nuij - 1 - cij + 2*iOR*iProb1)/(nuij*iProb12);
          iGrad(arma::span(2*Q + 1, 3*Q)) = bnow*iOR*
           (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
                (2*(iOR - 1)*(iOR - 1)*iProb12);
          iGrad(arma::span(3*Q + 1, 3*Q + pNT)) = 0.5*iXNTnow*iProb1*(1 - iProb1) * (nuij - 1 - cij + 2*iOR*iProb2)/(nuij*iProb12);
          iGrad(arma::span(3*Q + pNT + 1, 3*Q + pNT + pT)) = 0.5*iXTnow*iProb2*(1 - iProb2) * (nuij - 1 - cij + 2*iOR*iProb1)/(nuij*iProb12);
          iGrad(arma::span(3*Q + pNT + pT + 1, 3*Q + pNT + pT + pOR)) = iXORnow*iOR*
                    (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
                      (2*(iOR - 1)*(iOR - 1)*iProb12);
                }
        if (iYNT[j]==1 && iYT[j]==0) {
          iGrad(arma::span(1, Q)) =  ((bnow*iProb1*(1 - iProb1)/(iProb1 - iProb12))) *(1 -  0.5*(nuij - 1 - cij + 2*iOR*iProb2)/nuij);
          iGrad(arma::span(Q + 1, 2*Q)) = -0.5*bnow*iProb2*(1 - iProb2) * (nuij - 1 - cij + 2*iOR*iProb1)/(nuij*(iProb1 - iProb12));
          iGrad(arma::span(2*Q + 1, 3*Q)) = -bnow*iOR*
            (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
                (2*(iOR - 1)*(iOR - 1)*(iProb1 - iProb12));
          iGrad(arma::span(3*Q + 1, 3*Q + pNT)) = ((iXNTnow*iProb1*(1 - iProb1)/(iProb1 - iProb12))) *(1 -  0.5*(nuij - 1 - cij + 2*iOR*iProb2)/nuij);
          iGrad(arma::span(3*Q + pNT + 1, 3*Q + pNT + pT)) = -0.5*iXTnow*iProb2*(1 - iProb2) * (nuij - 1 - cij + 2*iOR*iProb1)/(nuij*(iProb1 - iProb12));
          iGrad(arma::span(3*Q + pNT + pT + 1, 3*Q + pNT + pT + pOR)) =  -iXORnow*iOR*
                    (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
                      (2*(iOR - 1)*(iOR - 1)*(iProb1 - iProb12));
          if (j < iJ) { 
            iRiskNT[j+1] = 0;}
              }
        if (iYNT[j]==0 && iYT[j]==1) {
          iGrad(arma::span(1, Q)) =  -0.5*bnow*iProb1*(1 - iProb1) * (nuij - 1 - cij + 2*iOR*iProb2)/(nuij*(iProb2 - iProb12));
          iGrad(arma::span(Q + 1, 2*Q)) = ((bnow*iProb2*(1 - iProb2)/(iProb2 - iProb12))) *(1 -  0.5*(nuij - 1 - cij + 2*iOR*iProb1)/nuij);
          iGrad(arma::span(2*Q + 1, 3*Q)) = -bnow*iOR*
            (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
                      (2*(iOR - 1)*(iOR - 1)*(iProb2 - iProb12));
          iGrad(arma::span(3*Q + 1, 3*Q + pNT)) = -0.5*iXNTnow*iProb1*(1 - iProb1) * (nuij - 1 - cij + 2*iOR*iProb2)/(nuij*(iProb2 - iProb12));
          iGrad(arma::span(3*Q + pNT + 1, 3*Q + pNT + pT)) = ((iXTnow*iProb2*(1 - iProb2)/(iProb2 - iProb12))) *(1 -  0.5*(nuij - 1 - cij + 2*iOR*iProb1)/nuij);
          iGrad(arma::span(3*Q + pNT + pT + 1, 3*Q + pNT + pT + pOR)) =  -iXORnow*iOR*
            (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
                    (2*(iOR - 1)*(iOR - 1)*(iProb2 - iProb12));
                }
        if (iYNT[j]==0 && iYT[j]==0) {
          iGrad(arma::span(1, Q)) =  ((bnow*iProb1*(1 - iProb1)/(1 - iProb1 - iProb2 + iProb12))) *(0.5*(nuij - 1 - cij + 2*iOR*iProb2)/nuij - 1);
          iGrad(arma::span(Q + 1, 2*Q)) = ((bnow*iProb2*(1 - iProb2)/(1 - iProb1 - iProb2 + iProb12))) *(0.5*(nuij - 1 - cij + 2*iOR*iProb1)/nuij - 1);
          iGrad(arma::span(2*Q + 1, 3*Q)) =  bnow*iOR*
            (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
                      (2*(iOR - 1)*(iOR - 1)*(1 - iProb1 - iProb2 + iProb12));
          iGrad(arma::span(3*Q + 1, 3*Q + pNT)) = ((iXNTnow*iProb1*(1 - iProb1)/(1 - iProb1 - iProb2 + iProb12))) *(0.5*(nuij - 1 - cij + 2*iOR*iProb2)/nuij - 1);
          iGrad(arma::span(3*Q + pNT + 1, 3*Q + pNT + pT)) = ((iXTnow*iProb2*(1 - iProb2)/(1 - iProb1 - iProb2 + iProb12))) *(0.5*(nuij - 1 - cij + 2*iOR*iProb1)/nuij - 1);
          iGrad(arma::span(3*Q + pNT + pT + 1, 3*Q + pNT + pT + pOR)) =  iXORnow*iOR*
                    (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
                      (2*(iOR - 1)*(iOR - 1)*(1 - iProb1 - iProb2 + iProb12));
        }}}
  }}
  iGrad(arma::span(1,Q)) -= penaltermNT.t()/n;
  iGrad(arma::span(Q + 1, 2*Q)) -= penaltermT.t()/n;
  iGrad(arma::span(2*Q + 1, 3*Q)) -= penaltermOR.t()/n;
  GradSquare += iGrad.t()*iGrad;
  return(GradSquare);
}