#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
//' @export
// [[Rcpp::export]]
arma::vec GradPenalLogLikTimeDepOLD(arma::vec param, arma::mat YNT, arma::mat riskNT, arma::mat riskT, arma::mat YT, 
                          arma::mat XNT, arma::mat XT, arma::mat XOR, 
                          arma::mat XNTtimeDep, arma::mat XTtimeDep, arma::mat XORtimeDep,
                          arma::mat TimeBase, arma::mat TimePen, arma::vec lambda, double epsOR)
     {
  //                        int pNTtimeDep, int pTtimeDep, int pORtimeDep,
  int n = YT.n_rows;
  int J = YT.n_cols;
  int pNT = XNT.n_cols;
  int pT = XT.n_cols;
  int pOR = XOR.n_cols;
  // Q is the the number of B-splines (number of rows in TimeBase should be J)
  // Decomposing "param" to individual paramters as in PenalLogLik
  int Q = TimeBase.n_cols;
  // Current verison of the code: only one time-dep variable (but throgout the code, I made preperations for generalizing it)
  int pNTtimeDep = 1;
  int pTtimeDep = 1;
  int pORtimeDep = 1;
  double gamma = param[0];
  arma::vec alphaNT = param.subvec(1, Q);
  arma::vec alphaT = param.subvec(Q + 1, 2*Q);
  arma::vec alphaOR = param.subvec(2*Q + 1, 3*Q);
  arma::vec penaltermNT = 2 * lambda[0] * TimePen * alphaNT;
  arma::vec penaltermT = 2 * lambda[1] * TimePen * alphaT;
  arma::vec penaltermOR = 2 * lambda[2]  * TimePen * alphaOR;
  // beta order is betaNT,betaNTtimeDep, betaT, betaTtimeDep, betaOR, betaORtimeDep
  arma::vec betaNT = param.subvec(3*Q + 1, 3*Q + pNT);
  arma::vec betaNTtimeDep = param.subvec(3*Q + pNT +  1, 3*Q + pNT + pNTtimeDep);
  arma::vec betaT = param.subvec(3*Q + pNT + pNTtimeDep + 1, 3*Q + pNT + pNTtimeDep + pT);
  arma::vec betaTtimeDep = param.subvec(3*Q + pNT + pNTtimeDep + pT + 1, 3*Q + pNT + pNTtimeDep + pT + pTtimeDep);
  arma::vec betaOR = param.subvec(3*Q + pNT + pNTtimeDep + pT + pTtimeDep + 1, 3*Q + pNT + pNTtimeDep + pT + pTtimeDep + pOR);
  arma::vec betaORtimeDep = param.subvec(3*Q + pNT + pNTtimeDep + pT + pTtimeDep + pOR + 1, 
                                         3*Q + pNT + pNTtimeDep + pT + pTtimeDep + pOR + pORtimeDep);
  arma::vec Grad(1 + 3*Q + pNT + pNTtimeDep + pT + pTtimeDep + pOR + pORtimeDep);
  arma::vec iGrad(1 + 3*Q + pNT + pNTtimeDep + pT + pTtimeDep + pOR + pORtimeDep);
  Grad.fill(0);
  iGrad.fill(0);
  double iProbTafterNT = 0;
  double iProb1 = 0;
  double iProb2 = 0;
  double iOR = 0;
  double iProb12 = 0;
  double cij = 0;
  double nuij = 0;
 // double x = 0;
  double b = 0;
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
    for(int m = 0; m < n; ++m)
    {
      arma::uvec indexNT =  as<arma::uvec>(wrap(Rcpp::Range((m - 1)* pNTtimeDep + 1, m*pNTtimeDep)));
      arma::uvec indexT =  as<arma::uvec>(wrap(Rcpp::Range((m - 1) * pTtimeDep + 1, m*pTtimeDep)));
      arma::uvec indexOR =  as<arma::uvec>(wrap(Rcpp::Range((m - 1)* pORtimeDep + 1, m*pORtimeDep)));
      arma::rowvec tempNT = XNTtimeDepNow(indexNT);
      arma::rowvec tempT = XTtimeDepNow(indexT);
      arma::rowvec tempOR = XORtimeDepNow(indexOR);
      // Rcpp::Rcout << "tempT:"<< tempT << std::endl;
      ExpXBetaNTtimeDep[m] = exp(as_scalar(tempNT * betaNTtimeDep));
      ExpXBetaTtimeDep[m] = exp(as_scalar(tempT * betaTtimeDep));
      ExpXBetaORtimeDep[m] = exp(as_scalar(tempOR * betaORtimeDep));
    }
    for (int i = 0; i < n; ++i)
    {
      iGrad.fill(0);
        if (riskT(i,j)==0) {
          iGrad.fill(0);
          //   no contribution to the gradient
        } else {
          if (riskNT(i,j)==0) {
            iProbTafterNT = (ExpAlphaTnow*ExpXBetaTfixed[i]*ExpXBetaTtimeDep[i]*exp(gamma))/
              (1 + (ExpAlphaTnow*ExpXBetaTfixed[i]*ExpXBetaTtimeDep[i]*exp(gamma)));
            if(YT(i,j)==1) {
              iGrad[0] = 1-iProbTafterNT;
              for (int q=0; q < Q; ++q)
              {
                b = TimeBase(j,q);
                iGrad[1 + Q + q] = b*(1 - iProbTafterNT);
              }
              for (int k =0; k < pT; ++k)
              {
                iGrad[3*Q + pNT + pNTtimeDep + k + 1] = XT(i,k)*(1 - iProbTafterNT);
              }
            //  for (int k =0; k < pTtimeDep; ++k)
            //  {
            
                iGrad[3*Q + pNT + pNTtimeDep + pT + 0 + 1] = XTtimeDepNow[i]*(1 - iProbTafterNT); // only one time-dep var
              // Rcpp::Rcout << "i:"<< i + 1 << std::endl;
              // Rcpp::Rcout << "XTtimeDepNow:"<< XTtimeDepNow[i+1] << std::endl;
              //  Rcpp::Rcout << "iGrad[3*Q + pNT + pNTtimeDep + pT + 0 + 1]:"<< iGrad[3*Q + pNT + pNTtimeDep + pT + 0 + 1] << std::endl;
              // // only one time-dep var
              //}
              }
            else {
              iGrad[0] = -iProbTafterNT;
              for (int q=0; q < Q; ++q)
              {
                b = TimeBase(j,q);
                iGrad[1 + Q + q] = -b*iProbTafterNT;
              }
              for (int k =0; k < pT; ++k)
              {
                iGrad[3*Q + pNT + pNTtimeDep + k + 1] = -XT(i,k)*iProbTafterNT;
              }
              //for (int k =0; k < pTtimeDep; ++k)
              //{
                iGrad[3*Q + pNT + pNTtimeDep + pT + 0 + 1] = -XTtimeDepNow[i]*iProbTafterNT;  // only one time-dep var
          //    }
            }
          } else {
            iProb1 = ExpAlphaNTnow*ExpXBetaNTfixed[i]*ExpXBetaNTtimeDep[i]/
              (1 + (ExpAlphaNTnow*ExpXBetaNTfixed[i]*ExpXBetaNTtimeDep[i]));
            iProb2 = ExpAlphaTnow*ExpXBetaTfixed[i]*ExpXBetaTtimeDep[i]/(1 + ExpAlphaTnow*ExpXBetaTfixed[i]*ExpXBetaTtimeDep[i]);
            iOR = ExpAlphaORnow*ExpXBetaORfixed[i]*ExpXBetaORtimeDep[i];
            if ((iOR > 1 - epsOR) & (iOR < 1 + epsOR))
            {
             iProb12 = iProb1*iProb2;
             if (YNT(i,j)==1 && YT(i,j)==1) {
                for (int q=0; q < Q; ++q)
                {
                  b = TimeBase(j,q);
                  iGrad[1 + q] = b*(1 - iProb1);
                  iGrad[Q + 1 + q] = b*(1 - iProb2);
                  iGrad[2*Q + 1 + q] = 0;
                }
                for (int k =0; k < pNT; ++k)
                {
                  iGrad[3*Q + k + 1] = XNT(i,k)*(1 - iProb1);
                }
                //for (int k =0; k < pNTtimeDep; ++k)
                //{
                  iGrad[3*Q + pNT + 0 + 1] = XNTtimeDepNow[i]*(1 - iProb1); // only one time-dep var
                //}
                for (int k =0; k < pT; ++k)
                {
                  iGrad[3*Q + pNT + pNTtimeDep + k + 1] = XT(i,k)*(1 - iProb2);
                }
                //for (int k =0; k < pTtimeDep; ++k)
                //{
                  iGrad[3*Q + pNT + pNTtimeDep + pT + 0 + 1] = XTtimeDepNow[i]*(1 - iProb2); // only one time-dep var
                //}
                for (int k =0; k < pOR; ++k)
                {
                  iGrad[3*Q + pNT + pNTtimeDep + pT + pTtimeDep + k + 1] = 0;
                }
                //for (int k =0; k < pORtimeDep; ++k)
                //{
                  iGrad[3*Q + pNT + pNTtimeDep + pT + pTtimeDep + pOR + 0 + 1] = 0;// only one time-dep var
                //}
                }
              if (YNT(i,j)==1 && YT(i,j)==0) {
                for (int q=0; q < Q; ++q)
                {
                  b = TimeBase(j,q);
                  iGrad[1 + q] =  b*(1 - iProb1);
                  iGrad[Q + 1 + q] = -b*iProb2;
                  iGrad[2*Q + 1 + q] = 0;
                }
                for (int k =0; k < pNT; ++k)
                {
                  iGrad[3*Q + k + 1] = XNT(i,k)*(1 - iProb1);
                }
                //for (int k =0; k < pNTtimeDep; ++k)
                //{
                  iGrad[3*Q + pNT + 0 + 1] = XNTtimeDepNow[i]*(1 - iProb1); // only one time-dep var
                //}
                for (int k =0; k < pT; ++k)
                {
                  iGrad[3*Q + pNT + pNTtimeDep + k + 1] = -XT(i,k)*iProb2;
                }
                //for (int k =0; k < pTtimeDep; ++k)
                //{
                  iGrad[3*Q + pNT + pNTtimeDep + pT + 0 + 1] = -XTtimeDepNow[i]*iProb2; // only one time-dep var
                //}
                for (int k = 0; k < pOR; ++k)
                {
                  iGrad[3*Q + pNT + pNTtimeDep + pT + pTtimeDep + k + 1] =  0;
                }
                //for (int k = 0; k < pORtimeDep; ++k)
                //{
                  iGrad[3*Q + pNT + pNTtimeDep + pT + pTtimeDep + pOR + 0 + 1] =  0;// only one time-dep var
                //}
                }
               if (YNT(i,j)==0 && YT(i,j)==1) {
                 for (int q = 0; q < Q; ++q)
                 {
                   b = TimeBase(j,q);
                   iGrad[1 + q] =  -b*iProb1;
                   iGrad[Q + 1 + q] = b*(1 - iProb2);
                   iGrad[2*Q + 1 + q] = 0;
                 }
                 for (int k =0; k < pNT; ++k)
                 {
                   iGrad[3*Q + k + 1] = -XNT(i,k)*iProb1;
                 }
                 //for (int k =0; k < pNTtimeDep; ++k)
                 //{
                   iGrad[3*Q + pNT + 0 + 1] = -XNTtimeDepNow[i]*iProb1; // only one time-dep var
                 //}
                 for (int k =0; k < pT; ++k)
                 {
                   iGrad[3*Q + pNT + pNTtimeDep + k + 1] = XT(i,k)*(1 - iProb2);
                 }
                 //for (int k =0; k < pTtimeDep; ++k)
                 //{
                   iGrad[3*Q + pNT + pNTtimeDep + pT + 0 + 1] = XTtimeDepNow[i]*(1 - iProb2); // only one time-dep var
                 //}
                 for (int k =0; k < pOR; ++k)
                 {
                   iGrad[3*Q + pNT + pNTtimeDep + pT + pTtimeDep + k + 1] =  0;
                 }
                 //for (int k =0; k < pORtimeDep; ++k)
                 //{
                   iGrad[3*Q + pNT + pNTtimeDep + pT + pTtimeDep + pOR + 0 + 1] =  0; // only one time-dep var
                 //}
                 }
               if (YNT(i,j)==0 && YT(i,j)==0) {
                 for (int q=0; q < Q; ++q)
                 {
                   b = TimeBase(j,q);
                   iGrad[1 + q] =  -b*iProb1;
                   iGrad[Q + 1 + q] = -b*iProb2;
                   iGrad[2*Q + 1 + q] =  0;
                 }
                 for (int k =0; k < pNT; ++k)
                 {
                 iGrad[3*Q + k + 1] = -XNT(i,k)*iProb1;
                 }
                 //for (int k =0; k < pNTtimeDep; ++k)
                 //{
                   iGrad[3*Q + pNT + 0 + 1] = -XNTtimeDepNow[i]*iProb1; // only one time-dep var
                 //}
                  for (int k =0; k < pT; ++k)
                  {
                   iGrad[3*Q + pNT + pNTtimeDep + k + 1] = -XT(i,k)*iProb2;
                  }
                  //for (int k =0; k < pTtimeDep; ++k)
                  //{
                    iGrad[3*Q + pNT + pNTtimeDep + pT + 0 + 1] = -XTtimeDepNow[i]*iProb2; // only one time-dep var
                  //}
                  for (int k =0; k < pOR; ++k)
                  {
                   iGrad[3*Q + pNT + pNTtimeDep + pT + pTtimeDep + k + 1] =  0;
                  }
                  //for (int k =0; k < pORtimeDep; ++k)
                  //{
                    iGrad[3*Q + pNT + pNTtimeDep + pT + pTtimeDep + pOR + 0 + 1] =  0; // only one time-dep var
                  //}
                  }
            } else {
              cij = (iProb1 + iProb2)*(iOR - 1);
              nuij = sqrt(pow(1 + cij, 2.0) - 4*iOR*(iOR - 1)*iProb1*iProb2);
              iProb12 = (1 + cij - nuij) / (2 * (iOR - 1));

            // Rcpp::Rcout << "iProb12  "<< iProb12 << std::endl;
            if (YNT(i,j)==1 && YT(i,j)==1) {
              for (int q=0; q < Q; ++q)
              {
                b = TimeBase(j,q);
                iGrad[1 + q] = 0.5*b*iProb1*(1 - iProb1) * (nuij - 1 - cij + 2*iOR*iProb2)/(nuij*iProb12);
                iGrad[Q + 1 + q] = 0.5*b*iProb2*(1 - iProb2) * (nuij - 1 - cij + 2*iOR*iProb1)/(nuij*iProb12);
                iGrad[2*Q + 1 + q] = b*iOR*
                  (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
                    (2*(iOR - 1)*(iOR - 1)*iProb12);
              }
              for (int k =0; k < pNT; ++k)
              {
                iGrad[3*Q + k + 1] = 0.5 * XNT(i,k) * iProb1 * (1 - iProb1) * 
                  (nuij - 1 - cij + 2 * iOR * iProb2) / (nuij * iProb12);
              }
              //for (int k =0; k < pNTtimeDep; ++k)
              //{
              // only one time-dep var
                iGrad[3*Q + pNT + 0 + 1] = 0.5 * XNTtimeDepNow[i] * iProb1 * (1 - iProb1) * 
                  (nuij - 1 - cij + 2 * iOR * iProb2) / (nuij * iProb12);
              //}
              for (int k =0; k < pT; ++k)
              {
                iGrad[3*Q + pNT + pNTtimeDep + k + 1] = 0.5 * XT(i,k) * iProb2 * (1 - iProb2) * 
                  (nuij - 1 - cij + 2 * iOR * iProb1) / (nuij * iProb12);
              }
              //for (int k =0; k < pTtimeDep; ++k)
              //{
              // only one time-dep var
                iGrad[3*Q + pNT + pNTtimeDep + pT + 0 + 1] = 0.5*XTtimeDepNow[i]*iProb2*(1 - iProb2) * 
                  (nuij - 1 - cij + 2*iOR*iProb1)/(nuij*iProb12);
              //}
              for (int k =0; k < pOR; ++k)
              {
                iGrad[3*Q + pNT + pNTtimeDep + pT + pTtimeDep + k + 1] = XOR(i,k)*iOR*
                  (((iOR-1)*((iProb1 + iProb2) - ((1 + cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
                  (2*(iOR - 1)*(iOR - 1)*iProb12);
              }
              //for (int k =0; k < pORtimeDep; ++k)
              //{
                iGrad[3*Q + pNT + pNTtimeDep + pT + pTtimeDep + pOR +  0 + 1] = XORtimeDepNow[i]*iOR*
                  (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
                    (2*(iOR - 1)*(iOR - 1)*iProb12);
              // only one time-dep var
              //}
              }
            if (YNT(i,j)==1 && YT(i,j)==0) {
              for (int q = 0; q < Q; ++q)
              {
                b = TimeBase(j,q);
                iGrad[1 + q] =  ((b*iProb1*(1 - iProb1)/(iProb1 - iProb12))) *(1 -  0.5*(nuij - 1 - cij + 2*iOR*iProb2)/nuij);
                iGrad[Q + 1 + q] = -0.5*b*iProb2*(1 - iProb2) * (nuij - 1 - cij + 2*iOR*iProb1)/(nuij*(iProb1 - iProb12));
                iGrad[2*Q + 1 + q] = -b*iOR*
                  (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
                    (2*(iOR - 1)*(iOR - 1)*(iProb1 - iProb12));
              }
      for (int k = 0; k < pNT; ++k)
      {
      iGrad[3*Q + k + 1] = ((XNT(i,k)*iProb1*(1 - iProb1)/(iProb1 - iProb12))) * 
        (1 -  0.5*(nuij - 1 - cij + 2*iOR*iProb2)/nuij);
      }
      //for (int k = 0; k < pNTtimeDep; ++k)
      //{
        iGrad[3*Q + pNT + 0 + 1] = ((XNTtimeDepNow[i]*iProb1*(1 - iProb1)/(iProb1 - iProb12))) * 
          (1 -  0.5*(nuij - 1 - cij + 2*iOR*iProb2)/nuij);
      // only one time-dep var
      //}
      for (int k =0; k < pT; ++k)
      {
      iGrad[3*Q + pNT + pNTtimeDep + k + 1] = -0.5*XT(i,k)*iProb2*(1 - iProb2) * 
        (nuij - 1 - cij + 2*iOR*iProb1)/(nuij*(iProb1 - iProb12));
      }
      //for (int k =0; k < pT; ++k)
      //{
        iGrad[3*Q + pNT + pNTtimeDep + pT + 0 + 1] = -0.5*XTtimeDepNow[i]*iProb2*(1 - iProb2) * 
          (nuij - 1 - cij + 2*iOR*iProb1)/(nuij*(iProb1 - iProb12));
      // only one time-dep var
      //}
      for (int k =0; k < pOR; ++k)
      {
      iGrad[3*Q + pNT + pNTtimeDep + pT + pTtimeDep + k + 1] =  -XOR(i,k)*iOR*
        (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
         (2*(iOR - 1)*(iOR - 1)*(iProb1 - iProb12));
      }
      //for (int k =0; k < pORtimeDep; ++k)
      //{
        iGrad[3*Q + pNT + pNTtimeDep + pT + pTtimeDep + pOR + 0 + 1] =  -XORtimeDepNow[i]*iOR*
          (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
            (2*(iOR - 1)*(iOR - 1)*(iProb1 - iProb12));
      //}
      // only one time-dep var
      }
      if (YNT(i,j)==0 && YT(i,j)==1) {
        for (int q=0; q < Q; ++q)
        {
          b = TimeBase(j,q);
          iGrad[1 + q] =  -0.5*b*iProb1*(1 - iProb1) * (nuij - 1 - cij + 2*iOR*iProb2)/(nuij*(iProb2 - iProb12));
          iGrad[Q + 1 + q] = ((b*iProb2*(1 - iProb2)/(iProb2 - iProb12))) *(1 -  0.5*(nuij - 1 - cij + 2*iOR*iProb1)/nuij);
          iGrad[2*Q + 1 + q] = -b*iOR*
            (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
              (2*(iOR - 1)*(iOR - 1)*(iProb2 - iProb12));
        }
      for (int k =0; k < pNT; ++k)
      {
      iGrad[3*Q + k + 1] = -0.5*XNT(i,k)*iProb1*(1 - iProb1) * 
        (nuij - 1 - cij + 2*iOR*iProb2)/(nuij*(iProb2 - iProb12));
      }
      //for (int k =0; k < pNTtimeDep; ++k)
      //{
        iGrad[3*Q + pNT + 0 + 1] = -0.5 * XNTtimeDepNow[i] * iProb1 * (1 - iProb1) * 
          (nuij - 1 - cij + 2*iOR*iProb2)/(nuij*(iProb2 - iProb12));
        // only one time-dep var
      //}
      for (int k =0; k < pT; ++k)
      {
      iGrad[3*Q + pNT + pNTtimeDep + k + 1] = ((XT(i,k)*iProb2*(1 - iProb2)/(iProb2 - iProb12))) * 
        (1 -  0.5*(nuij - 1 - cij + 2*iOR*iProb1)/nuij);
      }
      //for (int k =0; k < pTtimeDep; ++k)
      //{
        iGrad[3*Q + pNT + pNTtimeDep + pT + 0 + 1] = ((XTtimeDepNow[i]*iProb2*(1 - iProb2)/(iProb2 - iProb12))) * 
          (1 -  0.5*(nuij - 1 - cij + 2*iOR*iProb1)/nuij);
      // only one time-dep var
      //}
      for (int k =0; k < pOR; ++k)
      {
      iGrad[3*Q + pNT + pNTtimeDep + pT + pTtimeDep + k + 1] =  -XOR(i,k)*iOR*
        (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
          (2*(iOR - 1)*(iOR - 1)*(iProb2 - iProb12));
      }
      //for (int k =0; k < pORtimeDep; ++k)
      //{
        iGrad[3*Q + pNT + pNTtimeDep + pT + pTtimeDep + pOR + 0 + 1] =  -XORtimeDepNow[i]*iOR*
          (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
            (2*(iOR - 1)*(iOR - 1)*(iProb2 - iProb12));
      // only one time-dep var
      //}
      }
      if (YNT(i,j)==0 && YT(i,j)==0) {
      for (int q=0; q < Q; ++q)
        {
        b = TimeBase(j,q);
        iGrad[1 + q] =  ((b*iProb1*(1 - iProb1)/(1 - iProb1 - iProb2 + iProb12))) *(0.5*(nuij - 1 - cij + 2*iOR*iProb2)/nuij - 1);
        iGrad[Q + 1 + q] = ((b*iProb2*(1 - iProb2)/(1 - iProb1 - iProb2 + iProb12))) *(0.5*(nuij - 1 - cij + 2*iOR*iProb1)/nuij - 1);
                iGrad[2*Q + 1 + q] =  b*iOR*
                  (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
                    (2*(iOR - 1)*(iOR - 1)*(1 - iProb1 - iProb2 + iProb12));
        }
      for (int k =0; k < pNT; ++k)
      {
      iGrad[3*Q + k + 1] = ((XNT(i,k)*iProb1*(1 - iProb1)/(1 - iProb1 - iProb2 + iProb12))) *
        (0.5*(nuij - 1 - cij + 2*iOR*iProb2)/nuij - 1);
      }
      //for (int k =0; k < pNTtimeDep; ++k)
      //{
        iGrad[3*Q + pNT + 0 + 1] = ((XNTtimeDepNow[i]*iProb1*(1 - iProb1)/(1 - iProb1 - iProb2 + iProb12))) *
          (0.5*(nuij - 1 - cij + 2*iOR*iProb2)/nuij - 1);
        // only one time-dep var
      //}
      for (int k =0; k < pT; ++k)
      {
      iGrad[3*Q + pNT + pNTtimeDep + k + 1] = ((XT(i,k)*iProb2*(1 - iProb2)/(1 - iProb1 - iProb2 + iProb12))) * 
        (0.5*(nuij - 1 - cij + 2*iOR*iProb1)/nuij - 1);
      }
      //for (int k =0; k < pTtimeDep; ++k)
      //{
        iGrad[3*Q + pNT + pNTtimeDep + pT + 0 + 1] = ((XTtimeDepNow[i]*iProb2*(1 - iProb2)/
          (1 - iProb1 - iProb2 + iProb12))) * 
          (0.5*(nuij - 1 - cij + 2*iOR*iProb1)/nuij - 1);
        // only one time-dep var
      //}
      for (int k =0; k < pOR; ++k)
      {
      iGrad[3*Q + pNT + pNTtimeDep + pT + pTtimeDep + k + 1] =  XOR(i,k)*iOR*
        (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
          (2*(iOR - 1)*(iOR - 1)*(1 - iProb1 - iProb2 + iProb12));
      }
      //for (int k =0; k < pORtimeDep; ++k)
      //{
        iGrad[3*Q + pNT + pNTtimeDep + pT + pTtimeDep + pOR + 0 + 1] =  XORtimeDepNow[i]*iOR*
          (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
            (2*(iOR - 1)*(iOR - 1)*(1 - iProb1 - iProb2 + iProb12));
      // only one time-dep var
      //}
      }
      }}}
      // Rcpp::Rcout << "iContrib "<< iContrib << std::endl;
      // if (arma::is_finite(iContrib)) {} else {
      //   Rcpp::Rcout << "iContrib "<< iContrib << std::endl;
      //   Rcpp::Rcout << "j =  "<< j << std::endl;
      //   Rcpp::Rcout << "i =  "<< i << std::endl;
      //   Rcpp::Rcout << "iProb1  "<< iProb1 << std::endl;
      //   Rcpp::Rcout << "iProb2  "<< iProb2 << std::endl;
      //   Rcpp::Rcout << "iProb12  "<< iProb12 << std::endl;
      //   Rcpp::Rcout << "iOR  "<< iOR << std::endl;}
//        Rcpp::Rcout << "iORCOND  "<< iORcond << std::endl;}
// if (i==997)
// {
//   Rcpp::Rcout << "i:"<< i << std::endl;
//   Rcpp::Rcout << "ExpXBetaNTtimeDep:"<< ExpXBetaNTtimeDep[i] << std::endl;
//   Rcpp::Rcout << "ExpXBetaTtimeDep:"<< ExpXBetaTtimeDep[i] << std::endl;
//   Rcpp::Rcout << "ExpXBetaORtimeDep:"<< ExpXBetaORtimeDep[i] << std::endl;
//   Rcpp::Rcout << "iGrad:"<< iGrad << std::endl;
//   Rcpp::Rcout << "Location:"<< 3*Q + pNT + pNTtimeDep + pT + 0 + 1 << std::endl;
// }
      Grad += iGrad;
    }}
  for(int q=0; q<Q; ++q)
  {
    Grad[1 + q] -= penaltermNT[q];
    Grad[Q + 1 + q] -= penaltermT[q];
    Grad[2*Q + 1 + q] -= penaltermOR[q];
  }
  return(-Grad);
     }
