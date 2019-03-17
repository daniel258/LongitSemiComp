#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
//' @export
// [[Rcpp::export]]
arma::vec GradPenalLogLik(arma::vec param, arma::mat YNT, arma::mat riskNT, arma::mat riskT, arma::mat YT, 
                          arma::mat XNT, arma::mat XT, arma::mat XOR, arma::mat XinterMat, 
                          arma::mat TimeBase, arma::mat TimePen, arma::vec lambda, double epsOR)
     {
  int n = YT.n_rows;
  int J = YT.n_cols;
  int pNT = XNT.n_cols;
  int pT = XT.n_cols;
  int pOR = XOR.n_cols;
  int pGammaInter = XinterMat.n_cols;
  int Q = TimeBase.n_cols;
  // Q is the the number of B-splines (number of rows in TimeBase should be J)
  // Decomposing "param" to individual paramters as in PenalLogLik
  double gamma = param[0];
  arma::vec alphaNT = param.subvec(1, Q);
  arma::vec alphaT = param.subvec(Q + 1, 2*Q);
  arma::vec alphaOR = param.subvec(2*Q + 1, 3*Q);
  arma::vec penaltermNT = 2 * lambda[0] * TimePen * alphaNT;
  arma::vec penaltermT = 2 * lambda[1] * TimePen * alphaT;
  arma::vec penaltermOR = 2 * lambda[2]  * TimePen * alphaOR;
  arma::vec betaNT = param.subvec(3*Q + 1,3*Q + pNT);
  arma::vec betaT = param.subvec(3*Q + pNT + 1,3*Q + pNT + pT);
  arma::vec betaOR = param.subvec(3*Q + pNT + pT + 1,3*Q + pNT + pT + pOR);
  arma::vec gammaInt = param.subvec(3*Q + pNT + pT + pOR + 1, 3*Q +  pNT + pT + pOR + pGammaInter);
  arma::vec Grad(1 + 3*Q + pNT + pT + pOR + pGammaInter);
  arma::vec iGrad(1 + 3*Q + pNT + pT + pOR + pGammaInter);
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
  arma::vec ExpXBetaNT = exp(XNT * betaNT);
  arma::vec ExpXBetaT = exp(XT * betaT);
  arma::vec ExpXBetaOR = exp(XOR * betaOR);
  arma::vec ExpXGammaInt = exp(XinterMat * gammaInt);
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
      iGrad.fill(0);
        if (riskT(i,j)==0) {
          iGrad.fill(0);
          //   no contribution to the gradient
        } else {
          if (riskNT(i,j)==0) {
            iProbTafterNT = (ExpAlphaTnow*ExpXBetaT[i]*exp(gamma)*ExpXGammaInt[i])/
              (1 + (ExpAlphaTnow*ExpXBetaT[i]*exp(gamma)*ExpXGammaInt[i]));
            if(YT(i,j)==1) {
              iGrad[0] = 1-iProbTafterNT;
              for (int q=0; q < Q; ++q)
              {
                b = TimeBase(j,q);
                iGrad[1 + Q + q ] = b*(1-iProbTafterNT);
              }
              for (int k =0; k < pT; ++k)
              {
              iGrad[3*Q + pNT + k + 1] = XT(i,k)*(1-iProbTafterNT);
              }
              for (int k =0; k < pGammaInter; ++k)
              {
                iGrad[3*Q + pNT + pT + pOR + k + 1] = XinterMat(i,k)*(1-iProbTafterNT);
              }}
            else {
              iGrad[0] = -iProbTafterNT;
              for (int q=0; q < Q; ++q)
              {
                b = TimeBase(j,q);
                iGrad[1 + Q + q] = -b*iProbTafterNT;
              }
              for (int k =0; k < pT; ++k)
              {
                iGrad[3*Q + pNT + k + 1] = -XT(i,k)*iProbTafterNT;
              }
              for (int k =0; k < pGammaInter; ++k)
              {
                iGrad[3*Q + pNT + pT + pOR + k + 1] = -XinterMat(i,k)*iProbTafterNT;
              }}
          } else {
            iProb1 = ExpAlphaNTnow*ExpXBetaNT[i]/(1 + (ExpAlphaNTnow*ExpXBetaNT[i]));
            iProb2 = ExpAlphaTnow*ExpXBetaT[i]/(1 + ExpAlphaTnow*ExpXBetaT[i]);
            iOR = ExpAlphaORnow*ExpXBetaOR[i];
            // Rcpp::Rcout << "iProb1  "<< iProb1 << std::endl;
            // Rcpp::Rcout << "iProb2 "<< iProb2 << std::endl;
            // Rcpp::Rcout << "iOR "<< iOR << std::endl;
            if ((iOR > 1 - epsOR) & (iOR < 1 + epsOR))
            {
              iProb12 = iProb1*iProb2;
              // Rcpp::Rcout << "Small OR:  "<< iOR << std::endl;
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
                 for (int k =0; k < pT; ++k)
                 {
                   iGrad[3*Q + pNT + k + 1] = XT(i,k)*(1 - iProb2);
                 }
                 for (int k =0; k < pOR; ++k)
                 {
                   iGrad[3*Q + pNT + pT + k + 1] = 0;
                 }}
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
                 for (int k =0; k < pT; ++k)
                 {
                   iGrad[3*Q + pNT + k + 1] = -XT(i,k)*iProb2;
                 }
                 for (int k = 0; k < pOR; ++k)
                 {
                   iGrad[3*Q + pNT + pT + k + 1] =  0;
                 }}
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
                 for (int k =0; k < pT; ++k)
                 {
                   iGrad[3*Q + pNT + k + 1] = XT(i,k)*(1 - iProb2);
                 }
                 for (int k =0; k < pOR; ++k)
                 {
                   iGrad[3*Q + pNT + pT + k + 1] =  0;
                 }}
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
                   for (int k =0; k < pT; ++k)
                 {
                   iGrad[3*Q + pNT + k + 1] = -XT(i,k)*iProb2;
                  }
                   for (int k =0; k < pOR; ++k)
                   {
                   iGrad[3*Q + pNT + pT + k + 1] =  0;
                  }}
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
                iGrad[3*Q + k + 1] = 0.5*XNT(i,k)*iProb1*(1 - iProb1) * (nuij - 1 - cij + 2*iOR*iProb2)/(nuij*iProb12);
              }
              for (int k =0; k < pT; ++k)
              {
                iGrad[3*Q + pNT + k + 1] = 0.5*XT(i,k)*iProb2*(1 - iProb2) * (nuij - 1 - cij + 2*iOR*iProb1)/(nuij*iProb12);
              }
              for (int k =0; k < pOR; ++k)
              {
                iGrad[3*Q + pNT + pT + k + 1] = XOR(i,k)*iOR*
                  (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
                  (2*(iOR - 1)*(iOR - 1)*iProb12);
              }}
            if (YNT(i,j)==1 && YT(i,j)==0) {
              for (int q=0; q < Q; ++q)
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
      iGrad[3*Q + k + 1] = ((XNT(i,k)*iProb1*(1 - iProb1)/(iProb1 - iProb12))) *(1 -  0.5*(nuij - 1 - cij + 2*iOR*iProb2)/nuij);
      }
      for (int k =0; k < pT; ++k)
      {
      iGrad[3*Q + pNT + k + 1] = -0.5*XT(i,k)*iProb2*(1 - iProb2) * (nuij - 1 - cij + 2*iOR*iProb1)/(nuij*(iProb1 - iProb12));
      }
      for (int k =0; k < pOR; ++k)
      {
      iGrad[3*Q + pNT + pT + k + 1] =  -XOR(i,k)*iOR*
        (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
         (2*(iOR - 1)*(iOR - 1)*(iProb1 - iProb12));
      }}
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
      iGrad[3*Q + k + 1] = -0.5*XNT(i,k)*iProb1*(1 - iProb1) * (nuij - 1 - cij + 2*iOR*iProb2)/(nuij*(iProb2 - iProb12));
      }
      for (int k =0; k < pT; ++k)
      {
      iGrad[3*Q + pNT + k + 1] = ((XT(i,k)*iProb2*(1 - iProb2)/(iProb2 - iProb12))) *(1 -  0.5*(nuij - 1 - cij + 2*iOR*iProb1)/nuij);
      }
      for (int k =0; k < pOR; ++k)
      {
      iGrad[3*Q + pNT + pT + k + 1] =  -XOR(i,k)*iOR*
        (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
          (2*(iOR - 1)*(iOR - 1)*(iProb2 - iProb12));
      }}
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
      iGrad[3*Q + k + 1] = ((XNT(i,k)*iProb1*(1 - iProb1)/(1 - iProb1 - iProb2 + iProb12))) *(0.5*(nuij - 1 - cij + 2*iOR*iProb2)/nuij - 1);
      }
      for (int k =0; k < pT; ++k)
      {
      iGrad[3*Q + pNT + k + 1] = ((XT(i,k)*iProb2*(1 - iProb2)/(1 - iProb1 - iProb2 + iProb12))) *(0.5*(nuij - 1 - cij + 2*iOR*iProb1)/nuij - 1);
      }
      for (int k =0; k < pOR; ++k)
      {
      iGrad[3*Q + pNT + pT + k + 1] =  XOR(i,k)*iOR*
        (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
          (2*(iOR - 1)*(iOR - 1)*(1 - iProb1 - iProb2 + iProb12));
      }}
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
