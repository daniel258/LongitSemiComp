#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec GradPenalLogLik(arma::vec param, arma::mat YNT, arma::mat riskNT, arma::mat riskT, arma::mat YT, 
                         arma::mat XNT, arma::mat XT, arma::mat XOR, arma::mat XinterMat, 
                         arma::mat TimeBase, arma::mat TimePen, arma::vec lambda, double epsOR)
{
  int n = YT.n_rows;
  int K = YT.n_cols;
  int pNT = XNT.n_cols;
  int pT = XT.n_cols;
  int pOR = XOR.n_cols;
  int pInter = XinterMat.n_cols;
  int J = TimeBase.n_cols;
  // J is the the number of B-splines (number of rows in TimeBase should be K)
  // Decomposing "param" to individual paramters as in PenalLogLik
  arma::vec alphaNT = param.subvec(0, J - 1);
  arma::vec alphaT = param.subvec(J, 2*J -1 );
  arma::vec alphaOR = param.subvec(2*J, 3*J - 1);
  arma::vec penaltermNT = 2 * lambda[0] * TimePen * alphaNT;
  arma::vec penaltermT = 2 * lambda[1] * TimePen * alphaT;
  arma::vec penaltermOR = 2 * lambda[2]  * TimePen * alphaOR;
  arma::vec betaNT = param.subvec(3*J,3*J + pNT - 1);
  arma::vec betaT = param.subvec(3*J + pNT, 3*J + pNT + pT - 1);
  arma::vec betaOR = param.subvec(3*J + pNT + pT, 3*J + pNT + pT + pOR - 1);
  arma::vec betaInt = param.subvec(3*J + pNT + pT + pOR, 3*J +  pNT + pT + pOR + pInter - 1);
  arma::vec Grad(3*J + pNT + pT + pOR + pInter);
  arma::vec iGrad(3*J + pNT + pT + pOR + pInter);
  Grad.fill(0);
  iGrad.fill(0);
  double iProbTafterNT = 0;
  double iProb1 = 0;
  double iProb2 = 0;
  double iOR = 0;
  double iProb12 = 0;
  double cij = 0;
  double nuij = 0;
  double b = 0;
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
      iGrad.fill(0);
      if (riskT(i,j)==0) {
        iGrad.fill(0);
        //   no contribution to the gradient
      } else {
        if (riskNT(i,j)==0) {
          iProbTafterNT = (ExpAlphaTnow*ExpXBetaT[i]*ExpXBetaInt[i])/
            (1 + (ExpAlphaTnow*ExpXBetaT[i]*ExpXBetaInt[i]));
          if(YT(i,j)==1) {
            for (int q=0; q < J; ++q)
            {
              b = TimeBase(j,q);
              iGrad[J + q] = b*(1-iProbTafterNT);
            }
            for (int k =0; k < pT; ++k)
            {
              iGrad[3*J + pNT + k] = XT(i,k)*(1-iProbTafterNT);
            }
            for (int k =0; k < pInter; ++k)
            {
              iGrad[3*J + pNT + pT + pOR + k] = XinterMat(i,k)*(1-iProbTafterNT);
            }}
          else {
            for (int q = 0; q < J; ++q)
            {
              b = TimeBase(j,q);
              iGrad[J + q] = -b*iProbTafterNT;
            }
            for (int k = 0; k < pT; ++k)
            {
              iGrad[3*J + pNT + k] = -XT(i,k)*iProbTafterNT;
            }
            for (int k = 0; k < pInter; ++k)
            {
              iGrad[3*J + pNT + pT + pOR + k] = -XinterMat(i,k)*iProbTafterNT;
            }}
        } else {
          iProb1 = ExpAlphaNTnow*ExpXBetaNT[i]/(1 + (ExpAlphaNTnow*ExpXBetaNT[i]));
          iProb2 = ExpAlphaTnow*ExpXBetaT[i]/(1 + ExpAlphaTnow*ExpXBetaT[i]);
          iOR = ExpAlphaORnow*ExpXBetaOR[i];
          if ((iOR > 1 - epsOR) & (iOR < 1 + epsOR))
          {
            iProb12 = iProb1*iProb2;
            if (YNT(i,j)==1 && YT(i,j)==1) {
              for (int q=0; q < J; ++q)
              {
                b = TimeBase(j,q);
                iGrad[q] = b*(1 - iProb1);
                iGrad[J + q] = b*(1 - iProb2);
                iGrad[2*J + q] = 0;
              }
              for (int k =0; k < pNT; ++k)
              {
                iGrad[3*J + k] = XNT(i,k)*(1 - iProb1);
              }
              for (int k =0; k < pT; ++k)
              {
                iGrad[3*J + pNT + k] = XT(i,k)*(1 - iProb2);
              }
              for (int k =0; k < pOR; ++k)
              {
                iGrad[3*J + pNT + pT + k] = 0;
              }}
            if (YNT(i,j)==1 && YT(i,j)==0) {
              for (int q=0; q < J; ++q)
              {
                b = TimeBase(j,q);
                iGrad[q] =  b*(1 - iProb1);
                iGrad[J + q] = -b*iProb2;
                iGrad[2*J + q] = 0;
              }
              for (int k = 0; k < pNT; ++k)
              {
                iGrad[3*J + k] = XNT(i,k)*(1 - iProb1);
              }
              for (int k = 0; k < pT; ++k)
              {
                iGrad[3*J + pNT + k] = -XT(i,k)*iProb2;
              }
              for (int k = 0; k < pOR; ++k)
              {
                iGrad[3*J + pNT + pT + k] =  0;
              }}
            if (YNT(i,j)==0 && YT(i,j)==1) {
              for (int q = 0; q < J; ++q)
              {
                b = TimeBase(j,q);
                iGrad[q] =  -b*iProb1;
                iGrad[J + q] = b*(1 - iProb2);
                iGrad[2*J + q] = 0;
              }
              for (int k = 0; k < pNT; ++k)
              {
                iGrad[3*J + k] = -XNT(i,k)*iProb1;
              }
              for (int k =0; k < pT; ++k)
              {
                iGrad[3*J + pNT + k] = XT(i,k)*(1 - iProb2);
              }
              for (int k =0; k < pOR; ++k)
              {
                iGrad[3*J + pNT + pT + k] =  0;
              }}
            if (YNT(i,j)==0 && YT(i,j)==0) {
              for (int q=0; q < J; ++q)
              {
                b = TimeBase(j,q);
                iGrad[q] =  -b*iProb1;
                iGrad[J + q] = -b*iProb2;
                iGrad[2*J + q] =  0;
              }
              for (int k = 0; k < pNT; ++k)
              {
                iGrad[3*J + k] = -XNT(i,k)*iProb1;
              }
              for (int k = 0; k < pT; ++k)
              {
                iGrad[3*J + pNT + k] = -XT(i,k)*iProb2;
              }
              for (int k =0; k < pOR; ++k)
              {
                iGrad[3*J + pNT + pT + k] =  0;
              }}
          } else {
            cij = (iProb1 + iProb2)*(iOR - 1);
            nuij = sqrt(pow(1 + cij, 2.0) - 4*iOR*(iOR - 1)*iProb1*iProb2);
            iProb12 = (1 + cij - nuij) / (2 * (iOR - 1));
            if (YNT(i,j)==1 && YT(i,j)==1) {
              for (int q=0; q < J; ++q)
              {
                b = TimeBase(j,q);
                iGrad[q] = 0.5*b*iProb1*(1 - iProb1) * (nuij - 1 - cij + 2*iOR*iProb2)/(nuij*iProb12);
                iGrad[J + q] = 0.5*b*iProb2*(1 - iProb2) * (nuij - 1 - cij + 2*iOR*iProb1)/(nuij*iProb12);
                iGrad[2*J + q] = b*iOR*
                  (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
                    (2*(iOR - 1)*(iOR - 1)*iProb12);
              }
              for (int k = 0; k < pNT; ++k)
              {
                iGrad[3*J + k] = 0.5*XNT(i,k)*iProb1*(1 - iProb1) * (nuij - 1 - cij + 2*iOR*iProb2)/(nuij*iProb12);
              }
              for (int k = 0; k < pT; ++k)
              {
                iGrad[3*J + pNT + k] = 0.5*XT(i,k)*iProb2*(1 - iProb2) * (nuij - 1 - cij + 2*iOR*iProb1)/(nuij*iProb12);
              }
              for (int k = 0; k < pOR; ++k)
              {
                iGrad[3*J + pNT + pT + k] = XOR(i,k)*iOR*
                  (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
                    (2*(iOR - 1)*(iOR - 1)*iProb12);
              }}
            if (YNT(i,j)==1 && YT(i,j)==0) {
              for (int q = 0; q < J; ++q)
              {
                b = TimeBase(j,q);
                iGrad[q] =  ((b*iProb1*(1 - iProb1)/(iProb1 - iProb12))) *(1 -  0.5*(nuij - 1 - cij + 2*iOR*iProb2)/nuij);
                iGrad[J + q] = -0.5*b*iProb2*(1 - iProb2) * (nuij - 1 - cij + 2*iOR*iProb1)/(nuij*(iProb1 - iProb12));
                iGrad[2*J + q] = -b*iOR*
                  (((iOR - 1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
                    (2*(iOR - 1)*(iOR - 1)*(iProb1 - iProb12));
              }
              for (int k = 0; k < pNT; ++k)
              {
                iGrad[3*J + k] = ((XNT(i,k)*iProb1*(1 - iProb1)/(iProb1 - iProb12))) *(1 -  0.5*(nuij - 1 - cij + 2*iOR*iProb2)/nuij);
              }
              for (int k =0; k < pT; ++k)
              {
                iGrad[3*J + pNT + k] = -0.5*XT(i,k)*iProb2*(1 - iProb2) * (nuij - 1 - cij + 2*iOR*iProb1)/(nuij*(iProb1 - iProb12));
              }
              for (int k =0; k < pOR; ++k)
              {
                iGrad[3*J + pNT + pT + k] =  -XOR(i,k)*iOR*
                  (((iOR - 1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
                    (2*(iOR - 1)*(iOR - 1)*(iProb1 - iProb12));
              }}
            if (YNT(i,j)==0 && YT(i,j)==1) {
              for (int q=0; q < J; ++q)
              {
                b = TimeBase(j,q);
                iGrad[q] =  -0.5*b*iProb1*(1 - iProb1) * (nuij - 1 - cij + 2*iOR*iProb2)/(nuij*(iProb2 - iProb12));
                iGrad[J + q] = ((b*iProb2*(1 - iProb2)/(iProb2 - iProb12))) *(1 -  0.5*(nuij - 1 - cij + 2*iOR*iProb1)/nuij);
                iGrad[2*J + q] = -b*iOR*
                  (((iOR - 1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
                    (2*(iOR - 1)*(iOR - 1)*(iProb2 - iProb12));
              }
              for (int k =0; k < pNT; ++k)
              {
                iGrad[3*J + k] = -0.5*XNT(i,k)*iProb1*(1 - iProb1) * (nuij - 1 - cij + 2*iOR*iProb2)/(nuij*(iProb2 - iProb12));
              }
              for (int k = 0; k < pT; ++k)
              {
                iGrad[3*J + pNT + k] = ((XT(i,k)*iProb2*(1 - iProb2)/(iProb2 - iProb12))) *(1 -  0.5*(nuij - 1 - cij + 2*iOR*iProb1)/nuij);
              }
              for (int k = 0; k < pOR; ++k)
              {
                iGrad[3*J + pNT + pT + k] =  -XOR(i,k)*iOR*
                  (((iOR - 1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
                    (2*(iOR - 1)*(iOR - 1)*(iProb2 - iProb12));
              }}
            if (YNT(i,j)==0 && YT(i,j)==0) {
              for (int q=0; q < J; ++q)
              {
                b = TimeBase(j,q);
                iGrad[q] =  ((b*iProb1*(1 - iProb1)/(1 - iProb1 - iProb2 + iProb12))) *(0.5*(nuij - 1 - cij + 2*iOR*iProb2)/nuij - 1);
                iGrad[J + q] = ((b*iProb2*(1 - iProb2)/(1 - iProb1 - iProb2 + iProb12))) *(0.5*(nuij - 1 - cij + 2*iOR*iProb1)/nuij - 1);
                iGrad[2*J + q] =  b*iOR*
                  (((iOR - 1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
                    (2*(iOR - 1)*(iOR - 1)*(1 - iProb1 - iProb2 + iProb12));
              }
              for (int k = 0; k < pNT; ++k)
              {
                iGrad[3*J + k] = ((XNT(i,k)*iProb1*(1 - iProb1)/(1 - iProb1 - iProb2 + iProb12))) *(0.5*(nuij - 1 - cij + 2*iOR*iProb2)/nuij - 1);
              }
              for (int k =0; k < pT; ++k)
              {
                iGrad[3*J + pNT + k] = ((XT(i,k)*iProb2*(1 - iProb2)/(1 - iProb1 - iProb2 + iProb12))) *(0.5*(nuij - 1 - cij + 2*iOR*iProb1)/nuij - 1);
              }
              for (int k = 0; k < pOR; ++k)
              {
                iGrad[3*J + pNT + pT + k] =  XOR(i,k)*iOR*
                  (((iOR - 1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
                    (2*(iOR - 1)*(iOR - 1)*(1 - iProb1 - iProb2 + iProb12));
              }}
          }}}
      Grad += iGrad;
    }}
  for(int q=0; q<J; ++q)
  {
    Grad[q] -= penaltermNT[q];
    Grad[J + q] -= penaltermT[q];
    Grad[2*J + q] -= penaltermOR[q];
  }
  return(-Grad);
}
