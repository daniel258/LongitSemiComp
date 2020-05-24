#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec GradParamLogLikWithInter(arma::vec param, arma::mat YT, arma::mat YNT, arma::mat riskNT, arma::mat riskT, 
                      arma::mat XNT, arma::mat XT, arma::mat XOR, arma::mat XinterMat, double epsOR)
     {
  int n = YT.n_rows;
  int J = YT.n_cols;
  int pNT = XNT.n_cols;
  int pT = XT.n_cols;
  int pOR = XOR.n_cols;
  int pGammaInter = XinterMat.n_cols;
  // Decomposing "param" to individual paramters as in ParamLogLik
  double gamma = param[0];
  arma::vec alphaNT = param.subvec(1, J);
  arma::vec alphaT = param.subvec(J + 1, 2*J);
  arma::vec alphaOR = param.subvec(2*J + 1, 3*J);
  arma::vec betaNT = param.subvec(3*J + 1, 3*J + pNT);
  arma::vec betaT = param.subvec(3*J + pNT + 1, 3*J + pNT + pT);
  arma::vec betaOR = param.subvec(3*J + pNT + pT + 1, 3*J + pNT + pT + pOR);
  arma::vec gammaInt = param.subvec(3*J + pNT + pT + pOR + 1, 3*J +  pNT + pT + pOR + pGammaInter);
  arma::vec Grad(1 + 3*J + pNT + pT + pOR + pGammaInter);
  arma::vec iGrad(1 + 3*J + pNT + pT + pOR + pGammaInter);
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
  //double b = 0;
  double ExpAlphaNTnow = 0;
  double ExpAlphaTnow = 0;
  double ExpAlphaORnow = 0;
  arma::vec ExpXBetaNT = exp(XNT * betaNT);
  arma::vec ExpXBetaT = exp(XT * betaT);
  arma::vec ExpXBetaOR = exp(XOR * betaOR);
  arma::vec ExpXGammaInt = exp(XinterMat * gammaInt);
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
              iGrad[1 + J + j] = (1-iProbTafterNT);
              for (int k =0; k < pT; ++k)
              {
              iGrad[3*J + pNT + k + 1] = XT(i,k)*(1-iProbTafterNT);
              }
              for (int k =0; k < pGammaInter; ++k)
              {
                iGrad[3*J + pNT + pT + pOR + k + 1] = XinterMat(i,k)*(1-iProbTafterNT);
              }}
            else {
              iGrad[0] = -iProbTafterNT;
              iGrad[1 + J + j] = -iProbTafterNT;
              for (int k =0; k < pT; ++k)
              {
                iGrad[3*J + pNT + k + 1] = -XT(i,k)*iProbTafterNT;
              }
              for (int k =0; k < pGammaInter; ++k)
              {
                iGrad[3*J + pNT + pT + pOR + k + 1] = -XinterMat(i,k)*iProbTafterNT;
              }}
          } else {
            iProb1 = ExpAlphaNTnow*ExpXBetaNT[i]/(1 + (ExpAlphaNTnow*ExpXBetaNT[i]));
            iProb2 = ExpAlphaTnow*ExpXBetaT[i]/(1 + ExpAlphaTnow*ExpXBetaT[i]);
            iOR = ExpAlphaORnow*ExpXBetaOR[i];
            if ((iOR > 1 - epsOR) & (iOR < 1 + epsOR))
            {
             iProb12 = iProb1*iProb2;
             if (YNT(i,j)==1 && YT(i,j)==1) {
                iGrad[1 + j] = (1 - iProb1);
                iGrad[J + 1 + j] = (1 - iProb2);
                iGrad[2*J + 1 + j] = 0;
                for (int k =0; k < pNT; ++k)
                {
                  iGrad[3*J + k + 1] = XNT(i,k)*(1 - iProb1);
                }
                for (int k =0; k < pT; ++k)
                {
                  iGrad[3*J + pNT + k + 1] = XT(i,k)*(1 - iProb2);
                }
                for (int k =0; k < pOR; ++k)
                {
                  iGrad[3*J + pNT + pT + k + 1] = 0;
                }}
              if (YNT(i,j)==1 && YT(i,j)==0) {
                iGrad[1 + j] =  (1 - iProb1);
                iGrad[J + 1 + j] = -iProb2;
                iGrad[2*J + 1 + j] = 0;
                for (int k =0; k < pNT; ++k)
                {
                  iGrad[3*J + k + 1] = XNT(i,k)*(1 - iProb1);
                }
                for (int k =0; k < pT; ++k)
                {
                  iGrad[3*J + pNT + k + 1] = -XT(i,k)*iProb2;
                }
                for (int k = 0; k < pOR; ++k)
                {
                  iGrad[3*J + pNT + pT + k + 1] =  0;
                }}
               if (YNT(i,j)==0 && YT(i,j)==1) {
                iGrad[1 + j] =  -iProb1;
                iGrad[J + 1 + j] = (1 - iProb2);
                iGrad[2*J + 1 + j] = 0;
                for (int k =0; k < pNT; ++k)
                {
                  iGrad[3*J + k + 1] = -XNT(i,k)*iProb1;
                }
                for (int k =0; k < pT; ++k)
                {
                  iGrad[3*J + pNT + k + 1] = XT(i,k)*(1 - iProb2);
                }
                for (int k =0; k < pOR; ++k)
                {
                  iGrad[3*J + pNT + pT + k + 1] =  0;
                }}
               if (YNT(i,j)==0 && YT(i,j)==0) {
                 iGrad[1 + j] =  -iProb1;
                 iGrad[J + 1 + j] = -iProb2;
                 iGrad[2*J + 1 + j] =  0;
                 for (int k =0; k < pNT; ++k)
                 {
                  iGrad[3*J + k + 1] = -XNT(i,k)*iProb1;
                 }
                 for (int k =0; k < pT; ++k)
                 {
                  iGrad[3*J + pNT + k + 1] = -XT(i,k)*iProb2;
                 }
                 for (int k =0; k < pOR; ++k)
                 {
                  iGrad[3*J + pNT + pT + k + 1] =  0;
                 }}
            } else {
              cij = (iProb1 + iProb2)*(iOR - 1);
              nuij = sqrt(pow(1 + cij, 2.0) - 4*iOR*(iOR - 1)*iProb1*iProb2);
              iProb12 = (1 + cij - nuij) / (2 * (iOR - 1));

            // Rcpp::Rcout << "iProb12  "<< iProb12 << std::endl;
      if (YNT(i,j)==1 && YT(i,j)==1) {
              iGrad[1 + j] = 0.5*iProb1*(1 - iProb1) * (nuij - 1 - cij + 2*iOR*iProb2)/(nuij*iProb12);
              iGrad[J + 1 + j] = 0.5*iProb2*(1 - iProb2) * (nuij - 1 - cij + 2*iOR*iProb1)/(nuij*iProb12);
              iGrad[2*J + 1 + j] = iOR*
                (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
                    (2*(iOR - 1)*(iOR - 1)*iProb12);
              for (int k =0; k < pNT; ++k)
              {
                iGrad[3*J + k + 1] = 0.5*XNT(i,k)*iProb1*(1 - iProb1) * (nuij - 1 - cij + 2*iOR*iProb2)/(nuij*iProb12);
              }
              for (int k =0; k < pT; ++k)
              {
                iGrad[3*J + pNT + k + 1] = 0.5*XT(i,k)*iProb2*(1 - iProb2) * (nuij - 1 - cij + 2*iOR*iProb1)/(nuij*iProb12);
              }
              for (int k =0; k < pOR; ++k)
              {
                iGrad[3*J + pNT + pT + k + 1] = XOR(i,k)*iOR*
                  (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
                  (2*(iOR - 1)*(iOR - 1)*iProb12);
              }}
      if (YNT(i,j)==1 && YT(i,j)==0) {
              iGrad[1 + j] =  ((iProb1*(1 - iProb1)/(iProb1 - iProb12))) *(1 -  0.5*(nuij - 1 - cij + 2*iOR*iProb2)/nuij);
              iGrad[J + 1 + j] = -0.5*iProb2*(1 - iProb2) * (nuij - 1 - cij + 2*iOR*iProb1)/(nuij*(iProb1 - iProb12));
              iGrad[2*J + 1 + j] = -iOR*
                  (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
                    (2*(iOR - 1)*(iOR - 1)*(iProb1 - iProb12));
            for (int k = 0; k < pNT; ++k)
            {
            iGrad[3*J + k + 1] = ((XNT(i,k)*iProb1*(1 - iProb1)/(iProb1 - iProb12))) *(1 -  0.5*(nuij - 1 - cij + 2*iOR*iProb2)/nuij);
            }
            for (int k =0; k < pT; ++k)
            {
            iGrad[3*J + pNT + k + 1] = -0.5*XT(i,k)*iProb2*(1 - iProb2) * (nuij - 1 - cij + 2*iOR*iProb1)/(nuij*(iProb1 - iProb12));
            }
            for (int k =0; k < pOR; ++k)
            {
            iGrad[3*J + pNT + pT + k + 1] =  -XOR(i,k)*iOR*
              (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
              (2*(iOR - 1)*(iOR - 1)*(iProb1 - iProb12));
      }}
      if (YNT(i,j)==0 && YT(i,j)==1) {
        iGrad[1 + j] =  -0.5*iProb1*(1 - iProb1) * (nuij - 1 - cij + 2*iOR*iProb2)/(nuij*(iProb2 - iProb12));
        iGrad[J + 1 + j] = ((iProb2*(1 - iProb2)/(iProb2 - iProb12))) *(1 -  0.5*(nuij - 1 - cij + 2*iOR*iProb1)/nuij);
        iGrad[2*J + 1 + j] = -iOR*
          (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
              (2*(iOR - 1)*(iOR - 1)*(iProb2 - iProb12));
      for (int k =0; k < pNT; ++k)
      {
        iGrad[3*J + k + 1] = -0.5*XNT(i,k)*iProb1*(1 - iProb1) * (nuij - 1 - cij + 2*iOR*iProb2)/(nuij*(iProb2 - iProb12));
      }
      for (int k =0; k < pT; ++k)
      {
        iGrad[3*J + pNT + k + 1] = ((XT(i,k)*iProb2*(1 - iProb2)/(iProb2 - iProb12))) *(1 -  0.5*(nuij - 1 - cij + 2*iOR*iProb1)/nuij);
      }
      for (int k =0; k < pOR; ++k)
      {
        iGrad[3*J + pNT + pT + k + 1] =  -XOR(i,k)*iOR*
        (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
          (2*(iOR - 1)*(iOR - 1)*(iProb2 - iProb12));
      }}
      if (YNT(i,j)==0 && YT(i,j)==0) {
        iGrad[1 + j] =  ((iProb1*(1 - iProb1)/(1 - iProb1 - iProb2 + iProb12))) *(0.5*(nuij - 1 - cij + 2*iOR*iProb2)/nuij - 1);
        iGrad[J + 1 + j] = ((iProb2*(1 - iProb2)/(1 - iProb1 - iProb2 + iProb12))) *(0.5*(nuij - 1 - cij + 2*iOR*iProb1)/nuij - 1);
        iGrad[2*J + 1 + j] =  iOR*
              (((iOR-1)*((iProb1 + iProb2) - ((1+cij)*(iProb1 + iProb2) - 2*iProb1*iProb2*(2*iOR-1))/nuij)) - 1 - cij + nuij ) /
                    (2*(iOR - 1)*(iOR - 1)*(1 - iProb1 - iProb2 + iProb12));
      for (int k =0; k < pNT; ++k)
      {
      iGrad[3*J + k + 1] = ((XNT(i,k)*iProb1*(1 - iProb1)/(1 - iProb1 - iProb2 + iProb12))) *(0.5*(nuij - 1 - cij + 2*iOR*iProb2)/nuij - 1);
      }
      for (int k =0; k < pT; ++k)
      {
      iGrad[3*J + pNT + k + 1] = ((XT(i,k)*iProb2*(1 - iProb2)/(1 - iProb1 - iProb2 + iProb12))) *(0.5*(nuij - 1 - cij + 2*iOR*iProb1)/nuij - 1);
      }
      for (int k =0; k < pOR; ++k)
      {
      iGrad[3*J + pNT + pT + k + 1] =  XOR(i,k)*iOR*
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
  return(-Grad);
     }
