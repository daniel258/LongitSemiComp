// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// FirstOne
int FirstOne(NumericVector x);
RcppExport SEXP _LongitSemiComp_FirstOne(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(FirstOne(x));
    return rcpp_result_gen;
END_RCPP
}
// GradPenalLogLik
arma::vec GradPenalLogLik(arma::vec param, arma::mat X, arma::mat YNT, arma::mat riskNT, arma::mat riskT, arma::mat YT, arma::mat TimeBase, arma::mat TimePen, arma::vec lambda);
RcppExport SEXP _LongitSemiComp_GradPenalLogLik(SEXP paramSEXP, SEXP XSEXP, SEXP YNTSEXP, SEXP riskNTSEXP, SEXP riskTSEXP, SEXP YTSEXP, SEXP TimeBaseSEXP, SEXP TimePenSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type param(paramSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YNT(YNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskNT(riskNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskT(riskTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YT(YTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimeBase(TimeBaseSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimePen(TimePenSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(GradPenalLogLik(param, X, YNT, riskNT, riskT, YT, TimeBase, TimePen, lambda));
    return rcpp_result_gen;
END_RCPP
}
// GradPenalLogLikPers
arma::mat GradPenalLogLikPers(arma::vec param, arma::mat X, arma::mat YNT, arma::mat riskNT, arma::mat riskT, arma::mat YT, arma::mat TimeBase, arma::mat TimePen, arma::vec lambda);
RcppExport SEXP _LongitSemiComp_GradPenalLogLikPers(SEXP paramSEXP, SEXP XSEXP, SEXP YNTSEXP, SEXP riskNTSEXP, SEXP riskTSEXP, SEXP YTSEXP, SEXP TimeBaseSEXP, SEXP TimePenSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type param(paramSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YNT(YNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskNT(riskNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskT(riskTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YT(YTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimeBase(TimeBaseSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimePen(TimePenSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(GradPenalLogLikPers(param, X, YNT, riskNT, riskT, YT, TimeBase, TimePen, lambda));
    return rcpp_result_gen;
END_RCPP
}
// OddsRatioTime
List OddsRatioTime(NumericMatrix YNT, NumericMatrix YT, NumericMatrix riskNT, NumericMatrix riskT);
RcppExport SEXP _LongitSemiComp_OddsRatioTime(SEXP YNTSEXP, SEXP YTSEXP, SEXP riskNTSEXP, SEXP riskTSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type YNT(YNTSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type YT(YTSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type riskNT(riskNTSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type riskT(riskTSEXP);
    rcpp_result_gen = Rcpp::wrap(OddsRatioTime(YNT, YT, riskNT, riskT));
    return rcpp_result_gen;
END_RCPP
}
// ParamLogLik
double ParamLogLik(arma::vec param, arma::mat X, arma::mat YNT, arma::mat YT, arma::mat riskNT, arma::mat riskT);
RcppExport SEXP _LongitSemiComp_ParamLogLik(SEXP paramSEXP, SEXP XSEXP, SEXP YNTSEXP, SEXP YTSEXP, SEXP riskNTSEXP, SEXP riskTSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type param(paramSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YNT(YNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YT(YTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskNT(riskNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskT(riskTSEXP);
    rcpp_result_gen = Rcpp::wrap(ParamLogLik(param, X, YNT, YT, riskNT, riskT));
    return rcpp_result_gen;
END_RCPP
}
// PenalLogLik
double PenalLogLik(arma::vec param, arma::mat X, arma::mat YNT, arma::mat riskNT, arma::mat riskT, arma::mat YT, arma::mat TimeBase, arma::mat TimePen, arma::vec lambda);
RcppExport SEXP _LongitSemiComp_PenalLogLik(SEXP paramSEXP, SEXP XSEXP, SEXP YNTSEXP, SEXP riskNTSEXP, SEXP riskTSEXP, SEXP YTSEXP, SEXP TimeBaseSEXP, SEXP TimePenSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type param(paramSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YNT(YNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskNT(riskNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskT(riskTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YT(YTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimeBase(TimeBaseSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimePen(TimePenSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(PenalLogLik(param, X, YNT, riskNT, riskT, YT, TimeBase, TimePen, lambda));
    return rcpp_result_gen;
END_RCPP
}
// PenalLogLikNoX
double PenalLogLikNoX(arma::vec param, arma::mat YNT, arma::mat riskNT, arma::mat riskT, arma::mat YT, arma::mat TimeBase, arma::mat TimePen, arma::vec lambda);
RcppExport SEXP _LongitSemiComp_PenalLogLikNoX(SEXP paramSEXP, SEXP YNTSEXP, SEXP riskNTSEXP, SEXP riskTSEXP, SEXP YTSEXP, SEXP TimeBaseSEXP, SEXP TimePenSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type param(paramSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YNT(YNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskNT(riskNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskT(riskTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YT(YTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimeBase(TimeBaseSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimePen(TimePenSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(PenalLogLikNoX(param, YNT, riskNT, riskT, YT, TimeBase, TimePen, lambda));
    return rcpp_result_gen;
END_RCPP
}
// PenalLogLikORconst
double PenalLogLikORconst(arma::vec param, arma::mat X, arma::mat YNT, arma::mat riskNT, arma::mat riskT, arma::mat YT, arma::mat TimeBase, arma::mat TimePen, arma::vec lambda);
RcppExport SEXP _LongitSemiComp_PenalLogLikORconst(SEXP paramSEXP, SEXP XSEXP, SEXP YNTSEXP, SEXP riskNTSEXP, SEXP riskTSEXP, SEXP YTSEXP, SEXP TimeBaseSEXP, SEXP TimePenSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type param(paramSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YNT(YNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskNT(riskNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskT(riskTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YT(YTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimeBase(TimeBaseSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimePen(TimePenSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(PenalLogLikORconst(param, X, YNT, riskNT, riskT, YT, TimeBase, TimePen, lambda));
    return rcpp_result_gen;
END_RCPP
}
// PenalLogLikORconstNOBetaOR
double PenalLogLikORconstNOBetaOR(arma::vec param, arma::mat X, arma::mat YNT, arma::mat riskNT, arma::mat riskT, arma::mat YT, arma::mat TimeBase, arma::mat TimePen, arma::vec lambda);
RcppExport SEXP _LongitSemiComp_PenalLogLikORconstNOBetaOR(SEXP paramSEXP, SEXP XSEXP, SEXP YNTSEXP, SEXP riskNTSEXP, SEXP riskTSEXP, SEXP YTSEXP, SEXP TimeBaseSEXP, SEXP TimePenSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type param(paramSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YNT(YNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskNT(riskNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskT(riskTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YT(YTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimeBase(TimeBaseSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimePen(TimePenSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(PenalLogLikORconstNOBetaOR(param, X, YNT, riskNT, riskT, YT, TimeBase, TimePen, lambda));
    return rcpp_result_gen;
END_RCPP
}
// PenalLogLikORconstNoX
double PenalLogLikORconstNoX(arma::vec param, arma::mat YNT, arma::mat riskNT, arma::mat riskT, arma::mat YT, arma::mat TimeBase, arma::mat TimePen, arma::vec lambda);
RcppExport SEXP _LongitSemiComp_PenalLogLikORconstNoX(SEXP paramSEXP, SEXP YNTSEXP, SEXP riskNTSEXP, SEXP riskTSEXP, SEXP YTSEXP, SEXP TimeBaseSEXP, SEXP TimePenSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type param(paramSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YNT(YNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskNT(riskNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskT(riskTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YT(YTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimeBase(TimeBaseSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimePen(TimePenSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(PenalLogLikORconstNoX(param, YNT, riskNT, riskT, YT, TimeBase, TimePen, lambda));
    return rcpp_result_gen;
END_RCPP
}
// PenalLogLikORparam
double PenalLogLikORparam(arma::vec param, arma::mat X, arma::mat YNT, arma::mat riskNT, arma::mat riskT, arma::mat YT, arma::mat TimeBase, arma::mat TimePen, arma::vec lambda, arma::vec PieceWiseTimes);
RcppExport SEXP _LongitSemiComp_PenalLogLikORparam(SEXP paramSEXP, SEXP XSEXP, SEXP YNTSEXP, SEXP riskNTSEXP, SEXP riskTSEXP, SEXP YTSEXP, SEXP TimeBaseSEXP, SEXP TimePenSEXP, SEXP lambdaSEXP, SEXP PieceWiseTimesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type param(paramSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YNT(YNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskNT(riskNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskT(riskTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YT(YTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimeBase(TimeBaseSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimePen(TimePenSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type PieceWiseTimes(PieceWiseTimesSEXP);
    rcpp_result_gen = Rcpp::wrap(PenalLogLikORparam(param, X, YNT, riskNT, riskT, YT, TimeBase, TimePen, lambda, PieceWiseTimes));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_LongitSemiComp_FirstOne", (DL_FUNC) &_LongitSemiComp_FirstOne, 1},
    {"_LongitSemiComp_GradPenalLogLik", (DL_FUNC) &_LongitSemiComp_GradPenalLogLik, 9},
    {"_LongitSemiComp_GradPenalLogLikPers", (DL_FUNC) &_LongitSemiComp_GradPenalLogLikPers, 9},
    {"_LongitSemiComp_OddsRatioTime", (DL_FUNC) &_LongitSemiComp_OddsRatioTime, 4},
    {"_LongitSemiComp_ParamLogLik", (DL_FUNC) &_LongitSemiComp_ParamLogLik, 6},
    {"_LongitSemiComp_PenalLogLik", (DL_FUNC) &_LongitSemiComp_PenalLogLik, 9},
    {"_LongitSemiComp_PenalLogLikNoX", (DL_FUNC) &_LongitSemiComp_PenalLogLikNoX, 8},
    {"_LongitSemiComp_PenalLogLikORconst", (DL_FUNC) &_LongitSemiComp_PenalLogLikORconst, 9},
    {"_LongitSemiComp_PenalLogLikORconstNOBetaOR", (DL_FUNC) &_LongitSemiComp_PenalLogLikORconstNOBetaOR, 9},
    {"_LongitSemiComp_PenalLogLikORconstNoX", (DL_FUNC) &_LongitSemiComp_PenalLogLikORconstNoX, 8},
    {"_LongitSemiComp_PenalLogLikORparam", (DL_FUNC) &_LongitSemiComp_PenalLogLikORparam, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_LongitSemiComp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
