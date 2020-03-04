// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// GradParamLogLik
arma::vec GradParamLogLik(arma::vec param, arma::mat YT, arma::mat YNT, arma::mat riskNT, arma::mat riskT, arma::mat XNT, arma::mat XT, arma::mat XOR, double epsOR);
RcppExport SEXP _LongitSemiComp_GradParamLogLik(SEXP paramSEXP, SEXP YTSEXP, SEXP YNTSEXP, SEXP riskNTSEXP, SEXP riskTSEXP, SEXP XNTSEXP, SEXP XTSEXP, SEXP XORSEXP, SEXP epsORSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type param(paramSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YT(YTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YNT(YNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskNT(riskNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskT(riskTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XNT(XNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XT(XTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XOR(XORSEXP);
    Rcpp::traits::input_parameter< double >::type epsOR(epsORSEXP);
    rcpp_result_gen = Rcpp::wrap(GradParamLogLik(param, YT, YNT, riskNT, riskT, XNT, XT, XOR, epsOR));
    return rcpp_result_gen;
END_RCPP
}
// GradParamLogLikPersTimeDep
arma::mat GradParamLogLikPersTimeDep(arma::vec param, arma::vec ID, arma::uvec TM, arma::vec YT, arma::vec YNT, arma::mat XNT, arma::mat XT, arma::mat XOR, double epsOR);
RcppExport SEXP _LongitSemiComp_GradParamLogLikPersTimeDep(SEXP paramSEXP, SEXP IDSEXP, SEXP TMSEXP, SEXP YTSEXP, SEXP YNTSEXP, SEXP XNTSEXP, SEXP XTSEXP, SEXP XORSEXP, SEXP epsORSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type param(paramSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ID(IDSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type TM(TMSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type YT(YTSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type YNT(YNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XNT(XNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XT(XTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XOR(XORSEXP);
    Rcpp::traits::input_parameter< double >::type epsOR(epsORSEXP);
    rcpp_result_gen = Rcpp::wrap(GradParamLogLikPersTimeDep(param, ID, TM, YT, YNT, XNT, XT, XOR, epsOR));
    return rcpp_result_gen;
END_RCPP
}
// GradParamLogLikTimeDep
arma::rowvec GradParamLogLikTimeDep(arma::vec param, arma::vec ID, arma::uvec TM, arma::vec YT, arma::vec YNT, arma::mat XNT, arma::mat XT, arma::mat XOR, double epsOR);
RcppExport SEXP _LongitSemiComp_GradParamLogLikTimeDep(SEXP paramSEXP, SEXP IDSEXP, SEXP TMSEXP, SEXP YTSEXP, SEXP YNTSEXP, SEXP XNTSEXP, SEXP XTSEXP, SEXP XORSEXP, SEXP epsORSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type param(paramSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ID(IDSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type TM(TMSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type YT(YTSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type YNT(YNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XNT(XNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XT(XTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XOR(XORSEXP);
    Rcpp::traits::input_parameter< double >::type epsOR(epsORSEXP);
    rcpp_result_gen = Rcpp::wrap(GradParamLogLikTimeDep(param, ID, TM, YT, YNT, XNT, XT, XOR, epsOR));
    return rcpp_result_gen;
END_RCPP
}
// GradPenalLogLik
arma::vec GradPenalLogLik(arma::vec param, arma::mat YNT, arma::mat riskNT, arma::mat riskT, arma::mat YT, arma::mat XNT, arma::mat XT, arma::mat XOR, arma::mat TimeBase, arma::mat TimePen, arma::vec lambda, double epsOR);
RcppExport SEXP _LongitSemiComp_GradPenalLogLik(SEXP paramSEXP, SEXP YNTSEXP, SEXP riskNTSEXP, SEXP riskTSEXP, SEXP YTSEXP, SEXP XNTSEXP, SEXP XTSEXP, SEXP XORSEXP, SEXP TimeBaseSEXP, SEXP TimePenSEXP, SEXP lambdaSEXP, SEXP epsORSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type param(paramSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YNT(YNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskNT(riskNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskT(riskTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YT(YTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XNT(XNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XT(XTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XOR(XORSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimeBase(TimeBaseSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimePen(TimePenSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type epsOR(epsORSEXP);
    rcpp_result_gen = Rcpp::wrap(GradPenalLogLik(param, YNT, riskNT, riskT, YT, XNT, XT, XOR, TimeBase, TimePen, lambda, epsOR));
    return rcpp_result_gen;
END_RCPP
}
// GradPenalLogLikNullModelOR
arma::vec GradPenalLogLikNullModelOR(arma::vec param, arma::mat YNT, arma::mat riskNT, arma::mat riskT, arma::mat YT, arma::mat XNT, arma::mat XT, arma::mat TimeBase, arma::mat TimePen, arma::vec lambda, double epsOR);
RcppExport SEXP _LongitSemiComp_GradPenalLogLikNullModelOR(SEXP paramSEXP, SEXP YNTSEXP, SEXP riskNTSEXP, SEXP riskTSEXP, SEXP YTSEXP, SEXP XNTSEXP, SEXP XTSEXP, SEXP TimeBaseSEXP, SEXP TimePenSEXP, SEXP lambdaSEXP, SEXP epsORSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type param(paramSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YNT(YNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskNT(riskNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskT(riskTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YT(YTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XNT(XNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XT(XTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimeBase(TimeBaseSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimePen(TimePenSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type epsOR(epsORSEXP);
    rcpp_result_gen = Rcpp::wrap(GradPenalLogLikNullModelOR(param, YNT, riskNT, riskT, YT, XNT, XT, TimeBase, TimePen, lambda, epsOR));
    return rcpp_result_gen;
END_RCPP
}
// GradPenalLogLikPers
arma::mat GradPenalLogLikPers(arma::vec param, arma::mat YNT, arma::mat riskNT, arma::mat riskT, arma::mat YT, arma::mat XNT, arma::mat XT, arma::mat XOR, arma::mat TimeBase, arma::mat TimePen, arma::vec lambda, double epsOR);
RcppExport SEXP _LongitSemiComp_GradPenalLogLikPers(SEXP paramSEXP, SEXP YNTSEXP, SEXP riskNTSEXP, SEXP riskTSEXP, SEXP YTSEXP, SEXP XNTSEXP, SEXP XTSEXP, SEXP XORSEXP, SEXP TimeBaseSEXP, SEXP TimePenSEXP, SEXP lambdaSEXP, SEXP epsORSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type param(paramSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YNT(YNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskNT(riskNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskT(riskTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YT(YTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XNT(XNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XT(XTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XOR(XORSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimeBase(TimeBaseSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimePen(TimePenSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type epsOR(epsORSEXP);
    rcpp_result_gen = Rcpp::wrap(GradPenalLogLikPers(param, YNT, riskNT, riskT, YT, XNT, XT, XOR, TimeBase, TimePen, lambda, epsOR));
    return rcpp_result_gen;
END_RCPP
}
// GradPenalLogLikPersNullModelOR
arma::mat GradPenalLogLikPersNullModelOR(arma::vec param, arma::mat YNT, arma::mat riskNT, arma::mat riskT, arma::mat YT, arma::mat XNT, arma::mat XT, arma::mat TimeBase, arma::mat TimePen, arma::vec lambda, double epsOR);
RcppExport SEXP _LongitSemiComp_GradPenalLogLikPersNullModelOR(SEXP paramSEXP, SEXP YNTSEXP, SEXP riskNTSEXP, SEXP riskTSEXP, SEXP YTSEXP, SEXP XNTSEXP, SEXP XTSEXP, SEXP TimeBaseSEXP, SEXP TimePenSEXP, SEXP lambdaSEXP, SEXP epsORSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type param(paramSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YNT(YNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskNT(riskNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskT(riskTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YT(YTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XNT(XNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XT(XTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimeBase(TimeBaseSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimePen(TimePenSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type epsOR(epsORSEXP);
    rcpp_result_gen = Rcpp::wrap(GradPenalLogLikPersNullModelOR(param, YNT, riskNT, riskT, YT, XNT, XT, TimeBase, TimePen, lambda, epsOR));
    return rcpp_result_gen;
END_RCPP
}
// GradPenalLogLikPersTimeDep
arma::mat GradPenalLogLikPersTimeDep(arma::vec param, arma::vec ID, arma::uvec TM, arma::vec YT, arma::vec YNT, arma::mat XNT, arma::mat XT, arma::mat XOR, arma::mat TimeBase, arma::mat TimePen, arma::vec lambda, double epsOR);
RcppExport SEXP _LongitSemiComp_GradPenalLogLikPersTimeDep(SEXP paramSEXP, SEXP IDSEXP, SEXP TMSEXP, SEXP YTSEXP, SEXP YNTSEXP, SEXP XNTSEXP, SEXP XTSEXP, SEXP XORSEXP, SEXP TimeBaseSEXP, SEXP TimePenSEXP, SEXP lambdaSEXP, SEXP epsORSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type param(paramSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ID(IDSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type TM(TMSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type YT(YTSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type YNT(YNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XNT(XNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XT(XTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XOR(XORSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimeBase(TimeBaseSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimePen(TimePenSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type epsOR(epsORSEXP);
    rcpp_result_gen = Rcpp::wrap(GradPenalLogLikPersTimeDep(param, ID, TM, YT, YNT, XNT, XT, XOR, TimeBase, TimePen, lambda, epsOR));
    return rcpp_result_gen;
END_RCPP
}
// GradPenalLogLikTimeDep
arma::rowvec GradPenalLogLikTimeDep(arma::vec param, arma::vec ID, arma::uvec TM, arma::vec YT, arma::vec YNT, arma::mat XNT, arma::mat XT, arma::mat XOR, arma::mat TimeBase, arma::mat TimePen, arma::vec lambda, double epsOR);
RcppExport SEXP _LongitSemiComp_GradPenalLogLikTimeDep(SEXP paramSEXP, SEXP IDSEXP, SEXP TMSEXP, SEXP YTSEXP, SEXP YNTSEXP, SEXP XNTSEXP, SEXP XTSEXP, SEXP XORSEXP, SEXP TimeBaseSEXP, SEXP TimePenSEXP, SEXP lambdaSEXP, SEXP epsORSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type param(paramSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ID(IDSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type TM(TMSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type YT(YTSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type YNT(YNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XNT(XNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XT(XTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XOR(XORSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimeBase(TimeBaseSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimePen(TimePenSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type epsOR(epsORSEXP);
    rcpp_result_gen = Rcpp::wrap(GradPenalLogLikTimeDep(param, ID, TM, YT, YNT, XNT, XT, XOR, TimeBase, TimePen, lambda, epsOR));
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
double ParamLogLik(arma::vec param, arma::mat YT, arma::mat YNT, arma::mat riskNT, arma::mat riskT, arma::mat XNT, arma::mat XT, arma::mat XOR, double epsOR);
RcppExport SEXP _LongitSemiComp_ParamLogLik(SEXP paramSEXP, SEXP YTSEXP, SEXP YNTSEXP, SEXP riskNTSEXP, SEXP riskTSEXP, SEXP XNTSEXP, SEXP XTSEXP, SEXP XORSEXP, SEXP epsORSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type param(paramSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YT(YTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YNT(YNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskNT(riskNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskT(riskTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XNT(XNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XT(XTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XOR(XORSEXP);
    Rcpp::traits::input_parameter< double >::type epsOR(epsORSEXP);
    rcpp_result_gen = Rcpp::wrap(ParamLogLik(param, YT, YNT, riskNT, riskT, XNT, XT, XOR, epsOR));
    return rcpp_result_gen;
END_RCPP
}
// ParamLogLikTimeDep
double ParamLogLikTimeDep(arma::vec param, arma::vec ID, arma::uvec TM, arma::vec YT, arma::vec YNT, arma::mat XNT, arma::mat XT, arma::mat XOR, double epsOR);
RcppExport SEXP _LongitSemiComp_ParamLogLikTimeDep(SEXP paramSEXP, SEXP IDSEXP, SEXP TMSEXP, SEXP YTSEXP, SEXP YNTSEXP, SEXP XNTSEXP, SEXP XTSEXP, SEXP XORSEXP, SEXP epsORSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type param(paramSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ID(IDSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type TM(TMSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type YT(YTSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type YNT(YNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XNT(XNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XT(XTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XOR(XORSEXP);
    Rcpp::traits::input_parameter< double >::type epsOR(epsORSEXP);
    rcpp_result_gen = Rcpp::wrap(ParamLogLikTimeDep(param, ID, TM, YT, YNT, XNT, XT, XOR, epsOR));
    return rcpp_result_gen;
END_RCPP
}
// PenalLogLik
double PenalLogLik(arma::vec param, arma::mat YT, arma::mat YNT, arma::mat riskNT, arma::mat riskT, arma::mat XNT, arma::mat XT, arma::mat XOR, arma::mat TimeBase, arma::mat TimePen, arma::vec lambda, double epsOR);
RcppExport SEXP _LongitSemiComp_PenalLogLik(SEXP paramSEXP, SEXP YTSEXP, SEXP YNTSEXP, SEXP riskNTSEXP, SEXP riskTSEXP, SEXP XNTSEXP, SEXP XTSEXP, SEXP XORSEXP, SEXP TimeBaseSEXP, SEXP TimePenSEXP, SEXP lambdaSEXP, SEXP epsORSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type param(paramSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YT(YTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YNT(YNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskNT(riskNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskT(riskTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XNT(XNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XT(XTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XOR(XORSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimeBase(TimeBaseSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimePen(TimePenSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type epsOR(epsORSEXP);
    rcpp_result_gen = Rcpp::wrap(PenalLogLik(param, YT, YNT, riskNT, riskT, XNT, XT, XOR, TimeBase, TimePen, lambda, epsOR));
    return rcpp_result_gen;
END_RCPP
}
// PenalLogLikNoX
double PenalLogLikNoX(arma::vec param, arma::mat YNT, arma::mat riskNT, arma::mat riskT, arma::mat YT, arma::mat TimeBase, arma::mat TimePen, arma::vec lambda, double epsOR);
RcppExport SEXP _LongitSemiComp_PenalLogLikNoX(SEXP paramSEXP, SEXP YNTSEXP, SEXP riskNTSEXP, SEXP riskTSEXP, SEXP YTSEXP, SEXP TimeBaseSEXP, SEXP TimePenSEXP, SEXP lambdaSEXP, SEXP epsORSEXP) {
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
    Rcpp::traits::input_parameter< double >::type epsOR(epsORSEXP);
    rcpp_result_gen = Rcpp::wrap(PenalLogLikNoX(param, YNT, riskNT, riskT, YT, TimeBase, TimePen, lambda, epsOR));
    return rcpp_result_gen;
END_RCPP
}
// PenalLogLikNullModelOR
double PenalLogLikNullModelOR(arma::vec param, arma::mat YT, arma::mat YNT, arma::mat riskNT, arma::mat riskT, arma::mat XNT, arma::mat XT, arma::mat TimeBase, arma::mat TimePen, arma::vec lambda, double epsOR);
RcppExport SEXP _LongitSemiComp_PenalLogLikNullModelOR(SEXP paramSEXP, SEXP YTSEXP, SEXP YNTSEXP, SEXP riskNTSEXP, SEXP riskTSEXP, SEXP XNTSEXP, SEXP XTSEXP, SEXP TimeBaseSEXP, SEXP TimePenSEXP, SEXP lambdaSEXP, SEXP epsORSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type param(paramSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YT(YTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YNT(YNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskNT(riskNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type riskT(riskTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XNT(XNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XT(XTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimeBase(TimeBaseSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimePen(TimePenSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type epsOR(epsORSEXP);
    rcpp_result_gen = Rcpp::wrap(PenalLogLikNullModelOR(param, YT, YNT, riskNT, riskT, XNT, XT, TimeBase, TimePen, lambda, epsOR));
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
// PenalLogLikTimeDep
double PenalLogLikTimeDep(arma::vec param, arma::vec ID, arma::uvec TM, arma::vec YT, arma::vec YNT, arma::mat XNT, arma::mat XT, arma::mat XOR, arma::mat TimeBase, arma::mat TimePen, arma::vec lambda, double epsOR);
RcppExport SEXP _LongitSemiComp_PenalLogLikTimeDep(SEXP paramSEXP, SEXP IDSEXP, SEXP TMSEXP, SEXP YTSEXP, SEXP YNTSEXP, SEXP XNTSEXP, SEXP XTSEXP, SEXP XORSEXP, SEXP TimeBaseSEXP, SEXP TimePenSEXP, SEXP lambdaSEXP, SEXP epsORSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type param(paramSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ID(IDSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type TM(TMSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type YT(YTSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type YNT(YNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XNT(XNTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XT(XTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XOR(XORSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimeBase(TimeBaseSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TimePen(TimePenSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type epsOR(epsORSEXP);
    rcpp_result_gen = Rcpp::wrap(PenalLogLikTimeDep(param, ID, TM, YT, YNT, XNT, XT, XOR, TimeBase, TimePen, lambda, epsOR));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_LongitSemiComp_GradParamLogLik", (DL_FUNC) &_LongitSemiComp_GradParamLogLik, 9},
    {"_LongitSemiComp_GradParamLogLikPersTimeDep", (DL_FUNC) &_LongitSemiComp_GradParamLogLikPersTimeDep, 9},
    {"_LongitSemiComp_GradParamLogLikTimeDep", (DL_FUNC) &_LongitSemiComp_GradParamLogLikTimeDep, 9},
    {"_LongitSemiComp_GradPenalLogLik", (DL_FUNC) &_LongitSemiComp_GradPenalLogLik, 12},
    {"_LongitSemiComp_GradPenalLogLikNullModelOR", (DL_FUNC) &_LongitSemiComp_GradPenalLogLikNullModelOR, 11},
    {"_LongitSemiComp_GradPenalLogLikPers", (DL_FUNC) &_LongitSemiComp_GradPenalLogLikPers, 12},
    {"_LongitSemiComp_GradPenalLogLikPersNullModelOR", (DL_FUNC) &_LongitSemiComp_GradPenalLogLikPersNullModelOR, 11},
    {"_LongitSemiComp_GradPenalLogLikPersTimeDep", (DL_FUNC) &_LongitSemiComp_GradPenalLogLikPersTimeDep, 12},
    {"_LongitSemiComp_GradPenalLogLikTimeDep", (DL_FUNC) &_LongitSemiComp_GradPenalLogLikTimeDep, 12},
    {"_LongitSemiComp_OddsRatioTime", (DL_FUNC) &_LongitSemiComp_OddsRatioTime, 4},
    {"_LongitSemiComp_ParamLogLik", (DL_FUNC) &_LongitSemiComp_ParamLogLik, 9},
    {"_LongitSemiComp_ParamLogLikTimeDep", (DL_FUNC) &_LongitSemiComp_ParamLogLikTimeDep, 9},
    {"_LongitSemiComp_PenalLogLik", (DL_FUNC) &_LongitSemiComp_PenalLogLik, 12},
    {"_LongitSemiComp_PenalLogLikNoX", (DL_FUNC) &_LongitSemiComp_PenalLogLikNoX, 9},
    {"_LongitSemiComp_PenalLogLikNullModelOR", (DL_FUNC) &_LongitSemiComp_PenalLogLikNullModelOR, 11},
    {"_LongitSemiComp_PenalLogLikORconst", (DL_FUNC) &_LongitSemiComp_PenalLogLikORconst, 9},
    {"_LongitSemiComp_PenalLogLikORconstNoX", (DL_FUNC) &_LongitSemiComp_PenalLogLikORconstNoX, 8},
    {"_LongitSemiComp_PenalLogLikORparam", (DL_FUNC) &_LongitSemiComp_PenalLogLikORparam, 10},
    {"_LongitSemiComp_PenalLogLikTimeDep", (DL_FUNC) &_LongitSemiComp_PenalLogLikTimeDep, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_LongitSemiComp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
