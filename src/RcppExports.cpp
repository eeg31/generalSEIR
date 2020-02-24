// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// stepModel
IntegerVector stepModel(IntegerVector state, NumericVector r, IntegerVector spec, int nevents);
RcppExport SEXP _generalSEIR_stepModel(SEXP stateSEXP, SEXP rSEXP, SEXP specSEXP, SEXP neventsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type state(stateSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type spec(specSEXP);
    Rcpp::traits::input_parameter< int >::type nevents(neventsSEXP);
    rcpp_result_gen = Rcpp::wrap(stepModel(state, r, spec, nevents));
    return rcpp_result_gen;
END_RCPP
}
// tauleap
IntegerVector tauleap(IntegerVector state, NumericVector r, IntegerVector spec, int nevents, double tau);
RcppExport SEXP _generalSEIR_tauleap(SEXP stateSEXP, SEXP rSEXP, SEXP specSEXP, SEXP neventsSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type state(stateSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type spec(specSEXP);
    Rcpp::traits::input_parameter< int >::type nevents(neventsSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(tauleap(state, r, spec, nevents, tau));
    return rcpp_result_gen;
END_RCPP
}
// getRates
NumericVector getRates(IntegerVector var, NumericVector par, int nevents);
RcppExport SEXP _generalSEIR_getRates(SEXP varSEXP, SEXP parSEXP, SEXP neventsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type var(varSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< int >::type nevents(neventsSEXP);
    rcpp_result_gen = Rcpp::wrap(getRates(var, par, nevents));
    return rcpp_result_gen;
END_RCPP
}
// simulate1
NumericVector simulate1(IntegerVector state, NumericVector par, IntegerVector spec, int nevents, int tmax, int inc);
RcppExport SEXP _generalSEIR_simulate1(SEXP stateSEXP, SEXP parSEXP, SEXP specSEXP, SEXP neventsSEXP, SEXP tmaxSEXP, SEXP incSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type state(stateSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type spec(specSEXP);
    Rcpp::traits::input_parameter< int >::type nevents(neventsSEXP);
    Rcpp::traits::input_parameter< int >::type tmax(tmaxSEXP);
    Rcpp::traits::input_parameter< int >::type inc(incSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate1(state, par, spec, nevents, tmax, inc));
    return rcpp_result_gen;
END_RCPP
}
// simFixed
IntegerVector simFixed(IntegerVector var, NumericVector par, IntegerVector spec, int nevents, int tmax);
RcppExport SEXP _generalSEIR_simFixed(SEXP varSEXP, SEXP parSEXP, SEXP specSEXP, SEXP neventsSEXP, SEXP tmaxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type var(varSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type spec(specSEXP);
    Rcpp::traits::input_parameter< int >::type nevents(neventsSEXP);
    Rcpp::traits::input_parameter< int >::type tmax(tmaxSEXP);
    rcpp_result_gen = Rcpp::wrap(simFixed(var, par, spec, nevents, tmax));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_generalSEIR_stepModel", (DL_FUNC) &_generalSEIR_stepModel, 4},
    {"_generalSEIR_tauleap", (DL_FUNC) &_generalSEIR_tauleap, 5},
    {"_generalSEIR_getRates", (DL_FUNC) &_generalSEIR_getRates, 3},
    {"_generalSEIR_simulate1", (DL_FUNC) &_generalSEIR_simulate1, 6},
    {"_generalSEIR_simFixed", (DL_FUNC) &_generalSEIR_simFixed, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_generalSEIR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}