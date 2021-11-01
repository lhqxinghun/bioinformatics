// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// rundiffTADs
Rcpp::List rundiffTADs(Rcpp::List& matlistcon1, Rcpp::List& matlistcon2, Rcpp::String namecon1, Rcpp::String namecon2, Rcpp::List parnorm, Rcpp::List parhierTADs, Rcpp::List parFTRLreg);
RcppExport SEXP _TADfit_rundiffTADs(SEXP matlistcon1SEXP, SEXP matlistcon2SEXP, SEXP namecon1SEXP, SEXP namecon2SEXP, SEXP parnormSEXP, SEXP parhierTADsSEXP, SEXP parFTRLregSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type matlistcon1(matlistcon1SEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type matlistcon2(matlistcon2SEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type namecon1(namecon1SEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type namecon2(namecon2SEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type parnorm(parnormSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type parhierTADs(parhierTADsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type parFTRLreg(parFTRLregSEXP);
    rcpp_result_gen = Rcpp::wrap(rundiffTADs(matlistcon1, matlistcon2, namecon1, namecon2, parnorm, parhierTADs, parFTRLreg));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TADfit_rundiffTADs", (DL_FUNC) &_TADfit_rundiffTADs, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_TADfit(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
