// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// lhoodcpp
double lhoodcpp(SEXP eta, SEXP beta, SEXP doc_ct, SEXP mu, SEXP siginv);
RcppExport SEXP _stm_lhoodcpp(SEXP etaSEXP, SEXP betaSEXP, SEXP doc_ctSEXP, SEXP muSEXP, SEXP siginvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type doc_ct(doc_ctSEXP);
    Rcpp::traits::input_parameter< SEXP >::type mu(muSEXP);
    Rcpp::traits::input_parameter< SEXP >::type siginv(siginvSEXP);
    rcpp_result_gen = Rcpp::wrap(lhoodcpp(eta, beta, doc_ct, mu, siginv));
    return rcpp_result_gen;
END_RCPP
}
// gradcpp
arma::vec gradcpp(SEXP eta, SEXP beta, const arma::uvec& doc_ct, SEXP mu, SEXP siginv);
RcppExport SEXP _stm_gradcpp(SEXP etaSEXP, SEXP betaSEXP, SEXP doc_ctSEXP, SEXP muSEXP, SEXP siginvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type doc_ct(doc_ctSEXP);
    Rcpp::traits::input_parameter< SEXP >::type mu(muSEXP);
    Rcpp::traits::input_parameter< SEXP >::type siginv(siginvSEXP);
    rcpp_result_gen = Rcpp::wrap(gradcpp(eta, beta, doc_ct, mu, siginv));
    return rcpp_result_gen;
END_RCPP
}
// hpbcpp
SEXP hpbcpp(SEXP eta, SEXP beta, SEXP doc_ct, SEXP mu, SEXP siginv, SEXP sigmaentropy);
RcppExport SEXP _stm_hpbcpp(SEXP etaSEXP, SEXP betaSEXP, SEXP doc_ctSEXP, SEXP muSEXP, SEXP siginvSEXP, SEXP sigmaentropySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type doc_ct(doc_ctSEXP);
    Rcpp::traits::input_parameter< SEXP >::type mu(muSEXP);
    Rcpp::traits::input_parameter< SEXP >::type siginv(siginvSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sigmaentropy(sigmaentropySEXP);
    rcpp_result_gen = Rcpp::wrap(hpbcpp(eta, beta, doc_ct, mu, siginv, sigmaentropy));
    return rcpp_result_gen;
END_RCPP
}
// n_beta_sumcpp
void n_beta_sumcpp(SEXP sumc_, const arma::uvec& aw, SEXP input_);
RcppExport SEXP _stm_n_beta_sumcpp(SEXP sumc_SEXP, SEXP awSEXP, SEXP input_SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type sumc_(sumc_SEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type aw(awSEXP);
    Rcpp::traits::input_parameter< SEXP >::type input_(input_SEXP);
    n_beta_sumcpp(sumc_, aw, input_);
    return R_NilValue;
END_RCPP
}
// n_sigma_sumcpp
void n_sigma_sumcpp(SEXP sumc_, SEXP input_);
RcppExport SEXP _stm_n_sigma_sumcpp(SEXP sumc_SEXP, SEXP input_SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type sumc_(sumc_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type input_(input_SEXP);
    n_sigma_sumcpp(sumc_, input_);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_stm_lhoodcpp", (DL_FUNC) &_stm_lhoodcpp, 5},
    {"_stm_gradcpp", (DL_FUNC) &_stm_gradcpp, 5},
    {"_stm_hpbcpp", (DL_FUNC) &_stm_hpbcpp, 6},
    {"_stm_n_beta_sumcpp", (DL_FUNC) &_stm_n_beta_sumcpp, 3},
    {"_stm_n_sigma_sumcpp", (DL_FUNC) &_stm_n_sigma_sumcpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_stm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
