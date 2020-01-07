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
arma::vec gradcpp(SEXP eta, SEXP beta, SEXP doc_ct, SEXP mu, SEXP siginv);
RcppExport SEXP _stm_gradcpp(SEXP etaSEXP, SEXP betaSEXP, SEXP doc_ctSEXP, SEXP muSEXP, SEXP siginvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type doc_ct(doc_ctSEXP);
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
// sumcpp
arma::mat sumcpp(SEXP A_, SEXP B_);
RcppExport SEXP _stm_sumcpp(SEXP A_SEXP, SEXP B_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type A_(A_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type B_(B_SEXP);
    rcpp_result_gen = Rcpp::wrap(sumcpp(A_, B_));
    return rcpp_result_gen;
END_RCPP
}
// asumcpp
void asumcpp(SEXP SUM_, SEXP A_, SEXP B_);
RcppExport SEXP _stm_asumcpp(SEXP SUM_SEXP, SEXP A_SEXP, SEXP B_SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type SUM_(SUM_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type A_(A_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type B_(B_SEXP);
    asumcpp(SUM_, A_, B_);
    return R_NilValue;
END_RCPP
}
// estepcpp
SEXP estepcpp(SEXP docs_, SEXP beta_idx_, SEXP update_mu_, SEXP beta_, SEXP lambda_old_, SEXP mu_, SEXP sigma_, SEXP method_, SEXP verbose_);
RcppExport SEXP _stm_estepcpp(SEXP docs_SEXP, SEXP beta_idx_SEXP, SEXP update_mu_SEXP, SEXP beta_SEXP, SEXP lambda_old_SEXP, SEXP mu_SEXP, SEXP sigma_SEXP, SEXP method_SEXP, SEXP verbose_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type docs_(docs_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type beta_idx_(beta_idx_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type update_mu_(update_mu_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type lambda_old_(lambda_old_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type mu_(mu_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type sigma_(sigma_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type method_(method_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type verbose_(verbose_SEXP);
    rcpp_result_gen = Rcpp::wrap(estepcpp(docs_, beta_idx_, update_mu_, beta_, lambda_old_, mu_, sigma_, method_, verbose_));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_stm_lhoodcpp", (DL_FUNC) &_stm_lhoodcpp, 5},
    {"_stm_gradcpp", (DL_FUNC) &_stm_gradcpp, 5},
    {"_stm_hpbcpp", (DL_FUNC) &_stm_hpbcpp, 6},
    {"_stm_sumcpp", (DL_FUNC) &_stm_sumcpp, 2},
    {"_stm_asumcpp", (DL_FUNC) &_stm_asumcpp, 3},
    {"_stm_estepcpp", (DL_FUNC) &_stm_estepcpp, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_stm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
