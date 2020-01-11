// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]
#include "RcppArmadillo.h"
#include <R_ext/Applic.h>
#include "stm_types.h"

// [[Rcpp::export]]
double lhoodcpp(SEXP eta,
                   SEXP beta,
                   SEXP doc_ct,
                   SEXP mu,
                   SEXP siginv){
   
   Rcpp::NumericVector etav(eta); 
   arma::vec etas(etav.begin(), etav.size(), false);
   Rcpp::NumericMatrix betam(beta);
   arma::mat betas(betam.begin(), betam.nrow(), betam.ncol(), false);
   Rcpp::NumericVector doc_ctv(doc_ct);
   arma::vec doc_cts(doc_ctv.begin(), doc_ctv.size(), false);
   Rcpp::NumericVector muv(mu);
   arma::vec mus(muv.begin(), muv.size(), false);
   Rcpp::NumericMatrix siginvm(siginv);
   arma::mat siginvs(siginvm.begin(), siginvm.nrow(), siginvm.ncol(), false);
   
   arma::rowvec expeta(etas.size()+1); 
   expeta.fill(1);
   int neta = etav.size(); 
   for(int j=0; j <neta;  j++){
     expeta(j) = exp(etas(j));
   }
   double ndoc = sum(doc_cts);
   double part1 = arma::as_scalar(log(expeta*betas)*doc_cts - ndoc*log(sum(expeta)));
   arma::vec diff = etas - mus;
   double part2 = .5*arma::as_scalar(diff.t()*siginvs*diff);
   double out = part2 - part1;
   return out;
}

// [[Rcpp::export]]
arma::vec gradcpp(SEXP eta,
                   SEXP beta,
                   SEXP doc_ct,
                   SEXP mu,
                   SEXP siginv){
   
   Rcpp::NumericVector etav(eta); 
   arma::vec etas(etav.begin(), etav.size(), false);
   Rcpp::NumericMatrix betam(beta);
   arma::mat betas(betam.begin(), betam.nrow(), betam.ncol());
   Rcpp::NumericVector doc_ctv(doc_ct);
   arma::vec doc_cts(doc_ctv.begin(), doc_ctv.size(), false);
   Rcpp::NumericVector muv(mu);
   arma::vec mus(muv.begin(), muv.size(), false);
   Rcpp::NumericMatrix siginvm(siginv);
   arma::mat siginvs(siginvm.begin(), siginvm.nrow(), siginvm.ncol(), false);
   
    arma::colvec expeta(etas.size()+1); 
    expeta.fill(1);
    int neta = etas.size(); 
    for(int j=0; j <neta;  j++){
       expeta(j) = exp(etas(j));
    }
    betas.each_col() %= expeta;
    arma::vec part1 = betas*(doc_cts/arma::trans(sum(betas,0))) - (sum(doc_cts)/sum(expeta))*expeta;
    arma::vec part2 = siginvs*(etas - mus);
    part1.shed_row(neta);
    return part2-part1;
}

// [[Rcpp::export]]
SEXP hpbcpp(SEXP eta,
            SEXP beta,
            SEXP doc_ct,
            SEXP mu,
            SEXP siginv,
            SEXP sigmaentropy){
 
   Rcpp::NumericVector etav(eta); 
   arma::vec etas(etav.begin(), etav.size(), false);
   Rcpp::NumericMatrix betam(beta);
   arma::mat betas(betam.begin(), betam.nrow(), betam.ncol());
   Rcpp::NumericVector doc_ctv(doc_ct);
   arma::vec doc_cts(doc_ctv.begin(), doc_ctv.size(), false);
   Rcpp::NumericVector muv(mu);
   arma::vec mus(muv.begin(), muv.size(), false);
   Rcpp::NumericMatrix siginvm(siginv);
   arma::mat siginvs(siginvm.begin(), siginvm.nrow(), siginvm.ncol(), false);
   Rcpp::NumericVector sigmaentropym(sigmaentropy);
   arma::vec entropy(sigmaentropym);

   //Performance Nots from 3/6/2015
   //  I tried a few different variants and benchmarked this one as roughly twice as
   //  fast as the R code for a K=100 problem.  Key to performance was not creating
   //  too many objects and being selective in how things were flagged as triangular.
   //  Some additional notes in the code below.
   //
   //  Some things this doesn't have or I haven't tried
   //  - I didn't tweak the arguments much.  sigmaentropy is a double, and I'm still
   //    passing beta in the same way.  I tried doing a ", false" for beta but it didn't
   //    change much so I left it the same as in gradient.  
   //  - I tried treating the factors for doc_cts and colSums(EB) as a diagonal matrix- much slower.
   
   //  Haven't Tried/Done
   //  - each_row() might be much slower (not sure but arma is column order).  Maybe transpose in place?
   //  - depending on costs there are some really minor calculations that could be precomputed: 
   //     - sum(doc_ct)
   //     - sqrt(doc_ct)
   
   //  More on passing by reference here:
   //  - Hypothetically we could alter beta (because hessian is last thing we do) however down
   //    the road we may want to explore treating nonPD hessians by optimization at which point
   //    we would need it again.
   
   arma::colvec expeta(etas.size()+1); 
   expeta.fill(1);
   int neta = etas.size(); 
   for(int j=0; j <neta;  j++){
     expeta(j) = exp(etas(j));
   }
   arma::vec theta = expeta/sum(expeta);

   //create a new version of the matrix so we can mess with it
   arma::mat EB(betam.begin(), betam.nrow(), betam.ncol());
   //multiply each column by expeta
   EB.each_col() %= expeta; //this should be fastest as its column-major ordering
  
   //divide out by the column sums
   EB.each_row() %= arma::trans(sqrt(doc_cts))/sum(EB,0);
    
   //Combine the pieces of the Hessian which are matrices
   arma::mat hess = EB*EB.t() - sum(doc_cts)*(theta*theta.t());
  
   //we don't need EB any more so we turn it into phi
   EB.each_row() %= arma::trans(sqrt(doc_cts));
   
   //Now alter just the diagonal of the Hessian
   hess.diag() -= sum(EB,1) - sum(doc_cts)*theta;
   //Drop the last row and column
   hess.shed_row(neta);
   hess.shed_col(neta);
   //Now we can add in siginv
   hess = hess + siginvs;
   //At this point the Hessian is complete.
   
   //This next bit of code is from http://arma.sourceforge.net/docs.html#logging
   //It basically keeps arma from printing errors from chol to the console.
   std::ostream nullstream(0);
   //arma::set_stream_err2(nullstream);
   arma::set_cerr_stream(nullstream);
   
   ////
   //Invert via cholesky decomposition
   ////
   //Start by initializing an object
   arma::mat nu = arma::mat(hess.n_rows, hess.n_rows);
   //This version of chol generates a boolean which tells us if it failed.
   bool worked = arma::chol(nu,hess);
   if(!worked) {
     //It failed!  Oh Nos.
     // So the matrix wasn't positive definite.  In practice this means that it hasn't
     // converged probably along some minor aspect of the dimension.
     
     //Here we make it positive definite through diagonal dominance
     arma::vec dvec = hess.diag();
     //find the magnitude of the diagonal 
     arma::vec magnitudes = sum(abs(hess), 1) - abs(dvec);
     //iterate over each row and set the minimum value of the diagonal to be the magnitude of the other terms
     int Km1 = dvec.size();
     for(int j=0; j < Km1;  j++){
       if(arma::as_scalar(dvec(j)) < arma::as_scalar(magnitudes(j))) dvec(j) = magnitudes(j); //enforce diagonal dominance 
     }
     //overwrite the diagonal of the hessian with our new object
     hess.diag() = dvec;
     //that was sufficient to ensure positive definiteness so we now do cholesky
     nu = arma::chol(hess);
   }
   //compute 1/2 the determinant from the cholesky decomposition
   double detTerm = -sum(log(nu.diag()));
   
   //Now finish constructing nu
   nu = arma::inv(arma::trimatu(nu));
   nu = nu * nu.t(); //trimatu doesn't do anything for multiplication so it would just be timesink to signal here.
   
   //Precompute the difference since we use it twice
   arma::vec diff = etas - mus;
   //Now generate the bound and make it a scalar
   double bound = arma::as_scalar(log(arma::trans(theta)*betas)*doc_cts + detTerm - .5*diff.t()*siginvs*diff - entropy); 
   
   // Generate a return list that mimics the R output
   return Rcpp::List::create(
        Rcpp::Named("phis") = EB,
        Rcpp::Named("eta") = Rcpp::List::create(Rcpp::Named("lambda")=etas, Rcpp::Named("nu")=nu),
        Rcpp::Named("bound") = bound
        );
}


problem_data* problem_data_alloc(Rcpp::NumericMatrix* beta_,
                                 Rcpp::NumericVector* doc_ct_,
                                 Rcpp::NumericVector* mu_,
                                 Rcpp::NumericMatrix* siginv_) {
   problem_data* pr = new problem_data;
   pr->beta = new arma::mat(beta_->begin(), beta_->nrow(), beta_->ncol());
   pr->doc_ct = new arma::vec(doc_ct_->begin(), doc_ct_->size(), false);
   pr->mu = new arma::vec(mu_->begin(), mu_->size(), false);
   pr->siginv = new arma::mat(siginv_->begin(), siginv_->nrow(), siginv_->ncol(), false);
}

void* problem_data_free(problem_data* pr) {
   if(pr) {
      if(pr->beta)
         delete pr->beta;
      if(pr->doc_ct)
         delete pr->doc_ct;
      if(pr->mu)
         delete pr->mu;
      if(pr->siginv)
         delete pr->siginv;
   }
}


opt_method hash_opt_method (std::string const& method) {
   if (method == "BFGS") return BFGS;
   if (method == "L-BFGS-B") return LBFGSB;
   if (method == "CG") return CG;
}


double objlhood(int n, double *par, void *ex) {
   problem_data* pr = (problem_data*) ex;
   arma::vec eta(par, n, false, false);
   arma::rowvec expeta(n+1); 
   expeta.fill(1.0); 
   for(int j=0; j<n;  ++j)
      expeta(j) = exp(eta(j));
   
   double ndoc = sum(*pr->doc_ct);
   double first = arma::as_scalar(log(expeta*(*pr->beta))*(*pr->doc_ct) - ndoc*log(sum(expeta)));
   arma::vec diff = eta - *pr->mu;
   double second = .5*arma::as_scalar(diff.t()*(*pr->siginv)*diff);
   return (first - second);
}

/*
double objlhood(double* eta,
                long int K,
                long int V,
                double* beta,
                double* doc_ct,
                double* mu,
                double* siginv){
   
   std::cout<<"In objlhood..."<<std::endl;
   
   arma::vec etas(eta, K, false, false);
   std::cout<<"etas: "<<etas<<std::endl;
   
   arma::mat betas(beta, V, K, false, false);
   std::cout<<"betas: "<<betas<<std::endl;
   
   arma::vec doc_cts(doc_ct, V, false, false);
   std::cout<<"doc_cts: "<<doc_cts<<std::endl;
   
   arma::vec mus(mu, K, false, false);
   std::cout<<"mus: "<<mus<<std::endl;
   
   arma::mat siginvs(siginv, K, K, false, false);
   std::cout<<"siginvs: "<<siginvs<<std::endl;
   
   arma::rowvec expeta(etas.size()+1); 
   expeta.fill(1); 
   for(std::size_t j=0; j <K;  j++){
      expeta(j) = exp(etas(j));
   }
   std::cout<<"expeta: "<<expeta<<std::endl;
   
   arma::vec diff = etas - mus;
   std::cout<<"diff: "<<diff<<std::endl;
   
   double ndoc = sum(doc_cts);
   std::cout<<"ndoc: "<<ndoc<<std::endl;
   
   double part1 = arma::as_scalar(log(expeta*betas)*doc_cts - ndoc*log(sum(expeta)));
   double part2 = .5*arma::as_scalar(diff.t()*siginvs*diff);
   double out = part2 - part1;
   return out;
}
*/

void gradlhood(int n, double *par, double *gr, void *ex) {
   problem_data* pr = (problem_data*) ex;
   arma::vec agrad(n);
   arma::vec eta(par, n, false, false);
   arma::rowvec expeta(n+1); 
   expeta.fill(1); 
   for(int j=0; j<n;  j++) {
      expeta(j) = exp(eta(j));
   }
   
   pr->beta->each_col() %= expeta;
   arma::vec first = (*pr->beta)*((*pr->doc_ct)/arma::trans(sum((*pr->beta),0))) - (sum(*pr->doc_ct)/sum(expeta))*expeta;
   arma::vec second = *pr->siginv*(eta - *pr->mu);
   first.shed_row(n);
   agrad = second-first; 
   std::memcpy(gr, agrad.memptr(), sizeof(double)*n);
}

/*
void gradlhood(double* eta,
               double* grad,
               long int K,
               long int V,
               double* beta,
               double* doc_ct,
               double* mu,
               double* siginv){
   
   std::cout<<"In gradlhood..."<<std::endl;
   
   arma::vec etas(eta, K, false, false);
   std::cout<<"etas: "<<etas<<std::endl;
   
   arma::mat betas(beta, V, K, false, false);
   std::cout<<"betas: "<<betas<<std::endl;
   
   arma::vec doc_cts(doc_ct, V, false, false);
   std::cout<<"doc_cts: "<<doc_cts<<std::endl;
   
   arma::vec mus(mu, K, false, false);
   std::cout<<"mus: "<<mus<<std::endl;
   
   arma::mat siginvs(siginv, K, K, false, false);
   std::cout<<"siginvs: "<<siginvs<<std::endl;
   
   arma::vec agrad(K);
   
   arma::colvec expeta(etas.size()+1); 
   expeta.fill(1);
   for(int j=0; j<K;  ++j)
      expeta(j) = exp(etas(j));
   
   betas.each_col() %= expeta;
   arma::vec part1 = betas*(doc_cts/arma::trans(sum(betas,0))) - (sum(doc_cts)/sum(expeta))*expeta;
   arma::vec part2 = siginvs*(etas - mus);
   part1.shed_row(K);
   agrad = part2-part1;
   std::memcpy(grad, agrad.memptr(), sizeof(double)*K);
}
*/


// [[Rcpp::export]]
Rcpp::NumericVector sumcpp(SEXP A_, SEXP B_) {
   Rcpp::NumericVector A(A_);
   Rcpp::NumericVector B(B_);
   return (A + B);   
}

// [[Rcpp::export]]
void pluseqcpp(SEXP LHS_, SEXP RHS_) {
   Rcpp::NumericVector LHS(LHS_);
   Rcpp::NumericVector RHS(RHS_);
   LHS += RHS;   
}

// [[Rcpp::export]]
void pluseqcpp_idx(SEXP LHS_, SEXP RHS_, SEXP IDX_) {
   Rcpp::NumericMatrix LHS(LHS_);
   arma::mat ar_LHS(LHS.begin(), LHS.nrow(), LHS.ncol(), false);
   Rcpp::NumericMatrix RHS(RHS_);
   arma::mat ar_RHS(RHS.begin(), RHS.nrow(), RHS.ncol(), false);
   Rcpp::IntegerVector idx(IDX_);
   arma::irowvec ar_idx(idx);   
   for(int i=0; i<ar_idx.n_elem; ++i)
      ar_LHS.col((ar_idx[i]-1)) += ar_RHS.col(i);
}

// [[Rcpp::export]]
SEXP estepcpp(SEXP docs_, 
            SEXP beta_idx_, 
            SEXP update_mu_, // null allows for intercept only model  
            SEXP beta_, 
            SEXP lambda_old_, 
            SEXP mu_, 
            SEXP sigma_,  
            SEXP method_,
            SEXP verbose_) {
   
   ExactSum mySum;
   
   
   // According to this: https://thecoatlessprofessor.com/programming/cpp/unofficial-rcpp-api-documentation/#vmld
   // the code below should not do any deep copies of memory since the types are respected i.e. IntegerMatrix 
   // constructor actually gets an integer matrix passed to it etc.
   Rcpp::List docs(docs_);
   std::map<std::size_t, arma::imat> ar_docs;
   std::size_t ctr=0;
   for(auto&& le : docs) {
      Rcpp::IntegerMatrix temp((SEXP)le);
      ar_docs[ctr] = arma::imat(temp.begin(), temp.nrow(), temp.ncol(), false);
//      std::cout<<ar_docs[ctr]<<std::endl;
      ++ctr;
   }
   
   // No deep copies
   Rcpp::NumericVector beta_idx(beta_idx_);
   arma::vec ar_beta_idx(beta_idx.begin(), beta_idx.size(), false); // creates col vec 
   
   // No deep copies (no arma needed)
   Rcpp::LogicalVector update_mu(update_mu_);
   
   // No deep copies
   Rcpp::List beta(beta_);
   std::map<std::size_t, arma::mat> ar_beta;
   ctr=0;
   for(auto&& le : beta) {
      Rcpp::NumericMatrix temp((SEXP)le);
      ar_beta[ctr] = arma::mat(temp.begin(), temp.nrow(), temp.ncol(), false);
//      std::cout<<ar_beta[ctr]<<std::endl;
      ++ctr;
   }
   
   // No deep copies
   Rcpp::NumericMatrix lambda_old(lambda_old_);
   arma::mat ar_lambda_old(lambda_old.begin(), lambda_old.nrow(), lambda_old.ncol(), false);
   
   // No deep copies
   Rcpp::NumericVector mu(mu_);
   arma::vec ar_mu(mu.begin(), mu.size(), false); // creates col vec
   
   // No deep copies
   Rcpp::NumericMatrix sigma(sigma_);
   arma::mat ar_sigma(sigma.begin(), sigma.nrow(), sigma.ncol(), false); 
   arma::mat ar_sigmainv = arma::zeros<arma::mat>(ar_sigma.n_rows, ar_sigma.n_cols);
   
   // No deep copies (no arma needed)
   Rcpp::CharacterVector method(method_);
   Rcpp::LogicalVector verbose(verbose_);
   
   // Set up variables needed during estep comp
   std::size_t V = ar_beta[0].n_cols;
   std::size_t K = ar_beta[0].n_rows;
   std::size_t N = docs.length();
   std::size_t A = beta.length();
   
   
   // Print some debug info
   std::cout<<"Number of documents.length(): "<<docs.length()<<std::endl;
   std::cout<<"Number of beta indices: "<<beta_idx.length()<<std::endl;
   std::cout<<"Update mu: "<<update_mu<<std::endl;
   std::cout<<"Number of beta elements "<<beta.length()<<std::endl;
   std::cout<<"Number of old lambda elements "<<lambda_old.length()<<std::endl;
   std::cout<<"Number of elements of mu "<<mu.length()<<std::endl;
   std::cout<<"Number of elements of sigma "<<sigma.length()<<std::endl;
   std::cout<<"What method: "<<method<<std::endl;
   std::cout<<"Verbose: "<<verbose<<std::endl;
   std::cout<<"V: "<<V<<" K: "<<K<<" N: "<<N<<" A: "<<A<<std::endl;
   
   arma::vec ar_zero_vec = arma::zeros<arma::vec>(K-1);
   arma::mat ar_sigma_ss = diagmat(ar_zero_vec);
   
   std::map<std::size_t, arma::mat> ar_beta_ss;
   for(std::size_t i=0; i<A; ++i)
      ar_beta_ss[i] = arma::zeros<arma::mat>(K,V);
   
   arma::vec ar_bound = arma::zeros<arma::vec>(N);
   arma::mat ar_lambda_new = arma::zeros<arma::mat>(ar_lambda_old.n_rows, ar_lambda_old.n_cols);
   
   // Compute sigma inverse
   bool worked = arma::chol(ar_sigmainv, ar_sigma);
   if(!worked) {
      // Make matrix positive definite by diagonal dominance
      arma::vec dvec = ar_sigma.diag();
      //find the magnitude of the diagonal 
      arma::vec magnitudes = sum(abs(ar_sigma), 1) - abs(dvec);
      //iterate over each row and set the minimum value of the diagonal to be the magnitude of the other terms
      int Km1 = dvec.size();
      for(int j=0; j < Km1;  j++) {
         if(arma::as_scalar(dvec(j)) < arma::as_scalar(magnitudes(j))) dvec(j) = magnitudes(j); //enforce diagonal dominance 
      }
      //overwrite the diagonal of the hessian with our new object
      ar_sigma.diag() = dvec;
      //that was sufficient to ensure positive definiteness so we now do cholesky
      ar_sigmainv = arma::chol(ar_sigma);
   } else {
      std::cout<<"No diagonal dom adjustment needed..."<<std::endl;
   }
   
   std::cout<<"Solving the first optimization problem..."<<std::endl;
   
   optimfn objlhood;
   optimgr gradlhood;
   
   for(arma::uword d=0; ar_docs.size(); ++d) {
      arma::vec init = ar_lambda_old.row(d);
      problem_data* pr = problem_data_alloc();
      switch (hash_opt_method(collapse(method))) {
      case BFGS:
         vmmin(K-1, init.memptr(), 0.0, objlhood, gradlhood, 500, 0, NULL, 1e-8, 1e-8, 0, pr, );
         break;
      case LBFGSB:
         
         break;
      case CG:
         
         break;
         /*   case default:
          break; */
         
      }
   }
   
   /*cg_descent(ar_lambda_old.colptr(0), 
              K-1, 
              NULL, 
              NULL, 
              1.e-8, 
              objlhood, 
              gradlhood, 
              objgradlhood, 
              K-1, 
              V, 
              ar_beta[0].memptr(), 
              (double*) ar_docs[1].memptr(), 
              ar_mu.memptr(), 
              ar_sigmainv.memptr(), 
              NULL);
    */
   
   /*for(std::size_t i=0; i<N; ++i) {
      cg_descent(ar_lambda_old.colptr(i), 
                 K-1, 
                 NULL, 
                 NULL, 
                 1.e-8, 
                 objlhood, 
                 gradlhood, 
                 objgradlhood, 
                 V, 
                 K, 
                 ar_beta[0].memptr(), 
                 (double*) ar_docs[1].memptr(), 
                 ar_mu.memptr(), 
                 ar_sigmainv.memptr(), 
                 NULL);
   }*/


}


