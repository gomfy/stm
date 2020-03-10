
#pragma once

#include<armadillo>

class Likelihood
{
  std::size_t K;
  std::size_t V;
  arma::vec* expeta;
  arma::mat* beta;
  arma::vec* doc_ct;
  arma::vec* mu;
  arma::mat* siginv;
  
public:
  Likelihood(): K(0), V(0), beta(NULL), doc_ct(NULL), mu(NULL), siginv(NULL) {};
  Likelihood(std::size_t K_, std::size_t V_): K(K_), V(V_) {
    expeta = new arma::vec(K);
    expeta->fill(1.0);
    beta = new arma::mat(K,K);
    beta->fill(0.0);
    doc_ct = new arma::vec(V);
    doc_ct->fill(0.0);
    mu = new arma::vec(K-1);
    mu->fill(0.0);
    siginv = new arma::mat(K-1,K-1);
    siginv->fill(0.0);
  };
  ~Likelihood() {};

};