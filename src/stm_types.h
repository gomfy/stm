//
// Created by gomfy on 01/10/20.
//

#pragma once

#include "ExactSum.h"
#include "Likelihood.h"

extern "C" {
#include "cg_user.h"
}

enum opt_method {
  BFGS,
  LBFGSB,
  CG
};

typedef struct problem_data {
  arma::mat* beta;
  arma::vec* doc_ct;
  arma::vec* mu;
  arma::mat* siginv;
} problem_data;
