#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <cmath>
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat solve_betas(
    const arma::mat& design_full, 
    const arma::mat& lambda_mat, 
    const arma::mat& Ynorm
    ) {

    return arma::inv(design_full.t() * design_full + lambda_mat) * design_full.t() * Ynorm; 
}

// [[Rcpp::export]]
arma::mat merge_redundant_clusters(const arma::mat & R, float thresh = 0.8) {
  arma::mat cor_res = arma::cor(R.t());
  int K = R.n_rows;
  
  // define equivalence classes  
  arma::uvec equiv_classes = arma::linspace<arma::uvec>(0, K - 1, K);
  int new_idx;
  for (int i = 0; i < K - 1; i++) {
    for (int j = i + 1; j < K; j++) {
      if (cor_res(i, j) > thresh) {
          new_idx = std::min(equiv_classes(i), equiv_classes(j));
          equiv_classes(i) = new_idx;
          equiv_classes(j) = new_idx;
      }
    }
  }

  // sum over equivalence classes
  arma::uvec uclasses = arma::unique(equiv_classes);
  arma::mat R_new = arma::zeros<arma::mat>(uclasses.n_elem, R.n_cols); 
  for (int i = 0; i < R_new.n_rows; i++) {
      arma::uvec idx = arma::find(equiv_classes == uclasses(i)); 
      R_new.row(i) = arma::sum(R.rows(idx), 0);
  }
  return R_new;  
}

