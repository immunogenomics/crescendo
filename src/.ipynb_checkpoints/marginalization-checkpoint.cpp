#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <cmath>
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

// Returning dense matrix is slightly faster and better for memory
// [[Rcpp::export]]
arma::mat solve_mu1(const arma::sp_mat& design_full, const arma::mat& betas) {
    arma::mat result(design_full * betas);
    return(result);
}

// [[Rcpp::export]]
void solve_mu2(arma::mat& expected_logcounts_marginalized, 
               const arma::mat& R, const arma::mat& sigmas, const arma::mat& mu){
    expected_logcounts_marginalized = expected_logcounts_marginalized + R * (.5 * sigmas % sigmas + mu);
}

