// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
#include <RcppParallel.h>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;
using namespace RcppParallel;
using namespace std;

#include <boost/math/distributions/poisson.hpp>
//' @importFrom RcppParallel RcppParallelLibs

typedef boost::math::poisson_distribution<float, 
            boost::math::policies::policy<boost::math::policies::discrete_quantile<boost::math::policies::real> > > 
        poisson_real_quantile;

// set seed
// [[Rcpp::export]]
void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);  
}

// [[Rcpp::export]]
int draw_sample_cdf(
    const float qmin, 
    const float qmax, 
    const float mu
) {
  try {
    if (qmax == qmin) return floor(qmin); 
    
    int D = floor(qmax+1) - floor(qmin) + 1; // number of elements 
    arma::vec vals = arma::linspace<arma::vec>(floor(qmin), floor(qmax+1), D); 
    if (D == 1) return floor(qmin);
    
    arma::vec probs = arma::zeros(D);
    poisson_real_quantile p(mu);
    for (int j = 0; j < D; j++) {
        probs(j) = pdf(p, vals(j));
    }
    probs(0) *= floor(1 + qmin) - qmin; 
    probs(D-1) *= 1 - (floor(1 + qmax) - qmax);
    probs = arma::normalise(probs, 1);

    // set_seed(2.0);     
    // draw from the distribution
    return round(RcppArmadillo::sample(vals, 1L, false, probs)[0]);
  
  } catch (...) {
      return -1;
  }
}

// [[Rcpp::export]]
int icdfpois_scalar(
    const int k, 
    const float mu_source,
    const float mu_target
) {
    float cdf_ceil=0.99999;
    try {
        float cdf_min, cdf_max, qmin, qmax;
        if (k == 0) {
            qmin = 0;
        } else {
            cdf_min = min(cdf(poisson_real_quantile(mu_source), k-1), cdf_ceil);
            // poisson_real_quantile(mu_source) creates a poisson distribution with rate (lambda) of mu_source
            // k - 1 is the observed
            // cdf(poisson_real_quantile(mu_source), k - 1) gives the area under the curve of the full distribution
            qmin = boost::math::quantile(poisson_real_quantile(mu_target), cdf_min) + 1; 
            // poisson_real_quantile(mu_target) creates a poisson distribution with rate (lambda) of mu_target
            // qmin is the value of the distribution that gives the same area as cdf_min()
        }

        // upper bound on range        
        cdf_max = min(cdf(poisson_real_quantile(mu_source), k), cdf_ceil);
        qmax = boost::math::quantile(poisson_real_quantile(mu_target), cdf_max);
        if (qmin > qmax) qmax = qmin;

            // draw value from distribution
        return draw_sample_cdf(qmin, qmax, mu_target);
    } catch (...) {
        return -1;
    }
}

struct icdfpoisStruct : public Worker
{
    // source matrices
    const RMatrix<int> k;
    const RMatrix<double> mu_source;
    const RMatrix<double> mu_target;

    // destination matrix
    RMatrix<int> output;

    // Constructor
    icdfpoisStruct(
        const IntegerMatrix k, 
        const NumericMatrix mu_source, 
        const NumericMatrix mu_target, 
        IntegerMatrix output
    ) 
        : k(k), mu_source(mu_source), mu_target(mu_target), output(output) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (int i=begin; i<end; i++) {
            output[i] = icdfpois_scalar(k[i], mu_source[i], mu_target[i]);
        }
    }
};


// [[Rcpp::export]]
IntegerMatrix parallelicdfPoisson(
    IntegerMatrix k, 
    NumericMatrix mu_source, 
    NumericMatrix mu_target
) {
    IntegerMatrix output(k.nrow(), k.ncol());
    icdfpoisStruct icdfpoisInstance(k, mu_source, mu_target, output);
    parallelFor(0, k.length(), icdfpoisInstance);
    return output;
}
