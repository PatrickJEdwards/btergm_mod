#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;


// Negative log pseudo-likelihood for ERGM count models, using precomputed quantities.
//
// [[Rcpp::export]]
double ergmCntNLPL_cpp(NumericVector coef,
                       List obj,
                       int rtype,
                       NumericVector rparam) {
  double lpl = 0.0, ils, pd, reg = 0.0;
  int i, j, k;
  
  // Make sure that we were passed a legitimate input, lest we crash.
  if (!obj.inherits("ERGMCntPrep")) {
    stop("Must be called with an ERGMCntPrep object.");
  }
  
  // Coerce the importance weights.
  NumericVector iw = as<NumericVector>(obj["iw"]);
  
  // Add up the negative log pseudo-likelihood.
  for (i = 0; i < iw.size(); i++) {
    ils = 0.0;
    
    NumericVector rmr = as<NumericVector>(as<List>(obj["rmr"])[i]);
    NumericMatrix cs = as<NumericMatrix>(as<List>(obj["cs"])[i]);
    NumericVector ycwt = as<NumericVector>(as<List>(obj["ycwt"])[i]);
    
    // Find the inverse log sum.
    for (j = 0; j < rmr.size(); j++) {
      pd = 0.0;
      
      for (k = 0; k < coef.size(); k++) {
        pd += cs(j, k) * coef[k];
      }
      
      if (j == 0) {
        ils = rmr[j] + pd + std::log(ycwt[j]);
      } else {
        ils = R::logspace_add(ils, rmr[j] + pd + std::log(ycwt[j]));
      }
    }
    
    // Add to the total, multiplying by inverse inclusion weight.
    lpl -= iw[i] * ils;
  }
  
  // Compute the regularization penalty, if any.
  if (rtype == 1) {                              // L1 regularization
    for (i = 0; i < coef.length(); i++) {
      reg += std::fabs(coef[i]);
    }
    reg *= rparam[0];
    
  } else if (rtype == 2) {                       // L2 regularization
    for (i = 0; i < coef.length(); i++) {
      reg += coef[i] * coef[i];
    }
    reg *= rparam[0];
    
  } else if (rtype == 3) {                       // pseudo-Huber regularization
    for (i = 0; i < coef.length(); i++) {
      reg += rparam[1] *
        (std::sqrt(1.0 + coef[i] * coef[i] / (rparam[1] * rparam[1])) - 1.0);
    }
    reg *= rparam[0];
  }
  
  // Return the negative log pseudo-likelihood plus regularization penalty.
  return -lpl + reg;
}


// Negative log pseudo-likelihood and derivatives for ERGM count models,
// using precomputed quantities. Output is formatted for the trust package.
//
// [[Rcpp::export]]
List ergmCntNLPLDeriv_cpp(NumericVector coef,
                          List obj,
                          int rtype,
                          NumericVector rparam) {
  double lpl = 0.0, ils, pd, reg = 0.0;
  int i, j, k, l, p = coef.size();
  
  NumericVector gr(p), igr(p);             // Gradient of the negative LPL.
  NumericMatrix hess(p, p), ihess(p, p);   // Hessian of the negative LPL.
  
  // Make sure that we were passed a legitimate input, lest we crash.
  if (!obj.inherits("ERGMCntPrep")) {
    stop("Must be called with an ERGMCntPrep object.");
  }
  
  // Coerce the importance weights.
  NumericVector iw = as<NumericVector>(obj["iw"]);
  
  // Add up the likelihood, gradient, and Hessian.
  for (i = 0; i < iw.size(); i++) {
    ils = 0.0;
    
    NumericVector rmr = as<NumericVector>(as<List>(obj["rmr"])[i]);
    NumericMatrix cs = as<NumericMatrix>(as<List>(obj["cs"])[i]);
    NumericVector ycwt = as<NumericVector>(as<List>(obj["ycwt"])[i]);
    
    for (j = 0; j < p; j++) {
      igr[j] = 0.0;
      
      for (k = 0; k < p; k++) {
        ihess(j, k) = 0.0;
      }
    }
    
    // First pass: compute inverse log sum to avoid overflow.
    for (j = 0; j < rmr.size(); j++) {
      pd = 0.0;
      
      for (k = 0; k < p; k++) {
        pd += cs(j, k) * coef[k];
      }
      
      if (j == 0) {
        ils = rmr[j] + pd + std::log(ycwt[j]);
      } else {
        ils = R::logspace_add(ils, rmr[j] + pd + std::log(ycwt[j]));
      }
    }
    
    // Second pass: derivative elements.
    for (j = 0; j < rmr.size(); j++) {
      pd = 0.0;
      
      for (k = 0; k < p; k++) {
        pd += cs(j, k) * coef[k];
      }
      
      double weight = std::exp(rmr[j] + pd - ils) * ycwt[j];
      
      for (k = 0; k < p; k++) {
        // Local gradient contribution.
        igr[k] -= cs(j, k) * weight;
        
        // Local Hessian contribution.
        for (l = 0; l < p; l++) {
          ihess(k, l) += cs(j, k) * cs(j, l) * weight;
        }
      }
    }
    
    // Add to total, multiplying by inverse inclusion weight.
    lpl -= iw[i] * ils;
    
    for (j = 0; j < p; j++) {
      gr[j] -= iw[i] * igr[j];
      
      for (k = 0; k < p; k++) {
        hess(j, k) -= iw[i] * (igr[j] * igr[k] - ihess(j, k));
      }
    }
  }
  
  // Compute the regularization penalty, gradient, and Hessian, if any.
  if (rtype == 1) {                              // L1 regularization
    for (i = 0; i < p; i++) {
      reg += std::fabs(coef[i]);
      gr[i] += (coef[i] > 0.0 ? 1.0 : (coef[i] < 0.0 ? -1.0 : 0.0)) * rparam[0];
    }
    reg *= rparam[0];
    
  } else if (rtype == 2) {                       // L2 regularization
    for (i = 0; i < p; i++) {
      reg += coef[i] * coef[i];
      gr[i] += 2.0 * coef[i] * rparam[0];
      hess(i, i) += 2.0 * rparam[0];
    }
    reg *= rparam[0];
    
  } else if (rtype == 3) {                       // pseudo-Huber regularization
    for (i = 0; i < p; i++) {
      double delta = rparam[1];
      double denom = std::sqrt(1.0 + coef[i] * coef[i] / (delta * delta));
      
      reg += delta * (denom - 1.0);
      
      gr[i] += rparam[0] * coef[i] / (delta * denom);
      
      // NOTE: The original replication code appears to use rparam[i] here.
      // That is unsafe because rparam usually has length 2. This uses delta = rparam[1].
      hess(i, i) += rparam[0] * delta /
        ((delta * delta + coef[i] * coef[i]) * denom);
    }
    reg *= rparam[0];
  }
  
  // Return penalized negative log pseudo-likelihood, gradient, and Hessian.
  List out = List::create(
    Named("value") = -lpl + reg,
    Named("gradient") = gr,
    Named("hessian") = hess
  );
  
  return out;
}