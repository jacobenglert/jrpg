/**
 * This file implements the Polya-gamma sampler PG(1,z).
 * This is a C++ implementation of Algorithm 6 in PhD thesis of Jesse
 * Bennett Windle, 2013
 * URL: https://repositories.lib.utexas.edu/bitstream/handle/2152/21842/WINDLE-DISSERTATION-2013.pdf?sequence=1
 *
 * References:
 *
 *   Jesse Bennett Windle
 *   Forecasting High-Dimensional, Time-Varying Variance-Covariance Matrices
 *   with High-Frequency Data and Sampling Polya-Gamma Random Variates for
 *   Posterior Distributions Derived from Logistic Likelihoods
 *   PhD Thesis, 2013
 *
 *   Damien, P. & Walker, S. G. Sampling Truncated Normal, Beta, and Gamma Densities
 *   Journal of Computational and Graphical Statistics, 2001, 10, 206-215
 *
 *   Chung, Y.: Simulation of truncated gamma variables
 *   Korean Journal of Computational & Applied Mathematics, 1998, 5, 601-610
 *
 * (c) Copyright Enes Makalic and Daniel F Schmidt, 2018
 *
 * ! Changes for current work: Modified to work with RcppEigen
 */

#include <RcppArmadillo.h>
#include "rpg.h"
#include "rpg_sp_helpers.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using std::pow;

//' Multiple draw polya gamma latent variable for var c[i] with size b[i]
//'
//' @param b vector of binomial sizes
//' @param z vector of parameters
//' @returns Eigen::VectorXd
//' @export
// [[Rcpp::export]]
arma::vec jrpg(arma::vec b, arma::vec z) {
  int n = z.size();
  arma::vec y(n);

  for (int i = 0; i < n; ++i) {

    if (b(i) >= 170) {

      // Normal approximation
      y(i) = rpg_na(b(i), z(i));

    } else if (b(i) > 13) {

      y(i) = rpg_sp(b(i), z(i), 5000);

      if (y(i) == -999) {
        // Devroye method
        int bint = std::floor(b(i));
        y(i) = 0.0;
        y(i) += rpg_devroye(bint, z(i));

        // Add sum of gammas to cover remainder
        if (b(i) > bint) {
          y(i) += rpg_gamma(b(i) - bint, z(i));
        }
      }

    } else if (b(i) >= 1) {

      // Devroye method
      int bint = std::floor(b(i));
      y(i) = 0.0;
      y(i) += rpg_devroye(bint, z(i));

      // Add sum of gammas to cover remainder
      if (b(i) > bint) {
        y(i) += rpg_gamma(b(i) - bint, z(i));
      }

    } else if (b(i) > 0) {

      // Sum of gammas
      y(i) = rpg_gamma(b(i), z(i));

    } else if (b(i) == 0) {
      y(i) = 0.0;
    }

  }

  return y;
}

// SAMPLER FUNCTIONS

// Devroye Method: Sample PG(b,z)
// Usage Case: 1 < b < 13 and b is integer
// [[Rcpp:export]]
double rpg_devroye(int b, double z) {

  double sum = 0.0;

  double z_2 = 0.5 * fabs(z);
  double K = z_2*z_2 / 2.0 + MATH_PI2 / 8.0;
  double r = ratio(z_2);

  for (int j = 0; j < b; ++j) {
    sum += rpg_devroye_1(z_2, r, K);
  }

  return sum;
}

// Devroye Method: Sample PG(1,z)
// Based on Algorithm 6 in PhD thesis of Jesse Bennett Windle, 2013
// URL: https://repositories.lib.utexas.edu/bitstream/handle/2152/21842/WINDLE-DISSERTATION-2013.pdf?sequence=1
double rpg_devroye_1(double z, double r, double K) {

  // PG(b, z) = 0.25 * J*(b, z/2)

  double t = MATH_2_PI;
  double u, X;

  // Main sampling loop; page 130 of the Windle PhD thesis
  while (1) {
    // Step 1: Sample X from mixture proposal distribution g(x|z)
    u = R::runif(0.0,1.0);
    if (u < r) {
      // truncated exponential
      X = t + exprnd(1.0) / K;
    } else {
      // truncated Inverse Gaussian
      X = rrtinvgauss(1.0 / z, t);
    }

    // Step 2: Iteratively calculate Sn(X|z), starting at S1(X|z), until U ? Sn(X|z) for an odd n or U > Sn(X|z) for an even n
    int i     = 1;
    double Sn = aterm(0, X, t);
    double U  = R::runif(0.0,1.0) * Sn;
    int asgn  = -1;
    bool even = false;

    while (1) {
      Sn = Sn + asgn * aterm(i, X, t);

       // Accept if n is odd
       if (!even && (U <= Sn)) {
         X = X * 0.25; // Remember to divide by 4
         return X;
       }

       // Return to step 1 if n is even
       if (even && (U > Sn)) {
         break;
       }

       even = !even;
       asgn = -asgn;
       i++;
    }
  }
  return X;
}

// Normal Approximation Method
// Usage Case: b > 170
// [[Rcpp:export]]
double rpg_na(double b, double z) {

  double m1, m2;

  z = 0.5 * fabs(z);

  // Compute the first and second moments
  if (z > 1e-12) {
    m1 = b * tanh(z) / z; // off by scale of 1/4
    m2 = b * (b + 1) * pow(tanh(z) / z, 2) + b * ((tanh(z) - z) / pow(z, 3)); // off by scale of 1/16
  } else{
    m1 = b * (1 - (1.0/3) * pow(z,2) + (2.0/15) * pow(z,4) - (17.0/315) * pow(z,6));
    m2 = b * (b+1) * pow(1 - (1.0/3) * pow(z,2) + (2.0/15) * pow(z,4) - (17.0/315) * pow(z,6), 2) +
      b * ((-1.0/3) + (2.0/15) * pow(z,2) - (17.0/315) * pow(z,4));
  }

  // Generate a random variate based on a normal approximation
  return R::rnorm(0.25 * m1, 0.25 * sqrt(m2 - m1 * m1));

}

// Truncated Sum of Gammas Method
// Usage Case: 0 < b < 1 (and also for remainder term for 1 < b < 13)
// [[Rcpp::export]]
double rpg_gamma(double b, double z) {
  int trunc = 200;
  double z2 = 0.125 * z * z;
  double sum = 0.0;
  for (int i = 0; i < trunc; ++i) {
    sum += R::rgamma(b, 1) / (MATH_PI2_2 * (i + 0.5) * (i + 0.5) + z2);
  }
  return 0.25 * sum;
}

// Saddlepoint Approximation Method
// Usage Case: 13 < b < 170
double rpg_sp(double b, double z, int maxiter = 2000) {

  if (b < 1) {
    stop("Saddlepoint approximation only valid for b >= 1.");
  }

  z = 0.5 * fabs(z);

  // Select reference points
  double xl = y_eval(-1*z*z);   // Left point (mode of phi)
  double xc = 1.1 * xl;                    // Mid point
  double xr = 1.2 * xl;                    // Right point

  // Calculate inflation constants
  // I have no idea how this works - borrowed directly from BayesLogit
  double vxc  = v_eval(xc);
  double K2xc = 0.0;

  if (fabs(vxc) >= 1e-6)
    K2xc = xc*xc + (1-xc) / vxc;
  else
    K2xc = xc*xc - 1/3 - (2/15) * vxc;

  double xc2 = xc * xc;
  double al = xc2*xc / K2xc;
  double ar = xc2    / K2xc;

  // Calculate left and right tangent lines (slopes and intercepts)
  Line Ll, Lr;
  tangent_to_eta(xl, z, xc, Ll);
  tangent_to_eta(xr, z, xc, Lr);

  // Extract slopes and intercepts
  double rho_l = -1.0 * Ll.slope;
  double rho_r = -1.0 * Lr.slope;
  double b_l = Ll.intercept;
  double b_r = Lr.intercept;

  // Constants
  double lcn = 0.5 * log(0.5 * b / MATH_PI);
  double sqrt_rho_l = sqrt(2 * rho_l);

  double k_l, k_r, w_l, w_r, p_l;

  // (Error in BayesLogit? Changed to -0.5 * log(a) in both and + b * log(xc) in k_r)
  // double k_l = exp(0.5 * log(al) + b * (0.5 / xc + b_l - sqrt_rho_l));
  k_l = exp(0.5 * log(al) - b * sqrt_rho_l + b * b_l + 0.5 * b * 1.0 / xc);
  // k_r = exp(lcn + 0.5 * log(ar) - b * (log(xc) - b_r + log(b * rho_r)) + lgamma(b));
  k_r = exp(0.5 * log(ar) + lcn + (- b * log(b * rho_r) + b * b_r - b * log(xc) + lgamma(b)));

  // Weights
  w_l = k_l * pinvgauss(xc, 1.0 / sqrt_rho_l, b);
  w_r = k_r * (1.0 - R::pgamma(xc, b, 1.0 / (b * rho_r), true, false));
  p_l = w_l / (w_l + w_r);

  // Sample
  int iter = 0;
  double X = 2.0;
  double k = 0.0;
  bool go = true;
  while (go && iter < maxiter) {

    iter++;

    double phi_ev;
    if (R::runif(0.0,1.0) < p_l) {

      // Generate from right truncated inverse Gaussian
      X = rrtinvgauss(1.0 / sqrt_rho_l, b, xc);
      phi_ev = b * (b_l - rho_l * X) + 0.5 * b * ((1.0 - 1.0 / X) - (1.0 - 1.0 / xc));
      k = exp(0.5 * log(al) + lcn - 1.5 * log(X) + phi_ev);

    }
    else {

      // Generate from left-truncated gamma
      X = rltgamma(b, b * rho_r, xc);
      phi_ev = b * (b_r - rho_r * X) + b * (log(X) - log(xc));
      k = exp(0.5 * log(ar) + lcn + phi_ev) / X;

    }

    double spa = sp_approx(X, b, z);
    if (spa == -999) return spa;


    if (k * R::runif(0.0,1.0) < spa) go = false;

    // if (k * R::runif(0.0,1.0) < spa)
    //   break;

  }

  return b * 0.25 * X;

}
