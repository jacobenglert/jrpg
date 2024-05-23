#include <Rcpp.h>
#include "rpg_devroye_helpers.h"
using namespace Rcpp;

// Compute the ratio used to determine how to sample J*(1,z)
// Used for Devroye sampler
double ratio(double z) { //, double K) {

  // Point on the intersection IL = [0, 4/ log 3] and IR = [(log 3)/pi^2, \infty)
  double t = MATH_2_PI;

  // Compute the ratio q / (p + q)
  // (derived from scratch; derivation is not in the original paper)
  double K = z*z / 2.0 + MATH_PI2 / 8.0;

  double x0 =  (double)std::log(K) + K * t;
  double b = (double)std::sqrt(1.0 / t) * (t * z - 1);
  double a = (double)std::sqrt(1.0 / t) * (t * z + 1) * -1.0;

  double xb = x0 - z + R::pnorm(b, 0.0, 1.0, true, true);
  double xa = x0 + z + R::pnorm(a, 0.0, 1.0, true, true);

  double p_over_q = 4.0 / MATH_PI * (std::exp(xb) + std::exp(xa));

  // double K    = z*z/2.0 + MATH_PI2/8.0;
  // double logA = (double)std::log(4.0) - MATH_LOG_PI - z;
  // double logK = (double)std::log(K);
  // double Kt   = K * t;
  // double w    = (double)std::sqrt(MATH_PI_2);
  //
  // double logf1    = logA + R::pnorm(w*(t*z - 1),0.0,1.0,1,1) + logK + Kt;
  // double logf2    = logA + 2*z + R::pnorm(-w*(t*z+1),0.0,1.0,1,1) + logK + Kt;
  // double p_over_q = (double)std::exp(logf1) + (double)std::exp(logf2);

  double ratio = 1.0 / (1.0 + p_over_q);

  return ratio;
}

// Function a_n(x) defined in equations (12) and (13) of
// Bayesian inference for logistic models using Polya-Gamma latent variables
// Nicholas G. Polson, James G. Scott, Jesse Windle
// arXiv:1205.0310
//
// Also found in the PhD thesis of Windle (2013) in equations
// (2.14) and (2.15), page 24
double aterm(int n, double x, double t) {
  double a = 0;
  if (x <= t) {
    a = MATH_LOG_PI + std::log(n + 0.5) + 1.5*(M_LOG_2_PI-std::log(x)) - 2*(n + 0.5)*(n + 0.5)/x;
  } else {
    a = MATH_LOG_PI + std::log(n + 0.5) - x * MATH_PI2_2 * (n + 0.5)*(n + 0.5);
  }

  return exp(a);
}
