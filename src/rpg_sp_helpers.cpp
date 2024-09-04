#include <Rcpp.h>
#include "rpg_sp_helpers.h"
using namespace Rcpp;

// Compute delta as per Lemma 14 (Windle 2014 arXiv) or
// Lemma 2.22 (Windle 2013 thesis)
FD get_delta(double x, double mid) {

  FD delta;

  if (x >= mid) {

    // Compute the right side of delta
    delta.val = log(x) - log(mid);
    delta.der = 1.0 / x;

  } else {

    // Compute the left side of delta
    delta.val = 0.5 * (1 - 1.0 / x) - 0.5 * (1 - 1.0 / mid);
    delta.der = 0.5 / (x*x);

  }

  return delta;
}

// Find t(x) such that phi(x) = K(t) - tx
FD get_phi(double x, double z) {

  FD phi;

  // Determine the phi function by finding t(x)
  double v = v_eval(x);       // 2u from Windle paper/thesis
  double u = 0.5 * v;         // u from Windle paper/thesis
  double t = u + 0.5 * z*z;   // t from Windle paper/thesis

  // Evaluate the phi function and its derivative at t(x)
  phi.val = logcosh(z) - log_cos_rt(v) - t*x; // Windle 2014 arXiv Fact 9.4
  phi.der = -1.0 * t; // Windle 2014 arXiv Fact 9.5

  return phi;
}

// Calculate tangent line to eta at x
Line get_eta_tangent(double x, double z, double mid) {

  Line L;
  FD phi, delta, eta;

  // Compute phi and delta at x
  phi = get_phi(x, z);
  delta = get_delta(x, mid);

  // Evaluate eta and its derivative at x
  eta.val = phi.val - delta.val;
  eta.der = phi.der - delta.der;

  // Calculate slope and intercept for tangent line at x
  L.slope = eta.der;
  L.intercept = eta.val - eta.der * x;

  return L;
}

// Calculate the saddlepoint approximation given at the top of pg. 69 (Windle
// 2013 thesis) or the bottom of pg. 17 (Windle 2014 arXiv).
double sp_approx(double x, double b, double z) {

  // Find v
  double v = v_eval(x);

  // Occasionally (very rarely) the starting point is bad and we don't converge.
  // When this happens, v_eval returns -999 and we try again.
  if (v == -999) return v;

  // Define constants
  double u  = 0.5 * v;
  double z2 = z * z;
  double t  = u + 0.5 * z2;

  // Compute phi
  double phi = logcosh(z) - log_cos_rt(v) - t*x;

  // Evaluate second derivative of K at t(x)
  double K2  = 0.0;
  if (fabs(v) >= 1e-6)
    K2 = x*x + (1-x) / v;
  else
    K2 = x*x - 1/3 - (2/15) * v;

  // Compute SP approximation
  double log_spa = 0.5 * log(0.5 * b / MATH_PI) - 0.5 * log(K2) + b * phi;
  return exp(log_spa);
}

// More stable log of cosh
double logcosh(double x) {
  if (x < 0) {
    return -log(2) - x + log1pexp(2*x);
  } else {
    return -log(2) + x + log1pexp(-2*x);
  }
}

// Compute log of cos(sqrt(v))
// Windle 2013 thesis Fact 9 proof (and take the log)
double log_cos_rt(double v) {
  double r = sqrt(fabs(v));
  double y;
  if (v >= 0) y = log(cos(r));
  else y = logcosh(r);
  return y;
}

// Compute tan(sqrt(x)) / sqrt(x)
// Windle 2013 thesis Fact 2.19
double y_func(double v) {

  double tol = 1e-6;
  double y   = 0.0;
  double r   = sqrt(fabs(v));

  // If v > 0, compute directly
  // If v < 0, compute using hyperbolic tangent
  // If v close to 0, use a Taylor expansion around s = 0
  if (v > tol)
    y = tan(r) / r;
  else if (v < -1*tol)
    y = tanh(r) / r;
  else
    y = 1 + (1/3) * v + (2/15) * v * v + (17/315) * v * v * v;
  return y;
}

