#include <Rcpp.h>
#include "rpg_sp_helpers.h"
using namespace Rcpp;

// Compute delta as per Lemma 14 / 2.22
void delta_func(double x, double mid, FD& delta) {

  if (x >= mid) {

    // Compute the right side of delta
    delta.val = log(x) - log(mid);
    delta.der = 1.0 / x;

  } else {

    // Compute the left side of delta
    delta.val = 0.5 * (1 - 1.0 / x) - 0.5 * (1 - 1.0 / mid);
    delta.der = 0.5 / (x*x);

  }
}

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
void phi_func(double x, double z, FD& phi) {

  // Determine the phi function by finding t(x)
  double v = v_eval(x);       // 2u from paper
  double u = 0.5 * v;         // u from paper
  double t = u + 0.5 * z*z;   // t from paper

  // Evaluate the phi function and its derivative at t(x)
  phi.val = log(cosh(fabs(z))) - log(cos_rt(v)) - t * x; // Fact 9.4
  phi.der = -1.0 * t; // Fact 9.5

  // return v;

}

FD get_phi(double x, double z) {

  FD phi;

  // Determine the phi function by finding t(x)
  double v = v_eval(x);       // 2u from paper
  double u = 0.5 * v;         // u from paper
  double t = u + 0.5 * z*z;   // t from paper

  // Evaluate the phi function and its derivative at t(x)
  // phi.val = log(cosh(fabs(z))) - log(cos_rt(v)) - t * x; // Fact 9.4
  phi.val = logcosh(z) - log_cos_rt(v) - t*x;
  phi.der = -1.0 * t; // Fact 9.5

  // return v;

  return phi;
}

// Calculate tangent line to eta
void tangent_to_eta(double x, double z, double mid, Line& tl) {

  FD phi, delta, eta;
  // double v;

  // Compute v and update phi in place
  // v = phi_func(x, z, phi);
  phi_func(x, z, phi);

  // Compute and update delta in place
  delta_func(x, mid, delta);

  // Evaluate eta and its derivative at x
  eta.val = phi.val - delta.val;
  eta.der = phi.der - delta.der;

  // Calculate and update slope/intercept for line
  tl.slope = eta.der;
  tl.intercept = eta.val - eta.der * x;

  // return v;
}

Line get_eta_tangent(double x, double z, double mid) {

  Line L;
  FD phi, delta, eta;

  // Compute phi and delta at x
  phi = get_phi(x, z);
  delta = get_delta(x, mid);

  // Evaluate eta and its derivative at x
  eta.val = phi.val - delta.val;
  eta.der = phi.der - delta.der;

  // Calculate and update slope/intercept for line
  L.slope = eta.der;
  L.intercept = eta.val - eta.der * x;

  return L;
}



// Compute cos(sqrt(v))
// Windle thesis Fact 9 proof
double cos_rt(double v) {
  double y   = 0.0;
  double r   = sqrt(fabs(v));
  if (v >= 0)
    y = cos(r);
  else
    y = cosh(r);
  return y;
}





// Saddlepoint approximation given at the top of pg. 69
double sp_approx(double x, double b, double z) {

  double v = v_eval(x);
  // if (v == -999) return v;
  double u  = 0.5 * v;
  double z2 = z * z;
  double t  = u + 0.5 * z2;

  // Compute phi
  // double phi = log(cosh(z)) - log(cos_rt(v)) - t * x;
  double phi = logcosh(z) - log_cos_rt(v) - t*x;

  // Compute second derivative
  double K2  = 0.0;
  if (fabs(v) >= 1e-6)
    K2 = x*x + (1-x) / v;
  else
    K2 = x*x - 1/3 - (2/15) * v;

  // Compute SP approximation
  double log_spa = 0.5 * log(0.5 * b / MATH_PI) - 0.5 * log(K2) + b * phi;
  return exp(log_spa);
}


double logcosh(double x) {
  if (x < 0) {
    return -log(2) - x + log1pexp(2*x);
  } else {
    return -log(2) + x + log1pexp(-2*x);
  }
}

double log_cos_rt(double v) {
  double r = sqrt(fabs(v));
  double y;
  if (v >= 0) y = log(cos(r));
  else y = logcosh(r);
  return y;
}
