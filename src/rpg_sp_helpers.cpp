#include <Rcpp.h>
#include "rpg_sp_helpers.h"
using namespace Rcpp;

// // Calculate tangent line to eta
// void tangent_to_eta(double x, double z, double mid, Line& tl) {
//
//   FD phi, delta, eta;
//   // double v;
//
//   // Compute v and update phi in place
//   // v = phi_func(x, z, phi);
//   phi_func(x, z, phi);
//
//   // Compute and update delta in place
//   delta_func(x, mid, delta);
//
//   // Calculate eta
//   eta.val = phi.val - delta.val;
//   eta.der = phi.der - delta.der;
//
//   // Calculate and update slope/intercept for line
//   tl.slope = eta.der;
//   tl.intercept = eta.val - eta.der * x;
//
//   // return v;
// }
//
// // Evaluate the phi function
// void phi_func(double x, double z, FD& phi) {
//
//   double v = v_eval(x);       // 2u from paper
//   double u = 0.5 * v;         // u from paper
//   double t = u + 0.5 * z*z;   // t from paper
//
//   phi.val = log(cosh(fabs(z))) - log(cos_rt(v)) - t * x; // Fact 9.4
//   phi.der = -1.0 * t; // Fact 9.5
//
//   // return v;
//
// }
//
// // Compute cos(sqrt(v))
// // Windle thesis Fact 9 proof
// double cos_rt(double v) {
//   double y   = 0.0;
//   double r   = sqrt(fabs(v));
//   if (v >= 0)
//     y = cos(r);
//   else
//     y = cosh(r);
//   return y;
// }
//
// // Compute delta as per Lemma 14
// void delta_func(double x, double mid, FD& delta) {
//   if (x >= mid) {
//     delta.val = log(x) - log(mid);
//     delta.der = 1.0 / x;
//   }
//   else {
//     delta.val = 0.5 * (1 - 1.0 / x) - 0.5 * (1 - 1.0 / mid);
//     delta.der = 0.5 / (x*x);
//   }
// }
//
//
//
// // Saddlepoint approximation given at the top of pg. 69
// double sp_approx(double x, double b, double z) {
//
//   double v = v_eval(x);
//   double u  = 0.5 * v;
//   double z2 = z * z;
//   double t  = u + 0.5 * z2;
//
//   // Compute phi
//   double phi = log(cosh(z)) - log(cos_rt(v)) - t * x;
//
//   // Compute second derivative
//   double K2  = 0.0;
//   if (fabs(v) >= 1e-6)
//     K2 = x*x + (1-x) / v;
//   else
//     K2 = x*x - 1/3 - (2/15) * v;
//
//   // Compute SP approximation
//   double log_spa = 0.5 * log(0.5 * b / MATH_PI) - 0.5 * log(K2) + b * phi;
//   return exp(log_spa);
// }
