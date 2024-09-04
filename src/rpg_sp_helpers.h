
#ifndef __RPG_SP_HELPERS__
#define __RPG_SP_HELPERS__

#include "common.h"
#include "InvertY.h"

// Structure containing function and derivative
struct FD {
  double val;
  double der;
};

// Structure containing slope and intercept of a line
struct Line {
  double slope;
  double intercept;
};

// Functions for obtaining phi, delta, and line tangent to eta
double y_func(double v);
FD get_phi(double, double);
FD get_delta(double, double);
Line get_eta_tangent(double, double, double);

// Saddlepoint approximation
double sp_approx(double, double, double);

// Other useful functions
double logcosh(double x);
double log_cos_rt(double v);

#endif
