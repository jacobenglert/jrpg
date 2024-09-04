// (C) Nicholas Polson, James Scott, Jesse Windle, 2012-2019

// This file is part of BayesLogit.

// BayesLogit is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.

// BayesLogit is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with
// BayesLogit.  If not, see <https://www.gnu.org/licenses/>.

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
