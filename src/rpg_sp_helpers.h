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

// Function and derivative.
struct FD {
  double val;
  double der;
};

struct Line {
  double slope;
  double intercept;
};

// double tangent_to_eta(double x, double z, double mid, Line t1);
// double phi_func(double x, double z, FD phi);
// double delta_func(double x, double mid, FD delta);
// double cos_rt(double v);
// double sp_approx(double, double, double);

#endif
