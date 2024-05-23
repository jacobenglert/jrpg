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

#ifndef __RPG__
#define __RPG__

#include "common.h"
#include "rv.h"
#include "rpg_sp_helpers.h"
#include "InvertY.h"
#include "rpg_devroye_helpers.h"


// FCN prototypes
double rpg_devroye(int, double);
double rpg_devroye_1(double, double, double);
double rpg_na(double, double);
double rpg_gamma(double, double);
double rpg_sp(double, double, double);

#endif
