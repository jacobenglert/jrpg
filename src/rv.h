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

#ifndef __RV__
#define __RV__

// Exponential
double exprnd(double mu);

// Left-Truncated Gamma
double rltgamma(); // shape = 1/2, rate = 1/2, lower = MATH_PI_2
double rltgamma(double shape, double rate, double lower);

// Right-Truncated Inverse-Gamma
double rrtinvgamma(double shape, double scale, double upper);

// Inverse Gaussian
double rinvgauss(double mu); // lambda = 1
double rinvgauss(double mu, double lambda);
double pinvgauss(double x, double mu, double lambda);

// Right-Truncated Inverse Gaussian
double rrtinvgauss(double mu, double upper); // lambda = 1
double rrtinvgauss(double mu, double lambda, double upper);

// Right-Truncated Inverse Gamma
double rrtinvgamma(double shape, double scale, double upper);



// TESTING
double texpon_rate(double left, double rate);
double tnorm(double left);
double rrtinvchi2(double scale, double trunc);

#endif
