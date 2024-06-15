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

#include "rv.h"
#include "common.h"


// Generate exponential random variates
double exprnd(double mu) {
  // return R::rexp(mu);
  return -mu * std::log(1.0 - R::runif(0.0,1.0));
}

// Generate inverse Gaussian (mu, 1) random variate
// Similar to Windle (2013) algorithm 1
double rinvgauss(double mu) {

  // Sample
  double u = R::rnorm(0.0,1.0);
  double v = u * u; // differs in that we square u
  double X = mu * (1.0 + 0.5 * mu * v - 0.5 * sqrt(4.0*mu*v + mu*mu * v*v));

  if (R::runif(0.0,1.0) > mu / (mu + X))
    return mu * mu / X;
  else
    return X;
}

// Generate inverse Gaussian (mu, lambda) random variate
// Similar to Windle (2013) algorithm 1
double rinvgauss(double mu, double lambda) {

  // Sample
  double u = R::rnorm(0.0,1.0);
  double v = u * u; // differs in that we square u
  double X = mu * (1.0 + 0.5 * mu * v / lambda - 0.5 / lambda * sqrt(4.0*mu*lambda*v + mu*mu * v*v));

  if (R::runif(0.0,1.0) > mu / (mu + X))
    return mu * mu / X;
  else
    return X;
}

// Left-truncated gamma(1/2,1/2,pi/2) random variates
// Ref: Chung, Y.: Simulation of truncated gamma variables
// Korean Journal of Computational & Applied Mathematics, 1998, 5, 601-610
// Algorithm 3.2
double rltgamma() { // shape = 1/2, rate = 1/2, lower = PI/2

  // This version is used specifically to generate Gamma(1/2,1/2) random
  // variates left truncated at PI/2. Why PI/2? When we take the reciprocal of
  // the result of this function, the result will be an Inv-Gamma(1/2,1/2)
  // random variate that is right truncated at 2/PI. This is needed to generate
  // Inv-Gaussian(mu,lambda) random variates which are also right-truncated at
  // 2/PI, as outlined in Windle Algorithm 3.

  double t = MATH_PI_2;

  double X, gX;

  while (true) {

    // Generate an Exp(gamma) random variate with rate gamma
    // gamma = 1/2 because PI/2 = lower > rate = 1 / 2
    // X = R::rexp(2.0) + t; // rexp() takes the mean not the rate
    X = exprnd(1.0) * 2.0 + t; // Why adding t? I don't actually know...

    // Compute g(x)
    gX = M_SQRT_PI_2 / sqrt(X); // Exponentials cancel out - that's nice

    // Accept/Reject
    if (R::runif(0.0,1.0) <= gX)
      return X;
  }

}

// Left-truncated gamma(shape, rate, lower) random variates
// Windle Algorithm 5
// This is more general than the above function, which only works for shape <= 1
// This is needed directly for the RHS of the saddlepoint approximation, and
// indirectly for the LHS.
double rltgamma(double shape, double rate, double lower) {
  double a = shape;
  double b = rate * lower;

  if (lower <= 0) {
    return 0;
  }

  if (shape < 1) {
    return 0;
  }

  if (shape == 1) {
    return exprnd(1.0) / rate + lower;
  }


  double d1 = b - a;
  double d2 = a - 1;
  double c0 = 0.5 * (d1 + sqrt(d1*d1 + 4*b)) / b;
  double d3 = 1 - c0;

  double X = 0.0;
  bool accept = false;
  while (!accept) {
    X = exprnd(1.0) / c0 + b;
    double log_rho = d2 * log(X) - X * d3;
    double log_M = d2 * log(d2 / d3) - d2;

    double u = R::runif(0.0, 1.0);
    accept = log(u) <= (log_rho - log_M);

    // if (log(R::runif(0.0,1.0)) <= (log_rho - log_M))
    //   break;
  }

  return lower * (X / b);
}

// Right-truncated inverse-gamma(shape, scale, upper) random variates
// is rrtinvchisquare(2*shape, t) if scale = 1/2 (in Devroye sampler shape = 1/2)
// This is needed for the LHS of the saddlepoint approximation
double rrtinvgamma(double shape, double scale, double upper) {
  double lower = 1.0 / (upper + 0.0001);
  return 1.0 / rltgamma(shape, scale, lower);
}


// Sample truncated inverse-Gaussian random(mu, 1) variates
// Algorithm 4 in the Windle (2013) PhD thesis, page 129
// Used for the Devroye sampler
double rrtinvgauss(double mu, double upper){
  double X, u;

  // Pick sampler
  if (upper < mu) {
    // Sampler based on truncated gamma
    // Algorithm 3 in the Windle (2013) PhD thesis, page 128
    while (true) {
      u = R::runif(0.0, 1.0);
      X = 1.0 / rltgamma(); // Draw right-truncated inverse gamma (1/2, 1/2, 2/PI)

      if (log(u) < (-0.5 / mu / mu * X)) {
        break;
      }
    }
  } else {
    // Rejection sampler
    X = upper + 1.0;
    while (X >= upper) {
      X = rinvgauss(mu);
    }
  }
  return X;
}

// Sample right-truncated inverse-Gaussian(mu, lambda) random variates
double rrtinvgauss(double mu, double lambda, double upper) {
  double X, u;

  if (upper < mu) {
    // Sampler based on truncated gamma
    while (true) {
      u = R::runif(0.0, 1.0);
      // Sample right-truncated Inv-Gamma(1/2,lambda/2,upper)
      // X = rrtinvgamma(0.5, 0.5 * lambda, upper);
      X = rrtinvgamma(0.5 * lambda, 0.5, upper);
      //X = rrtinvchi2(lambda, upper);

      if (log(u) <= (-0.5 * lambda / (mu * mu) * X)) {
        break;
      }
    }
  } else {
    // Rejection sampler
    X = upper + 1.0;
    while (X >= upper) {
      X = rinvgauss(mu, lambda);
    }
  }
  return X;
}

// Evaluate CDF of inverse-Gaussian(mu, lambda)
double pinvgauss(double x, double mu, double lambda) {

  double z = 1.0 / mu;
  double b = std::sqrt(lambda / x) * (x * z - 1);
  double a = std::sqrt(lambda / x) * (x * z + 1) * -1.0;

  double y = R::pnorm(b, 0.0, 1.0, true, false) +
    exp(2 * lambda * z + R::pnorm(a, 0.0, 1.0, true, true));

  return y;
}




// TESTING

double texpon_rate(double left, double rate){
  // return left - log(unif()) / rate;
  return exprnd(1.0/rate) + left;
}

double tnorm(double left) {
  double rho, ppsl;

  if (left < 0) { // Accept/Reject Normal
    while (true) {
      ppsl = R::rnorm(0.0, 1.0);
      if (ppsl > left) return ppsl;
    }
  }
  else { // Accept/Reject Exponential
    // return tnorm_tail(left); // Use Devroye.
    double astar = 0.5 * (left + sqrt(left*left + 4));
    while (true) {
      ppsl = texpon_rate(left, astar);
      rho  = exp( -0.5 * (ppsl - astar) * (ppsl - astar) );
      if (R::runif(0.0,1.0) < rho) return ppsl;
    }
  }

}

double rrtinvchi2(double scale, double trunc) {
  double R = trunc / scale;
  // double X = 0.0;
  // // I need to consider using a different truncated normal sampler.
  // double E1 = r.expon_rate(1.0); double E2 = r.expon_rate(1.0);
  // while ( (E1*E1) > (2 * E2 / R)) {
  //   // printf("E %g %g %g %g\n", E1, E2, E1*E1, 2*E2/R);
  //   E1 = r.expon_rate(1.0); E2 = r.expon_rate(1.0);
  // }
  // // printf("E %g %g \n", E1, E2);
  // X = 1 + E1 * R;
  // X = R / (X * X);
  // X = scale * X;
  double E = tnorm(1/sqrt(R));
  double X = scale / (E*E);
  return X;
}
