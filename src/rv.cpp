
#include "rv.h"
#include "common.h"


// Generate exponential random variates
double exprnd(double mu) {
  return R::rexp(mu);
  // return -mu * std::log(1.0 - R::runif(0.0,1.0));
}

// Generate inverse Gaussian (mu, 1) random variate
// Michael et al. (1976) in The American Statistician
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
// Michael et al. (1976) in The American Statistician
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

// Generate left-truncated gamma(1/2, 1/2, pi/2) random variates
// Implements Algorithm 3.2 of
// Ref: Chung, Y.: Simulation of truncated gamma variables
// Korean Journal of Computational & Applied Mathematics, 1998, 5, 601-610
// Borrowed from pgdraw
double rltgamma() { // shape = 1/2, rate = 1/2, lower = PI/2

  // This version is used specifically to generate Gamma(1/2, 1/2) random
  // variates left truncated at PI/2. Why PI/2? When we take the reciprocal of
  // the result of this function, the result will be an Inv-Gamma(1/2, 1/2)
  // random variate that is right truncated at 2/PI. This is needed to generate
  // Inv-Gaussian(mu, lambda) random variates which are also right-truncated at
  // 2/PI, as outlined in Windle 2013 thesis Algorithm 3. 2/PI is the t we
  // want to use for the J*(1,z) sampler (not necessarily the other samplers).

  double t = MATH_PI_2;

  double X, gX;

  while (true) {

    // Generate an Exp(gamma) random variate with rate gamma
    // gamma = 1/2 because PI/2 = lower > rate = 1 / 2
    X = exprnd(2.0) + t; // Add t to left truncate

    // Compute g(x)
    gX = M_SQRT_PI_2 / sqrt(X); // Exponentials cancel out - that's nice

    // Accept/Reject
    if (R::runif(0.0,1.0) <= gX)
      return X;
  }

}

// Generate left-truncated gamma(shape, rate, lower) random variates
// Implements Windle 2014 thesis Algorithm 5
// Modified from BayesLogit
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

  }

  return lower * (X / b);

}


// Generate right-truncated inv-gamma(shape, scale, upper) random variates
// This is needed for generating right-truncated inv-Gaussian random variates.
// Note: not currently using. Instead, using rrtinvchisq() since the inv-gamma
// scale is always 1 / 2 in Windle 2013 thesis Algorithm 3.
double rrtinvgamma(double shape, double scale, double upper) {
  double lower = 1.0 / (upper + 0.0001);
  return 1.0 / rltgamma(shape, scale, lower);
}


// Generate right-truncated inverse-Gaussian (mu, 1) random variates at 2/PI
// Algorithm 4 in the Windle (2013) PhD thesis, page 129
// Used for the Devroye sampler
double rrtinvgauss(double mu){
  double X, u;
  double t = MATH_2_PI;

  // Pick sampler
  if (t < mu) {

    // Sampler based on right-truncated inverse-gamma
    // Algorithm 3 in the Windle (2013) PhD thesis, page 128
    while (true) {
      u = R::runif(0.0, 1.0);

      // Code from pgdraw
      X = 1.0 / rltgamma(); // Draw right-truncated inverse gamma (1/2, 1/2, 2/PI)

      if (log(u) < (-0.5 / mu / mu * X)) {
        break;
      }

    }
  } else {

    // Rejection sampler
    X = t + 1.0;
    while (X >= t) {
      X = rinvgauss(mu);
    }
  }

  return X;
}


// Generate right-truncated inverse-Gaussian (mu, 1) random variates
// Algorithm 4 in the Windle (2013) PhD thesis, page 129
// Used for the Devroye sampler
double rrtinvgauss(double mu, double upper){
  double X, u;

  // Pick sampler
  if (upper < mu) {

    // Sampler based on right-truncated inverse-gamma
    // Algorithm 3 in the Windle (2013) PhD thesis, page 128
    while (true) {
      u = R::runif(0.0, 1.0);

      // Code from BayesLogit
      double E1 = exprnd(1.0); double E2 = exprnd(1.0);
      while ( E1*E1 > 2 * E2 / upper) {
        E1 = exprnd(1.0); E2 = exprnd(1.0);
      }
      X = 1 + E1 * upper;
      X = upper / (X * X);

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

// Generate right-truncated inverse-Gaussian(mu, lambda) random variates
double rrtinvgauss(double mu, double lambda, double upper) {
  double X, u;

  // Pick sampler
  if (upper < mu) {

    // Sampler based on right-truncated inverse-gamma (i.e. chi-square)
    // Windle 2013 thesis algorithm 3
    while (true) {

      u = R::runif(0.0, 1.0);

      // Sample right-truncated Inv-Gamma(lambda / 2, 1 / 2)
      // X = rrtinvgamma(0.5 * lambda, 0.5, upper);

      // Sample right-truncated Inverse chi-square (lambda)
      X = rrtinvchi2(lambda, upper);

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



// Generate right-truncated inverse Chi-squared (scale) random variates
// Borrowed from BayesLogit
double rrtinvchi2(double scale, double upper) {
  double R = upper / scale;
  double E = rltnorm(1.0 / sqrt(R));
  double X = scale / (E*E);
  return X;
}

// Used for rrtinvchi2
double rltnorm(double lower) {
  double rho, ppsl;

  if (lower < 0) { // Accept/Reject Normal
    while (true) {
      ppsl = R::rnorm(0.0, 1.0);
      if (ppsl > lower) return ppsl;
    }
  }
  else { // Accept/Reject Exponential
    double astar = 0.5 * (lower + sqrt(lower*lower + 4));
    while (true) {
      ppsl = exprnd(1.0 / astar) + lower;
      rho  = exp( -0.5 * (ppsl - astar) * (ppsl - astar) );
      if (R::runif(0.0,1.0) < rho) return ppsl;
    }
  }

}



