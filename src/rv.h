
#ifndef __RV__
#define __RV__

// Exponential
double exprnd(double mu);

// Left-Truncated Gamma
double rltgamma(); // shape = 1/2, rate = 1/2, lower = MATH_PI_2
double rltgamma(double shape, double rate, double lower);

// Right-Truncated Inverse Gamma
double rrtinvgamma(double shape, double scale, double upper);

// Right-Truncated Inverse Chi-Squared
double rrtinvchi2(double scale, double upper);

// Left-Truncated Gaussian (standard)
double rltnorm(double lower);

// Inverse Gaussian
double rinvgauss(double mu); // lambda = 1
double rinvgauss(double mu, double lambda);
double pinvgauss(double x, double mu, double lambda);

// Right-Truncated Inverse Gaussian
double rrtinvgauss(double mu); // lambda = 1, upper = 2/PI
double rrtinvgauss(double mu, double upper); // lambda = 1
double rrtinvgauss(double mu, double lambda, double upper);


#endif
