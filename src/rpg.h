
#ifndef __RPG__
#define __RPG__

#include "common.h"
#include "rv.h"
#include "rpg_sp_helpers.h"
#include "InvertY.h"
#include "rpg_devroye_helpers.h"


// Hybrid Sampler
arma::vec jrpg(arma::vec b, arma::vec z);

// Devroye Method
double rpg_devroye(int b, double z);
double rpg_devroye_1(double z, double r, double K);

// Normal Approximation
double rpg_na(double b, double z);

// Sum-of-Gammas
double rpg_gamma(double b, double z, int trunc = 200);

// Saddlepoint Approximation
double rpg_sp(double b , double z, int maxiter = 200);

#endif
