#ifndef MERADGEN_GLOBALS_HPP
#define MERADGEN_GLOBALS_HPP

#include "meradgen_constants.hpp"

namespace meradgen {

// Fixed constants, masses and barn conversion constant
extern const double pi, alfa, m, m2, mu, mu2, tau, tau2, barn;

// Event specific values
extern double En, s, als, coeb, coer, Egmin, t, pl;
extern double vprad[4], phirad[4], weight, xs0_save, sirad_out, sinonr_out;
extern int ich, itest;

// Grid variables
extern double grv[NV + 1];
extern double grt1[NT1 + 1];
extern double grz[NZ + 1];

extern double vmin, vmax, vgen;
extern double t1min, t1max, t1gen;
extern double zmin, zmax, zgen;

extern double az, bz, cz, az1, bz1, cz1, az2, bz2, cz2;

extern double bin[NBIN + 1];
extern double argbin[NBIN + 1];
extern double sigbin[NBIN + 1];
extern double step;

// dcanc() variables
extern double xxs, xxt, xxu;

} // namespace meradgen

#endif
