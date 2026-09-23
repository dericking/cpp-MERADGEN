#ifndef MERADGEN_API_HPP
#define MERADGEN_API_HPP

#include <cstdint>

namespace meradgen {

// RNG (Fortran URAND - updates seed in place, returns [0,1))
// OLD: float urand(int32_t& iy);

// Initialization
void grid_init();
void merad_init(double elab);

// Main event generation: ppl = polarization, vpgen = 4-momentum of virtual photon (input)
// Outputs in globals: vprad, phirad, weight, ich, xs0_save, vgen, t1gen, zgen
// OLD: void meradgen(double ppl, const float vpgen[4]);
void meradgen(double ppl, const float vpgen[4], const double rand4[4]);
// OLD: void vectrec(const float vpgen[4]);
void vectrec(const float vpgen[4], double r4);

// Kinematics helper (sets gr.inc az, bz, cz, ...)
void zd(double tt, double t1, double v);

// Cross sections and physics (used internally; exposed for tests)
double sig(double tt, double pl_in, int i);
double xsBt(double pl_in, double xs, double xt, double xu);
double dcanc(double xvmin, double xs, double xt, double xu);
double fsir(double tt, double t1, double v, double z, double pl_in, int& nn, int ikey);

// Integrand wrappers for simpsx (take single double, return double)
double fsirv(double v);
double fsirv1(double v);

// Integration
void simps(double a1, double b1, double h1, double reps1, double aeps1,
           double (*funct)(double), double x, double& ai, double& aih, double& aiabs);
void simpsx(double a, double b, int np, double ep, double (*func)(double), double& res);
// Copies of simps/simpsx in the FORTRAN sources; unused by meradgen itself.
void simpt(double a1, double b1, double h1, double reps1, double aeps1,
           double (*funct)(double), double x, double& ai, double& aih, double& aiabs);
void simpu(double a1, double b1, double h1, double reps1, double aeps1,
           double (*funct)(double), double x, double& ai, double& aih, double& aiabs);
void simptx(double a, double b, int np, double ep, double (*func)(double), double& res);
void simpux(double a, double b, int np, double ep, double (*func)(double), double& res);

} // namespace meradgen

#endif
