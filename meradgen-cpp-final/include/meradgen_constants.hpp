#ifndef MERADGEN_CONSTANTS_HPP
#define MERADGEN_CONSTANTS_HPP

namespace meradgen {

// Grid dimensions (from Fortran parameter statements)
constexpr int NV = 60;
constexpr int NT1 = 30;
constexpr int NZ = 60;
constexpr int NBIN = 20;

// Integer powers as multiply trees (gfortran -O2 association).
// std::pow(x, 3|4|5) is a libm call and can round 1 ULP away.
inline double ipow3(double x) { return (x * x) * x; }
inline double ipow4(double x) { const double t = x * x; return t * t; }
inline double ipow5(double x) { const double t2 = x * x; return (x * t2) * t2; }

} // namespace meradgen

#endif
