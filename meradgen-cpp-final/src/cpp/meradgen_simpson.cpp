#include "meradgen_api.hpp"
#include "meradgen_globals.hpp"
#include <cmath>
#include <cstdlib>

namespace meradgen {

static void simps_impl(double a1, double b1, double h1, double reps1, double aeps1,
                      double (*funct)(double), double& x, double& ai, double& aih, double& aiabs);

void simps(double a1, double b1, double h1, double reps1, double aeps1,
           double (*funct)(double), double x, double& ai, double& aih, double& aiabs) {
  double xx = x;
  simps_impl(a1, b1, h1, reps1, aeps1, funct, xx, ai, aih, aiabs);
}

static void simps_impl(double a1, double b1, double h1, double reps1, double aeps1,
                       double (*funct)(double), double& x, double& ai, double& aih, double& aiabs) {
  double f[7], p[5];
  double h = (b1 > a1) ? h1 : -h1;
  double s = (h > 0) ? 1.0 : -1.0;
  double a = a1, b = b1;
  ai = 0.0;
  aih = 0.0;
  aiabs = 0.0;
  p[1] = 4.0;
  p[3] = 4.0;
  p[2] = 2.0;
  p[4] = 1.0;

  if (b <= a) return;

  double reps = std::abs(reps1);
  double aeps = std::abs(aeps1);
  for (int k = 0; k < 7; k++) f[k] = 1e16;
  x = a;
  double c = 0.0;
  f[0] = funct(x) / 3.0;

  for (;;) {
    double x0 = x;
    if ((x0 + 4.0 * h - b) * s > 0) {
      h = (b - x0) / 4.0;
      if (h == 0) return;
      for (int k = 1; k < 7; k++) f[k] = 1e16;
      c = 1.0;
    }
    double di2 = f[0];
    double di3 = std::abs(f[0]);
    for (int k = 1; k <= 4; k++) {
      x += h;
      if ((x - b) * s >= 0) x = b;
      if (f[k] >= 1e16) f[k] = funct(x) / 3.0;
      di2 += p[k] * f[k];
      di3 += p[k] * std::abs(f[k]);
    }
    double di1 = (f[0] + 4.0 * f[2] + f[4]) * 2.0 * h;
    di2 *= h;
    di3 *= h;
    double eps = (reps > 0) ? std::abs((aiabs + di3) * reps) : aeps;
    if (eps < aeps) eps = aeps;
    double delta = std::abs(di2 - di1);
    if (delta <= eps) {
      if (delta <= eps / 8.0) {
        h *= 2.0;
        f[0] = f[4];
        f[1] = f[5];
        f[2] = f[6];
        for (int k = 3; k < 7; k++) f[k] = 1e16;
      } else {
        f[0] = f[4];
        f[2] = f[5];
        f[4] = f[6];
        f[1] = f[3] = f[5] = f[6] = 1e16;
      }
      di1 = di2 + (di2 - di1) / 15.0;
      ai += di1;
      aih += di2;
      aiabs += di3;
      if (c == 0) continue;
      return;
    }
    h /= 2.0;
    f[6] = f[4];
    f[5] = f[3];
    f[4] = f[2];
    f[2] = f[1];
    f[1] = f[3] = 1e16;
    x = x0;
    c = 0.0;
  }
}

void simpsx(double a, double b, int np, double ep, double (*func)(double), double& res) {
  double step_val = (b - a) / np;
  double ra = 0, r2 = 0, r3 = 0;
  simps(a, b, step_val, ep, 1e-18, func, ra, res, r2, r3);
}

// FORTRAN meradgen10.f also ships simpt/simpu and simptx/simpux. Their bodies
// are copies of simps/simpsx (verified 2026-08-31). Nothing in the generator
// calls them; they are here so the port is complete. A later pass may drop them.
void simpt(double a1, double b1, double h1, double reps1, double aeps1,
           double (*funct)(double), double x, double& ai, double& aih, double& aiabs) {
  simps(a1, b1, h1, reps1, aeps1, funct, x, ai, aih, aiabs);
}

void simpu(double a1, double b1, double h1, double reps1, double aeps1,
           double (*funct)(double), double x, double& ai, double& aih, double& aiabs) {
  simps(a1, b1, h1, reps1, aeps1, funct, x, ai, aih, aiabs);
}

void simptx(double a, double b, int np, double ep, double (*func)(double), double& res) {
  simpsx(a, b, np, ep, func, res);
}

void simpux(double a, double b, int np, double ep, double (*func)(double), double& res) {
  simpsx(a, b, np, ep, func, res);
}

} // namespace meradgen
