// Cross-section kernels ported from meradgen10.f:
//   sig, vacpol, L1f, xsBt, dgg1, dgg2, fspen, fspens, dcanc
//
// Paper: Afanasev, Ilyichev, Merenkov, Comput. Phys. Commun. (2007)
//        https://doi.org/10.1016/j.cpc.2006.10.002
//        https://arxiv.org/abs/hep-ph/0603027

#include "meradgen_api.hpp"
#include "meradgen_globals.hpp"
#include "meradgen_constants.hpp"
#include "meradgen_parity_trace.hpp"
#include <cmath>
#include <cstdlib>

namespace meradgen {

  static double L1f(double xt, double xs, double xm2);
  static double vacpol(double tt);
  static double dgg1(double xs, double xt, double xu);
  static double dgg2(double xs, double xt, double xu);
  static double fspens(double x);

  double sig(double tt, double pl_in, int i) {
    const double u = 4 * m2 - s - tt;
    const double ss = s - 2.0 * m2;
    const double u1 = (ss * ss + u * u) / 2.0 + 2.0 * m2 * (s + 2.0 * tt - 3.0 * m2);
    const double u2 = (ss * ss + tt * tt) / 2.0 + 2.0 * m2 * (s + 2.0 * u - 3.0 * m2);
    const double u3 = (s - 2.0 * m2) * (s - 6.0 * m2);
    const double pl1 = -tt * (-tt * s * s / 2.0 / als - ss);
    const double pl2 = -tt * (-tt * s * s / 2.0 / als - 2.0 * m2) + als / 2.0 - ss * s;
    const double pl3 = -(ss * ss * ss * ss + 4.0 * m2 * (ss * ss * tt + m2 * (-ss * ss + 2.0 * tt * tt - 4.0 * m2 * tt))) / als;

    if (i == 0) return coeb * (u1 / (tt * tt) + u2 / (u * u) + u3 / u / tt + pl_in * (pl1 / (tt * tt) + pl2 / (u * u) + pl3 / u / tt));

    const double dsvt = vacpol(-tt) + L1f(tt, s, m2);
    const double dsvu = vacpol(-u) + L1f(u, s, m2);
    const double du1 = dsvt, du2 = dsvu, du3 = dsvt + dsvu;
    const double dp1 = dsvt, dp2 = dsvu, dp3 = dsvt + dsvu;
    const double out = alfa / pi * coeb * (du1 * u1 / (tt * tt) + du2 * u2 / (u * u) + du3 * u3 / u / tt / 2.0 + pl_in * (dp1 * pl1 / (tt * tt) + dp2 * pl2 / (u * u) + dp3 * pl3 / u / tt / 2.0));
    ptrc_d("sig1.t", tt);
    ptrc_d("sig1.u", u);
    ptrc_d("sig1.dsvt", dsvt);
    ptrc_d("sig1.dsvu", dsvu);
    ptrc_d("sig1.out", out);
    return out;
  }


  double vacpol(double tt) {
    // Fortran vacpol DATA am2 (lepton mass^2). Not the same as m2 / mu2 / tau2.
    const double am2[3] = { 0.26110e-6, 0.111637e-1, 3.18301 }; // [source MERADGEN Fortran]
    // const double am2[3] = { m2, mu2, tau2 };                 // [source PDG, via globals.cpp]
    ptrc_d("vp.in", tt);
    double suml = 0.0;
    for (int i = 0; i < 3; i++) {
      const double a2 = 2.0 * am2[i];
      const double sqlmi = std::sqrt(tt * tt + 2.0 * a2 * tt);
      const double allmi = std::log((sqlmi + tt) / (sqlmi - tt)) / sqlmi;
      ptrc_d("vp.a2", a2);
      ptrc_d("vp.sqlmi", sqlmi);
      ptrc_d("vp.allmi", allmi);
      // FORTRAN default REAL*4 literals (2., 10./9., 3., 4., 1.) then promote.
      // 10./9. is the inexact one; do not use 10.0/9.0 here (ParityWork Stage 4).
      // TODO-FINAL.md: restore double literals after parity.
      suml += static_cast<double>(2.f) * (tt + a2) * allmi / static_cast<double>(3.f)
           - static_cast<double>(10.f / 9.f)
           + static_cast<double>(4.f) * a2
                 * (static_cast<double>(1.f) - a2 * allmi)
                 / static_cast<double>(3.f) / tt;
    }
    double aaa, bbb, ccc;
    if (tt < 1.0) {
      aaa = -1.345e-9;
      bbb = -2.302e-3;
      ccc = static_cast<double>(4.091f);  // FORTRAN `ccc = 4.091` (REAL*4 literal)
    } else if (tt < 64.0) {
      aaa = -1.512e-3;
      bbb = -2.822e-3;
      ccc = static_cast<double>(1.218f);
    } else {
      aaa = -1.1344e-3;
      bbb = -3.0680e-3;
      ccc = static_cast<double>(9.9992e-1f);
    }
    const double sumh = -(aaa + bbb * std::log(1.0 + ccc * tt)) * 2.0 * pi / alfa;
    const double out = suml + sumh;
    ptrc_d("vp.suml", suml);
    ptrc_d("vp.sumh", sumh);
    ptrc_d("vp.out", out);
    return out;
  }


  static double L1f(double xt, double xs, double xm2) {
    const double a = std::abs(xt);
    const double out = -2.0 * std::log(a / xs) * (std::log(a / xm2) - 1.0)
        + std::log(a / xm2) + std::log(a / xm2) * std::log(a / xm2)
        + 4.0 * (pi * pi / 12.0 - 1.0);
    ptrc_d("L1.xt", xt);
    ptrc_d("L1.xs", xs);
    ptrc_d("L1.xm2", xm2);
    ptrc_d("L1.absxt", a);
    ptrc_d("L1.log_xs", std::log(a / xs));
    ptrc_d("L1.log_xm2", std::log(a / xm2));
    ptrc_d("L1.log_xm2_sq", std::log(a / xm2) * std::log(a / xm2));
    ptrc_d("L1.pi2term", 4.0 * (pi * pi / 12.0 - 1.0));
    ptrc_d("L1.out", out);
    return out;
  }


  double xsBt(double pl_in, double xs, double xt, double xu) {  // called from meradgen_main
    return 2.0 * alfa * alfa * alfa / (xt * xt) * barn
        * ((1.0 + pl_in) * xu * xu / xs * dgg1(xs, xt, xu)
            - (1.0 - pl_in) * xs * xs / xu * dgg2(xs, xt, xu));
  }


  static double dgg1(double xs, double xt, double xu) {
    const double LS = std::log(xs / std::abs(xt));
    const double LX = std::log(xu / xt);
    const double dgg = LS * LS * (xs * xs + xu * xu) / 2.0 / xt - LS * xu - (LX * LX + pi * pi) * xu * xu / xt;
    return 2.0 * std::log(xs / std::abs(xu)) * std::log(std::sqrt(std::abs(xu / xs))) - xt / (xu * xu) * dgg;
  }


  static double dgg2(double xs, double xt, double xu) {
    const double LS = std::log(xs / std::abs(xt));
    const double LX = std::log(xu / xt);
    const double dgg = LS * LS * xs * xs / xt + LX * xs - (LX * LX + pi * pi) * (xs * xs + xu * xu) / 2.0 / xt;
    return 2.0 * std::log(xs / std::abs(xu)) * std::log(std::sqrt(std::abs(xu / xs))) - xt / (xs * xs) * dgg;
  }


  static double fspen(double x) {
    const double f1 = 1.644934;
    // Boundaries match FORTRAN arithmetic IF: if(x-1d0)4,4,5 and if(x-2d0)6,6,7.
    // x==1 must take the f1 branch (C++ `x < 1` sent it into 0*(-inf) → NaN).
    if (x < -1.0) return -0.5 * std::log(1.0 - x) * std::log(x * x / (1.0 - x)) - f1 + fspens(1.0 / (1.0 - x));
    if (x < 0)    return -0.5 * std::log(1.0 - x) * std::log(1.0 - x) - fspens(x / (x - 1.0));
    if (x <= 0.5) return fspens(x);
    if (x <= 1.0) return f1 - std::log(x) * std::log(1.0 - x + 1e-10) - fspens(1.0 - x);
    if (x <= 2.0) return f1 - 0.5 * std::log(x) * std::log((x - 1.0) * (x - 1.0) / x) + fspens(1.0 - 1.0 / x);
    return 2.0 * f1 - 0.5 * std::log(x) * std::log(x) - fspens(1.0 / x);
  }


  static double fspens(double x) {
    double f = 0.0, a = 1.0, an = 0.0;
    const double tch = 1e-16;
    for (;;) {
      an += 1.0;
      a *= x;
      const double b = a / (an * an);
      f += b;
      if (b <= tch) break;
    }
    return f;
  }


  double dcanc(double xvmin, double xs, double xt, double xu) {
    xxs = xs;
    xxt = xt;
    xxu = xu;
    const double lm = std::log(-xt / m2);
    const double lr = std::log(-xu / xs);
    const double del1s = -2.5 * lm * lm + (3.0 - 2.0 * lr) * lm - lr * lr / 2.0
        - (lm - 1.0) * std::log(xs * (xs + xt) / (xt * xt)) - pi * pi / 3.0 + 1.0;
    const double del1h = -lm * lm / 2.0
        + (std::log(xt * xt * (xs + xt) * (xs + xt) * (xs - xvmin) / xs / (xvmin - xt) / xvmin / (xs + xt - xvmin) / (xs + xt - xvmin)) + 1.0) * lm
        - std::log(-xvmin / xt) * std::log(-xvmin / xt) / 2.0 - std::log(1.0 - xvmin / xt) * std::log(1.0 - xvmin / xt)
        + std::log((xs + xt) / (xs + xt - xvmin)) * std::log((xs + xt) * (xs + xt - xvmin) / (xt * xt))
        + std::log((xvmin - xs) / xt) * std::log(1.0 - xvmin / xs) + std::log(-xvmin / xt)
        + fspen((xs - xvmin) / xs) - fspen((xt - xvmin) / xt)
        + 2.0 * (fspen(xvmin / xs) - fspen(xvmin / xt) - fspen(xvmin / (xt + xs)))
        - pi * pi / 6.0;
    return alfa / pi * (4.0 * std::log(xvmin / std::sqrt(m2 * xs)) * (std::log(xt * xu / m2 / xs) - 1.0) + del1s + del1h);
  }


} // namespace meradgen
