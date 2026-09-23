#include "meradgen_api.hpp"
#include "meradgen_globals.hpp"
#include "meradgen_constants.hpp"
#include "meradgen_parity_trace.hpp"
#include <cmath>
#include <fstream>
#include <cstdint>

namespace meradgen {

  void grid_init() {
    for (int i = 1;  i <= 30; i++) grv[i] = static_cast<double>(i - 1) / 29.0 / 4.0;
    for (int i = 31; i <= 45; i++) grv[i] = 0.25 + static_cast<double>(i - 30) / (45 - 30) / 4.0;
    for (int i = 46; i <= NV; i++) grv[i] = 0.5 + static_cast<double>(i - 45) / (NV - 45) / 2.0;
    for (int i = 1; i <= 7; i++) {
      grt1[i] = 0.1 * static_cast<double>(i * i) / 49.0 / 2.0;
      grt1[31 - i] = 1.0 + grt1[1] - grt1[i];
    }
    for (int i = 8; i <= 15; i++) {
      grt1[i] = (0.1 + 0.9 * static_cast<double>(i - 7) / 8.0) / 2.0;
      grt1[31 - i] = 1.0 + grt1[1] - grt1[i];
    }
    for (int i = 1; i <= 30; i++) grz[i] = 0.5 * static_cast<double>(i * i) / (30.0 * 30.0);
    for (int i = 1; i <= 30; i++) grz[61 - i] = 1.0 - 0.49 * static_cast<double>(i * i) / (30.0 * 30.0);
  }

  void merad_init(double elab) {
    // FORTRAN merad_init takes REAL*4. Do not round here; parity drivers pass
    // a float-rounded elab (ParityWork.md §6). Production may pass full double.
    En = elab;
    s = 2.0 * (En * m + m2);
    als = s * (s - 4.0 * m2);
    coeb = 4.0 * pi * alfa * alfa * barn * (s - 2.0 * m2) / als;
    coer = alfa * alfa * alfa * barn * (s - 2.0 * m2) / als / pi / 4.0;
    Egmin = En * 1e-2;
    grid_init();
  }

  void zd(double tt, double t1, double v) {
    const double u = v - s - tt + 4.0 * m2;
    az = (v - tt) * (v - tt) - 4.0 * m2 * tt;
    bz = -(v * (2.0 * m2 * (tt + t1) + t1 * (v - tt))) + s * (-tt * tt + t1 * v + tt * (t1 + v));
    cz = (s * (tt - t1) + t1 * v) * (s * (tt - t1) + t1 * v) - 4.0 * m2 * (s * (tt - t1) * (tt - t1) + t1 * v * v);
    az1 = az;
    bz1 = -(tt * (s + tt - 4.0 * m2) * (tt - t1)) + (tt * (2 * tt - t1) - 2.0 * m2 * (tt + t1) + s * (tt + t1)) * v - tt * v * v;
    cz1 = ((s + tt) * (tt - t1) - tt * v) * ((s + tt) * (tt - t1) - tt * v) + 4.0 * m2 * (-((s + tt) * (tt - t1) * (tt - t1)) + (tt - t1) * (tt + t1) * v - t1 * v * v);
    az2 = az;
    bz2 = (4.0 * m2 - s - tt) * tt * (4.0 * m2 - t1) + (6.0 * m2 * tt - s * tt - 2.0 * m2 * t1 + s * t1 - tt * t1) * v + (-4.0 * m2 + s) * v * v;
    cz2 = u * (-4.0 * m2 * (s + t1 - 4.0 * m2) * (s + t1 - 4.0 * m2) + (16.0 * m2 * m2 + t1 * t1 - 4.0 * m2 * (s + 2.0 * t1)) * u)
          - 2.0 * (2.0 * m2 - t1) * (4.0 * m2 - s - t1) * u * v + (s + t1 - 4.0 * m2) * (s + t1 - 4.0 * m2) * v * v;
  }

  double fsirv(double v) {
    int nn = 0;
    return fsir(t, 0.0, v, 0.0, pl, nn, -1);
  }

  double fsirv1(double v) {
    int nn = 0;
    return fsir(t, 0.0, v, 0.0, pl, nn, 2);
  }

  // OLD: void vectrec(const float vpgen[4]) {
  void vectrec(const float vpgen[4], double r4) {
    t = static_cast<double>(vpgen[3] * vpgen[3] - vpgen[0] * vpgen[0] - vpgen[1] * vpgen[1] - vpgen[2] * vpgen[2]);
    
    // phi: FORTRAN evaluates atan2 on REAL*4 vpgen, then uses the result as real*8.
    // Double atan2 is more accurate; restore it in the final tree (TODO-FINAL.md).
    // const double phi = std::atan2(static_cast<double>(vpgen[1]), static_cast<double>(vpgen[0])); // [source double, final port]
    const float phi_f = std::atan2(vpgen[1], vpgen[0]);           // [source MERADGEN Fortran vectrec]
    const double phi  = static_cast<double>(phi_f);

    const double al1 = (s - vgen) * (s - vgen) - 4.0 * m2 * s;
    const double al2 = s + 2.0 * t - vgen - 4.0 * m2;
    const double al3 = -s * t * (s + t - vgen - 4.0 * m2) - m2 * vgen * vgen;
    const double al4 = s * (s - vgen - 4.0 * m2) - (s + vgen) * zgen;
    const double al5 = vgen * zgen * (s - vgen - zgen) - m2 * (vgen + zgen) * (vgen + zgen);
    const double al6 = s * (vgen - zgen) - vgen * (vgen + zgen);
    const double al7 = (s + 2.0 * t1gen - zgen - 4.0 * m2) * al1 - al2 * al4;
    const double al8 = 16.0 * al3 * al5 - al7 * al7;
    const double sls = std::sqrt(als);
    const double sl1 = std::sqrt(al1);
    const double sl3 = std::sqrt(al3);
    double sl8 = std::sqrt(al8);

    // OLD: if (urand(iy) > 0.5f) sl8 = -sl8;
    if (r4 > 0.5) sl8 = -sl8;

    const double sp = std::sin(phi);
    const double cp = std::cos(phi);
    const double denom = 4.0 * al1 * sls * sl3;

    vprad[0] = -static_cast<float>((sls * sl1 * sl8 * sp + (4.0 * al3 * al4 - s * al2 * al7) * cp) / denom);
    vprad[1] =  static_cast<float>((sls * sl1 * sl8 * cp - (4.0 * al3 * al4 - s * al2 * al7) * sp) / denom);
    vprad[2] =  static_cast<float>((als * al1 - s * (al7 + al2 * al4)) / (2.0 * std::sqrt(s) * al1 * sls));
    vprad[3] = -static_cast<float>(zgen / 2.0 / std::sqrt(s));
    phirad[0] = static_cast<float>((sls * sl1 * sl8 * sp + (4.0 * al3 * al6 - s * al2 * al7) * cp) / denom);
    phirad[1] = static_cast<float>((-sls * sl1 * sl8 * cp + (4.0 * al3 * al6 - s * al2 * al7) * sp) / denom);
    phirad[2] = static_cast<float>(std::sqrt(s) * (al7 + al2 * al6) / (2.0 * al1 * sls));
    phirad[3] = static_cast<float>((vgen + zgen) / 2.0 / std::sqrt(s));
  }

  // OLD: void meradgen(double ppl, const float vpgen[4]) {
  void meradgen(double ppl, const float vpgen[4], const double rand4[4]) {
    pl = ppl;
    t = static_cast<double>(vpgen[3] * vpgen[3] - vpgen[0] * vpgen[0] - vpgen[1] * vpgen[1] - vpgen[2] * vpgen[2]);
    const double vmax = (s + t) / 2.0;
    vmin = 2.0 * Egmin * m;
#ifdef MERADGEN_PARITY_TRACE
    ptrc_event();
    ptrc_d("t", t);
    ptrc_d("vmin", vmin);
    ptrc_d("vmax", vmax);
#endif

    static int ikey = 0;
    static double distsiv[NV + 1], distarv[NV + 1];
    static double distsit1[4 * NT1 + 1], distart1[4 * NT1 + 1];
    static double distsiz[NZ + 1], distarz[NZ + 1];
    static double sinonr_save = 0, sirad_save = 0;
    static double t1p[5];

    double sirad = 0.0, sitot = 0.0, sirand = 0.0;

    if (itest == 2) goto L22;
    if (itest == 3) goto L44;

    if (ikey == 0) {
      if (itest == 1) {
        step = (vmax - vmin) / static_cast<double>(NBIN);
        for (int i = 0; i <= NBIN; i++)
          bin[i] = vmin + static_cast<double>(i) * step;
      }
      double vvn = vmin;
      double sin_val = 0.0;
      int nn = 0;
      distsiv[0] = 0.0;
      distarv[0] = vvn;
#ifdef MERADGEN_PARITY_TRACE
      ptrc_row('V', 0, vvn, sin_val, distsiv[0], distarv[0]);
#endif
      xs0_save = sig(t, pl, 0);
      for (int iv = 1; iv <= NV; iv++) {
        const double vvo = vvn;
        const double sio = sin_val;
        vvn = vmin + (vmax - vmin) * grv[iv];
        sin_val = fsir(t, 0.0, vvn, 0.0, pl, nn, 2);
        distsiv[iv] = distsiv[iv - 1] + (sin_val + sio) * (vvn - vvo) / 2.0;
        distarv[iv] = vvn;
#ifdef MERADGEN_PARITY_TRACE
        ptrc_row('V', iv, vvn, sin_val, distsiv[iv], distarv[iv]);
#endif
      }
      if (itest == 1) ikey = 1;
      const double u0 = -s - t;
      xs0_save = sig(t, pl, 0);
      const double xsvr = sig(t, pl, 1);
      const double xsB = xsBt(pl, s, t, u0) + xsBt(pl, s, u0, t);
      const double xsF = xs0_save * dcanc(vmin, s, t, u0);
      double xsadd = 0;
      simpsx(1e-22, vmin, 10000, 1e-3, fsirv, xsadd);
      sinonr_save = xs0_save + xsvr + xsB + xsF + xsadd;
#ifdef MERADGEN_PARITY_TRACE
      ptrc_d("xs0", xs0_save);
      ptrc_d("xsvr", xsvr);
      ptrc_d("xsB", xsB);
      ptrc_d("xsF", xsF);
      ptrc_d("xsadd", xsadd);
      ptrc_d("sinonr", sinonr_save);
#endif
    }

    sirad = distsiv[NV];
    sitot = sirad + sinonr_save;
#ifdef MERADGEN_PARITY_TRACE
    ptrc_d("sirad", sirad);
    ptrc_d("sitot", sitot);
#endif
    if (itest == 0) weight = static_cast<float>(sitot / xs0_save);
    if (itest == 1) weight = static_cast<float>(sirad);

    // OLD: sirand = static_cast<double>(urand(iy)) * sitot;
    sirand = rand4[0] * sitot;
#ifdef MERADGEN_PARITY_TRACE
    ptrc_d("sirand_v", sirand);
#endif
    if (sirand <= sinonr_save) {
      vgen = 0.0;
#ifdef MERADGEN_PARITY_TRACE
      ptrc_d("vgen", vgen);
#endif
      t1gen = t;
      zgen = 0.0;
      ich = 0;
      for (int i = 0; i < 4; i++) {
        vprad[i] = vpgen[i];
        phirad[i] = 0.0f;
      }
      return;
    }
    ich = 1;
    for (int iv = 1; iv <= NV; iv++) {
      if (distsiv[iv] > sirand - sinonr_save) {
        vgen = distarv[iv - 1] + (distarv[iv] - distarv[iv - 1]) * (sirand - sinonr_save - distsiv[iv - 1]) / (distsiv[iv] - distsiv[iv - 1]);
#ifdef MERADGEN_PARITY_TRACE
        ptrc_i("v_bin", iv);
        ptrc_d("vgen", vgen);
#endif
        if (itest == 1) return;
        goto L22;
      }
    }


  L22:
    if (ikey == 0) {
      const double u = 4.0 * m2 + vgen - s - t;
      t1min = (2.0 * m2 * t + vgen * (t - vgen - std::sqrt((t - vgen) * (t - vgen) - 4.0 * m2 * t))) / (2.0 * (m2 + vgen));
      t1max = m2 * t * t / (m2 + vgen) / t1min;
      if (itest == 2) ikey = 1;
      step = (t1max - t1min) / static_cast<double>(NBIN);
      for (int i = 0; i <= NBIN; i++) bin[i] = t1min + static_cast<double>(i) * step;
      const double t1z = (-s * t * (t + u) + 2.0 * m2 * vgen * vgen) / ((s - vgen) * (s - vgen) - 4.0 * m2 * s);
      const double t1z1 = (-u * t * (t + s) + 2.0 * m2 * vgen * vgen) / ((u - vgen) * (u - vgen) - 4.0 * m2 * u);
      const double t1z2 = (-s * vgen * (t + s) - 2.0 * m2 * (2.0 * t * u + (u - 2.0 * (s + vgen)) * vgen)) / ((u - vgen) * (u - vgen) - 4.0 * m2 * u);
      t1p[0] = t1min;
      t1p[1] = t1z;
      t1p[2] = std::min(t1z1, t1z2);
      t1p[3] = std::max(t1z1, t1z2);
      t1p[4] = t1max;
      double tt1n = t1min;
      double sin_val = 0.0;
      distsit1[0] = 0.0;
      distart1[0] = tt1n;
#ifdef MERADGEN_PARITY_TRACE
      ptrc_row('T', 0, tt1n, sin_val, distsit1[0], distart1[0]);
#endif
      int nn_t1 = 0;
      for (int i = 0; i < 4; i++) {
        for (int it1 = 1; it1 <= NT1; it1++) {
          const double tt1o = tt1n;
          const double sio = sin_val;
          tt1n = t1p[i] + (t1p[i + 1] - t1p[i]) * grt1[it1];
          zd(t, tt1n, vgen);
          sin_val = fsir(t, tt1n, vgen, 0.0, pl, nn_t1, 1);
          const int idx = i * NT1 + it1;
          distsit1[idx] = distsit1[idx - 1] + (sin_val + sio) * (tt1n - tt1o) / 2.0;
          distart1[idx] = tt1n;
#ifdef MERADGEN_PARITY_TRACE
          ptrc_row('T', idx, tt1n, sin_val, distsit1[idx], distart1[idx]);
#endif
        }
      }
      sirad_save = distsit1[4 * NT1];
#ifdef MERADGEN_PARITY_TRACE
      ptrc_d("sirad_t1", sirad_save);
#endif
      if (itest == 2) weight = static_cast<float>(sirad_save);
    }

    // OLD: sirand = static_cast<double>(urand(iy)) * sirad_save;
    sirand = rand4[1] * sirad_save;
#ifdef MERADGEN_PARITY_TRACE
    ptrc_d("sirand_t1", sirand);
#endif
    ich = 1;
    for (int it1 = 1; it1 <= 4 * NT1; it1++) {
      if (distsit1[it1] > sirand) {
        t1gen = distart1[it1 - 1] + (distart1[it1] - distart1[it1 - 1]) * (sirand - distsit1[it1 - 1]) / (distsit1[it1] - distsit1[it1 - 1]);
#ifdef MERADGEN_PARITY_TRACE
        ptrc_i("t1_bin", it1);
        ptrc_d("t1gen", t1gen);
#endif
        if (itest == 2) return;
        goto L44;
      }
    }

  L44:
    if (ikey == 0) {
      ich = 1;
      zd(t, t1gen, vgen);
      const double det = bz * bz - az * cz;
      const double zmax_val = (-bz + std::sqrt(det)) / az;
      zmin = cz / az / zmax_val;
      if (itest == 3) ikey = 1;
      step = (zmax_val - zmin) / static_cast<double>(NBIN);
      for (int i = 0; i <= NBIN; i++) bin[i] = zmin + static_cast<double>(i) * step;
      double zzn = zmin;
      double sin_val = 0.0;
      int nn = 0;
      distsiz[0] = 0.0;
      distarz[0] = zzn;
#ifdef MERADGEN_PARITY_TRACE
      ptrc_row('Z', 0, zzn, sin_val, distsiz[0], distarz[0]);
#endif
      for (int iz = 1; iz <= NZ; iz++) {
        const double zzo = zzn;
        const double sio = sin_val;
        zzn = zmin + (zmax_val - zmin) * grz[iz];
        sin_val = fsir(t, t1gen, vgen, zzn, pl, nn, 0);
        distsiz[iz] = distsiz[iz - 1] + (sin_val + sio) * (zzn - zzo) / 2.0;
        distarz[iz] = zzn;
#ifdef MERADGEN_PARITY_TRACE
        ptrc_row('Z', iz, zzn, sin_val, distsiz[iz], distarz[iz]);
#endif
      }
      sirad_save = distsiz[NZ];
#ifdef MERADGEN_PARITY_TRACE
      ptrc_d("sirad_z", sirad_save);
#endif
      if (itest == 3) weight = static_cast<float>(sirad_save);
    }

    // OLD: sirand = static_cast<double>(urand(iy)) * sirad_save;
    sirand = rand4[2] * sirad_save;
#ifdef MERADGEN_PARITY_TRACE
    ptrc_d("sirand_z", sirand);
#endif
    for (int iz = 1; iz <= NZ; iz++) {
      if (distsiz[iz] > sirand) {
        zgen = distarz[iz - 1] + (distarz[iz] - distarz[iz - 1]) * (sirand - distsiz[iz - 1]) / (distsiz[iz] - distsiz[iz - 1]);
#ifdef MERADGEN_PARITY_TRACE
        ptrc_i("z_bin", iz);
        ptrc_d("zgen", zgen);
#endif
        if (itest == 3) return;
        goto L55;
      }
    }

  L55:
    if (itest != 0 && itest != 4) return;
    // OLD: vectrec(vpgen);
    vectrec(vpgen, rand4[3]);
  }

} // namespace meradgen
