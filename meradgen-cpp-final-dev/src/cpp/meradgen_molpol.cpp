#include "meradgen_molpol.hpp"
#include "meradgen_api.hpp"
#include "meradgen_globals.hpp"

#include <cmath>
#include <limits>

namespace meradgen {

namespace {

constexpr double kDensRefEps = 1e-300;

double mandelstam_t(const double vp[4]) {
  return vp[3] * vp[3] - vp[0] * vp[0] - vp[1] * vp[1] - vp[2] * vp[2];
}

// Soft+virtual integrated piece (sinonr). Mutates globals t, pl, vmin, xs0_save.
// Caller must have called merad_init (s, coer, Egmin, …).
double compute_sinonr(double pl_in, double t_in, double& xs0_out) {
  pl = pl_in;
  t = t_in;
  vmin = 2.0 * Egmin * m;
  xs0_out = sig(t_in, pl_in, 0);
  xs0_save = xs0_out;
  const double u0 = -s - t_in;
  const double xsvr = sig(t_in, pl_in, 1);
  const double xsB = xsBt(pl_in, s, t_in, u0) + xsBt(pl_in, s, u0, t_in);
  const double xsF = xs0_out * dcanc(vmin, s, t_in, u0);
  double xsadd = 0.0;
  simpsx(1e-22, vmin, 10000, 1e-3, fsirv, xsadd);
  return xs0_out + xsvr + xsB + xsF + xsadd;
}

// Differential hard density fsir(..., ikey=0). Sets zd globals and t, pl.
double compute_fsir_z(double pl_in, double t_in, double t1, double v, double z) {
  t = t_in;
  pl = pl_in;
  zd(t_in, t1, v);
  int nn = 0;
  return fsir(t_in, t1, v, z, pl_in, nn, 0);
}

bool pieces_finite(const WeightPieces& w) {
  return std::isfinite(w.xs0) && std::isfinite(w.sinonr) && std::isfinite(w.fsir_z)
      && std::isfinite(w.dens) && std::isfinite(w.dens_ref) && std::isfinite(w.lr)
      && std::isfinite(w.weight) && std::isfinite(w.sitot_ref);
}

} // namespace

bool event_is_finite(const MolPolEvent& e) {
  for (int i = 0; i < 4; i++) {
    if (!std::isfinite(e.vprad[i]) || !std::isfinite(e.phirad[i]))
      return false;
  }
  return std::isfinite(e.vgen) && std::isfinite(e.t1gen) && std::isfinite(e.zgen)
      && std::isfinite(e.weight);
}

void snapshot_event(MolPolEvent& e) {
  for (int i = 0; i < 4; i++) {
    e.vprad[i] = vprad[i];
    e.phirad[i] = phirad[i];
  }
  e.vgen = vgen;
  e.t1gen = t1gen;
  e.zgen = zgen;
  e.weight = weight;
  e.xs0 = xs0_save;
  e.sirad = sirad_out;
  e.sinonr = sinonr_out;
  e.ich = ich;
}

bool generate(double ppl, const double vp_pxpypzE[4], const double rand4[4],
              MolPolEvent& out) {
  itest = 0;
  meradgen::meradgen(ppl, vp_pxpypzE, rand4);
  snapshot_event(out);
  return event_is_finite(out);
}

bool generate_pair(double ppl, const double vp_pxpypzE[4], const double rand4[4],
                   MolPolEvent& plus, MolPolEvent& minus) {
  const bool a = generate(ppl, vp_pxpypzE, rand4, plus);
  const bool b = generate(-ppl, vp_pxpypzE, rand4, minus);
  return a && b;
}

void vpgen_from_angles(double elab, double thetacm, double phi, double vp_pxpypzE[4]) {
  merad_init(elab);
  const double ecm = std::sqrt(s) / 2.0;
  const double pcm = std::sqrt(ecm * ecm - m2);
  vp_pxpypzE[0] = -pcm * std::sin(thetacm) * std::cos(phi);
  vp_pxpypzE[1] = -pcm * std::sin(thetacm) * std::sin(phi);
  vp_pxpypzE[2] = pcm * (1.0 - std::cos(thetacm));
  vp_pxpypzE[3] = 0.0;
}

bool sample_reference(const double vp_pxpypzE[4], const double rand4[4],
                      MolPolEvent& kin, double pl_ref) {
  return generate(pl_ref, vp_pxpypzE, rand4, kin);
}

bool weight_at(double pl_in, const double vp_pxpypzE[4],
               double v, double t1, double z, int ich_in,
               double sitot_ref, double dens_ref, WeightPieces& out) {
  out = WeightPieces{};
  out.sitot_ref = sitot_ref;
  out.dens_ref = dens_ref;

  if (!std::isfinite(sitot_ref) || !std::isfinite(dens_ref)
      || std::fabs(dens_ref) < kDensRefEps) {
    out.lr = std::numeric_limits<double>::quiet_NaN();
    out.weight = std::numeric_limits<double>::quiet_NaN();
    return false;
  }

  const double t_in = mandelstam_t(vp_pxpypzE);

  if (ich_in == 0) {
    // Soft: dens = sinonr (needs soft Simpson). xs0 from the same path.
    double xs0_pl = 0.0;
    out.sinonr = compute_sinonr(pl_in, t_in, xs0_pl);
    out.xs0 = xs0_pl;
    out.fsir_z = 0.0;
    out.dens = out.sinonr;
  } else {
    // Hard: dens = fsir(ikey=0). Born only — skip soft integral.
    out.xs0 = sig(t_in, pl_in, 0);
    xs0_save = out.xs0;
    out.sinonr = 0.0;  // not computed for hard; see DESIGN_WEIGHT_AT.md
    out.fsir_z = compute_fsir_z(pl_in, t_in, t1, v, z);
    out.dens = out.fsir_z;
  }

  out.lr = out.dens / dens_ref;
  out.weight = out.lr * sitot_ref / out.xs0;

  if (!pieces_finite(out) || std::fabs(out.xs0) < kDensRefEps)
    return false;
  return true;
}

bool weight_at(double pl_in, const double vp_pxpypzE[4], const MolPolEvent& kin,
               double pl_ref, WeightPieces& out) {
  const double sitot_ref = kin.sirad + kin.sinonr;
  double dens_ref = 0.0;

  if (kin.ich == 0) {
    dens_ref = kin.sinonr;
  } else {
    const double t_in = mandelstam_t(vp_pxpypzE);
    dens_ref = compute_fsir_z(pl_ref, t_in, kin.t1gen, kin.vgen, kin.zgen);
  }

  return weight_at(pl_in, vp_pxpypzE, kin.vgen, kin.t1gen, kin.zgen, kin.ich,
                   sitot_ref, dens_ref, out);
}

} // namespace meradgen
