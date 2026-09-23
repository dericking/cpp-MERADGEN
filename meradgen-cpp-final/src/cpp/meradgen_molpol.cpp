#include "meradgen_molpol.hpp"
#include "meradgen_api.hpp"
#include "meradgen_globals.hpp"

#include <cmath>

namespace meradgen {

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

} // namespace meradgen
