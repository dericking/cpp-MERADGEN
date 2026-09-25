#include "MolPolMeradgen.hh"

#include "meradgen_api.hpp"
#include "meradgen_globals.hpp"

#include <algorithm>
#include <cmath>

namespace {

void boost_cm_to_lab(double elab, double pxpypzE_cm[4]) {
  using meradgen::m;
  using meradgen::s;
  const double gamma = (elab + m) / std::sqrt(s);
  const double beta = std::sqrt(std::max(0.0, 1.0 - 1.0 / (gamma * gamma)));
  const double px = pxpypzE_cm[0];
  const double py = pxpypzE_cm[1];
  const double pz = pxpypzE_cm[2];
  const double e = pxpypzE_cm[3];
  pxpypzE_cm[0] = px;
  pxpypzE_cm[1] = py;
  pxpypzE_cm[2] = gamma * (pz + beta * e);
  pxpypzE_cm[3] = gamma * (e + beta * pz);
}

void cm_from_event(const meradgen::MolPolEvent& ev, double ecm, double pcm,
                   double k2[4], double p2[4], double gam[4]) {
  p2[0] = ev.vprad[0];
  p2[1] = ev.vprad[1];
  p2[2] = ev.vprad[2] - pcm;
  p2[3] = ev.vprad[3] + ecm;

  for (int i = 0; i < 4; ++i)
    gam[i] = ev.phirad[i];

  k2[0] = -p2[0] - gam[0];
  k2[1] = -p2[1] - gam[1];
  k2[2] = -p2[2] - gam[2];
  k2[3] = 2.0 * ecm - p2[3] - gam[3];
}

} // namespace

void MolPolMeradgen::InitBeam(double elab_GeV) {
  elab_ = elab_GeV;
  meradgen::merad_init(elab_);
}

bool MolPolMeradgen::Generate(double thetacm_rad, double phi_rad,
                              const double rand4[4]) {
  weights_ready_ = false;
  for (int i = 0; i < 4; ++i)
    rand4_[i] = rand4[i];

  meradgen::vpgen_from_angles(elab_, thetacm_rad, phi_rad, vp_);

  if (!meradgen::sample_reference(vp_, rand4_, kin_, kPlRef))
    return false;

  ReconstructLab(kin_, elab_, lab_);
  return FillWeights();
}

bool MolPolMeradgen::FillWeights() {
  meradgen::WeightPieces wp{}, wm{};
  if (!meradgen::weight_at(+1.0, vp_, kin_, kPlRef, wp))
    return false;
  if (!meradgen::weight_at(-1.0, vp_, kin_, kPlRef, wm))
    return false;

  // Shared-track analyzing power uses the density LR (see
  // kinematics_examination/METHOD.md). Do NOT put WeightPieces::weight here:
  // that is LR*sitot_ref/xs0(P) and yields A ~ O(0.02), not Born ~ -7/9.
  weights_.w_plus = wp.lr;
  weights_.w_minus = wm.lr;
  weights_.w0 = 0.5 * (wp.lr + wm.lr);
  weights_ready_ = true;
  return true;
}

void MolPolMeradgen::ReconstructLab(const meradgen::MolPolEvent& ev,
                                   double elab_GeV, LabParticles& out) {
  using meradgen::m2;
  using meradgen::s;
  const double ecm = std::sqrt(s) / 2.0;
  const double pcm = std::sqrt(ecm * ecm - m2);

  cm_from_event(ev, ecm, pcm, out.scattered, out.recoil, out.photon);
  boost_cm_to_lab(elab_GeV, out.scattered);
  boost_cm_to_lab(elab_GeV, out.recoil);
  boost_cm_to_lab(elab_GeV, out.photon);
}
