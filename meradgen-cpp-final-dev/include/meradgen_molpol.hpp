#ifndef MERADGEN_MOLPOL_HPP
#define MERADGEN_MOLPOL_HPP

#include "meradgen_api.hpp"

namespace meradgen {

// One helicity's output. Four-vectors are (px, py, pz, E) GeV in the CM frame
// (native meradgen order). MolPol's old (E, px, py, pz) packing is to_Epxpypz().
struct MolPolEvent {
  double vprad[4];
  double phirad[4];
  double vgen;
  double t1gen;
  double zgen;
  double weight;
  double xs0;
  double sirad;
  double sinonr;
  int ich;
};

// Polarized pieces at a frozen (v, t1, z, ich) for Approach A reweighting.
// See DESIGN_WEIGHT_AT.md for the likelihood-ratio convention.
struct WeightPieces {
  double xs0;       // Born σ(t, pl, 0)
  double sinonr;    // soft+virtual at this pl (ich==0 only; 0 if ich==1)
  double fsir_z;    // dσ/dt/dv/dt1/dz = fsir(..., ikey=0); 0 if ich==0
  double dens;      // sampling density numerator at pl: sinonr (ich=0) or fsir_z
  double dens_ref;  // same at pl_ref
  double lr;        // dens / dens_ref
  double sitot_ref; // echo of reference sitot = sirad_ref + sinonr_ref
  double weight;    // lr * sitot_ref / xs0  (E_ref[weight] → sitot(pl)/xs0(pl))
};

// Copy (px,py,pz,E) into (E,px,py,pz).
inline void to_Epxpypz(const double pxpypzE[4], double Epxpypz[4]) {
  Epxpypz[0] = pxpypzE[3];
  Epxpypz[1] = pxpypzE[0];
  Epxpypz[2] = pxpypzE[1];
  Epxpypz[3] = pxpypzE[2];
}

inline void to_pxpypzE(const double Epxpypz[4], double pxpypzE[4]) {
  pxpypzE[0] = Epxpypz[1];
  pxpypzE[1] = Epxpypz[2];
  pxpypzE[2] = Epxpypz[3];
  pxpypzE[3] = Epxpypz[0];
}

// True if vprad/phirad are finite. NaN/Inf means vectrec hit a negative
// radicand or a vanishing denominator — reject/reroll in Geant4, do not clamp.
bool event_is_finite(const MolPolEvent& e);

// Fill MolPolEvent from the live globals after meradgen().
void snapshot_event(MolPolEvent& e);

// One meradgen call. vp_pxpypzE is (px,py,pz,E) GeV. Returns false if the
// reconstructed 4-momenta are not finite (caller should reroll rand4).
bool generate(double ppl, const double vp_pxpypzE[4], const double rand4[4],
              MolPolEvent& out);

// +ppl then -ppl with the same rand4 (old MolPolMERADGEN Create() pattern).
bool generate_pair(double ppl, const double vp_pxpypzE[4], const double rand4[4],
                   MolPolEvent& plus, MolPolEvent& minus);

// Build CM vpgen (px,py,pz,E) from lab energy and CM scattering angles (radians).
void vpgen_from_angles(double elab, double thetacm, double phi, double vp_pxpypzE[4]);

// --- Approach A: sample once, reweight helicities at frozen kinematics -------

// Sample (v,t1,z,ich,4-vectors) from the pl_ref CDFs (default unpolarized).
// Freezes radiation; caller then calls weight_at for pl=±1.
bool sample_reference(const double vp_pxpypzE[4], const double rand4[4],
                      MolPolEvent& kin, double pl_ref = 0.0);

// Evaluate polarized density / LR / recommended weight at fixed kinematics.
// Does not redraw CDFs or call vectrec.
//   sitot_ref = kin.sirad + kin.sinonr from sample_reference
//   dens_ref  = kin.sinonr if ich==0; else fsir at pl_ref (recomputed)
// Returns false if dens_ref is non-finite / too small, or outputs non-finite.
bool weight_at(double pl, const double vp_pxpypzE[4],
               double v, double t1, double z, int ich,
               double sitot_ref, double dens_ref, WeightPieces& out);

// Convenience: dens_ref / sitot_ref taken from kin; hard dens_ref recomputed
// at pl_ref via fsir(ikey=0).
bool weight_at(double pl, const double vp_pxpypzE[4], const MolPolEvent& kin,
               double pl_ref, WeightPieces& out);

} // namespace meradgen

#endif
