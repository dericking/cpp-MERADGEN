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

} // namespace meradgen

#endif
