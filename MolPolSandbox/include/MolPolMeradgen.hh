#ifndef MOLPOL_MERADGEN_HH
#define MOLPOL_MERADGEN_HH

#include "meradgen_molpol.hpp"

// Geant4-facing adapter around meradgen-cpp-final-2 (Approach A).
// Units inside: GeV, radians, four-vectors (px, py, pz, E).
// See MERADGEN_INTEGRATION.md.

class MolPolMeradgen {
public:
  struct LabParticles {
    double scattered[4]; // beam-side k2, (px,py,pz,E) GeV lab
    double recoil[4];    // target-side p2
    double photon[4];    // real photon (zero if non-radiative); not tracked
  };

  struct HelicityWeights {
    double w0;      // unpolarized / rate weight (mean of ± or pl=0 sample)
    double w_plus;  // pl = +1
    double w_minus; // pl = -1
  };

  void InitBeam(double elab_GeV);

  // Approach A: sample radiation at pl_ref=0, freeze lab kinematics,
  // then weight_at(+1) and weight_at(-1). Returns false if kinematics are
  // non-finite or either weight_at fails (caller should reroll rand4).
  bool Generate(double thetacm_rad, double phi_rad, const double rand4[4]);

  double Elab() const { return elab_; }

  const LabParticles& Lab() const { return lab_; }
  const meradgen::MolPolEvent& Kinematics() const { return kin_; }
  const HelicityWeights& Weights() const { return weights_; }
  const double* Rand4() const { return rand4_; }
  const double* Vpgen() const { return vp_; }

  bool WeightsReady() const { return weights_ready_; }

  static void ReconstructLab(const meradgen::MolPolEvent& ev, double elab_GeV,
                             LabParticles& out);

private:
  bool FillWeights();

  double elab_ = 0.0;
  double rand4_[4] = {0, 0, 0, 0};
  double vp_[4] = {0, 0, 0, 0};
  meradgen::MolPolEvent kin_{};
  LabParticles lab_{};
  HelicityWeights weights_{};
  bool weights_ready_ = false;
  static constexpr double kPlRef = 0.0;
};

#endif
