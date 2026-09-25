#ifndef __MOLPOLEVENT_HH
#define __MOLPOLEVENT_HH

/*!
  Event information class.  This needs to
  contain all the information we need to
  generate particles and provide sufficient
  output.
*/

#include <vector>
#include "G4ThreeVector.hh"

class G4ParticleDefinition;

class MolPolEvent {
public:
  MolPolEvent();
  ~MolPolEvent();

  void ProduceNewParticle( G4ThreeVector, G4ThreeVector, G4String );
  void SetEffCrossSection( G4double xs ){ fEffXs = xs; }
  void SetAsymmetry( G4double Asym ){ fAsym = Asym; }
  void SetThCoM( G4double th ){fThCoM = th;}
  void SetPhCoM( G4double ph ){fPhCoM = ph;}
  void SetTargetMomentum( G4double pt ){fTargMom = pt;}

  void Reset();
  void UndoLastParticle();

  G4bool EventIsSane();
  void   Print();

private:

public:
  // Particles to be produced
  std::vector<G4ThreeVector>    fPartPos;
  std::vector<G4ThreeVector>    fPartMom;  // Generated direction (no ms)
  std::vector<G4ParticleDefinition *> fPartType;

  G4double fThCoM;
  G4double fPhCoM;
  G4double fEffXs;
  G4double fAsym;
  G4double fUnpolWght;
  G4double fpolPlusWghtX;
  G4double fpolPlusWghtY;
  G4double fpolPlusWghtZ;
  G4double fpolMinusWghtX;
  G4double fpolMinusWghtY;
  G4double fpolMinusWghtZ;
  G4double fTargMom;

  // MERADGEN Approach A diagnostics (filled when gentype == "meradgen")
  G4bool   fMeradUsed;
  G4bool   fMeradWeightsReady;
  G4double fMeradW0;
  G4double fMeradWplus;
  G4double fMeradWminus;
  G4int    fMeradIch;
  G4double fMeradV;
  G4double fMeradT1;
  G4double fMeradZ;
  G4double fMeradE1;   // lab energy of scattered e- (GeV)
  G4double fMeradE2;   // lab energy of recoil e- (GeV)
  G4double fMeradRand4[4];
};

#endif//__MOLPOLEVENT_HH
