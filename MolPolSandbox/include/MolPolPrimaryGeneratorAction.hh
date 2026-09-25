#ifndef MolPolPrimaryGeneratorAction_h
#define MolPolPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include "MolPolIO.hh"
#include "G4ThreeVector.hh"
#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <vector>

class MolPolDetectorConstruction;
class G4ParticleGun;
class G4Event;
class PrimaryGeneratorMessenger;
class Simulation;
class remollMultScatt;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class MolPolPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  // Shared moller vertex / beam sample (MS, external eloss, Levchuk, angles).
  struct MollerVertex {
    G4double beamE;
    G4double thcom;
    G4double phcom;
    G4double xpos;
    G4double ypos;
    G4double zpos;
    G4ThreeVector direction;
    G4double zLum;
  };

  MolPolPrimaryGeneratorAction();
  virtual ~MolPolPrimaryGeneratorAction();

  void GeneratePrimaries(G4Event*);
  void SetAngle(double deg){
    double rad = deg*M_PI/180.;
    angle = deg;
    particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,sin(rad),cos(rad)));
  }
  void SetIO( MolPolIO *io ){ fIO = io; }
  void SetDetector( MolPolDetectorConstruction *det ){ fDet = det; }
  void rand();
  double xPos(){return particleGun->GetParticlePosition().x();}
  double yPos(){return particleGun->GetParticlePosition().y();}
  double GetAngle(){return angle;}
  G4double fXmin, fXmax, fYmin, fYmax;
  G4double fXsmear, fYsmear;
  G4double fX, fY, fZ;
  G4double fTargLen;
  G4double fTargPol;//[0,1]
  G4double fBeamE;
  G4double fEmin, fEmax;
  G4double fthetaMin, fthetaMax;
  G4double fthetaComMin, fthetaComMax;
  G4double fphiMin, fphiMax;
  G4String fBeamPol;
  void SourceModeSet(G4int );
  void SetGenerator(G4String genname){ gentype = genname; }
  G4String GetGenerator() const { return gentype; }

  // Radiative-correction model for gentype "moller":
  //   0 = none (Born)
  //   1 = Alexander / Peking structure-function approximation
  //   2 = MERADGEN Approach A (sample P=0, freeze pair, reweight ±)
  void SetRCtype(G4int t);
  G4int GetRCtype() const { return fRCtype; }

  G4bool fLevchukFlag;
  G4int  fRCtype;
  G4bool fRemollMSFlag;
  G4double x1,x2,x3,x4,u1,u2,u3,u4,s;

  G4double fBeamRotZX = 0.00 * rad;
  G4double fBeamRotZY = 0.00 * rad;

private:
  MolPolEvent   *fDefaultEvent;
  G4ParticleGun *particleGun;
  MolPolIO      *fIO;
  MolPolDetectorConstruction *fDet;
  class MolPolMeradgen *fMeradgen;

  G4String rndmFlag;
  G4double angle;
  G4String gentype;

  void InitTargetMomentum();
  G4String fLevFilename = "../levchukMomentaBulkIron.tab";
  Int_t lexExpLines;
  std::vector< std::vector< G4double > > fLevchukData;

  void LevchukEffect();
  G4double fLEcorFac, fLEtgtPol;
  G4double SampleTargetMomentum(G4bool);

  G4double GetElectronStructFct(G4double&, const G4double);

  // Shared preamble, then one RC branch each (readable like gentype == "beam").
  void SampleMollerVertex(MollerVertex& v);
  void GenerateMollerBorn(G4Event* anEvent, const MollerVertex& v);
  void GenerateMollerAlexander(G4Event* anEvent, const MollerVertex& v);
  void GenerateMollerMeradgen(G4Event* anEvent, const MollerVertex& v);

  // Shared Born / Alexander kinematics+guns; applyStructureRC selects model.
  void GenerateMollerStructurePath(G4Event* anEvent, const MollerVertex& v,
                                   G4bool applyStructureRC);

  G4bool CheckLUNDFile(G4String LUNDfile_name);
  G4int    fNLUNDLines;
  std::ifstream LUNDfile;

  G4double fTargetA;
  G4double fTargetZ;
  G4double fTargetDensity;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
