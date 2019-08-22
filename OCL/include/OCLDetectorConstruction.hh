
#ifndef OCLDetectorConstruction_h
#define OCLDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4SystemOfUnits.hh"
#include "G4PVPlacement.hh"

/// Detector construction class to define materials and geometry.
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4PVPlacement;
class G4Box;

class OCLLaBr3;
class OCLLaBr3Messenger;
class OCLCollimator;
class OCLDetectorMessenger;

const G4int numberOf_OCLLaBr3=32;

class OCLDetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    OCLDetectorConstruction();
    ~OCLDetectorConstruction();

  public:
    G4VPhysicalVolume* Construct();

    // void UpdateGeometry();
    void SetLaBr3_Distance(G4int i, G4double dist);
    void SetUseLaBr3(G4int i, G4bool use);
    void SetUseCSGOldTargetChamber(G4bool use);
    void SetUseCSGOldTarget(G4bool use);
    void SetUseCSGRadSource(G4bool use);
    void SetUseCSGSiRi(G4bool use);
    void SetUseCSGNiff(G4bool use);
    //virtual void ConstructSDandField();

  private:


    void SetPlacementParameters();

    G4VPhysicalVolume* ConstructVolumes();

    OCLDetectorMessenger* fMessenger;
    OCLLaBr3Messenger* fMessenger_labr;

    G4Box* solidWorld;
    G4LogicalVolume* WorldLog;
    G4PVPlacement* physiWorld;
    G4double fworld_sizeXYZ;

    OCLLaBr3* labr3[numberOf_OCLLaBr3];
    OCLCollimator* collimator[numberOf_OCLLaBr3];

    bool fOCLLaBr3_presence[numberOf_OCLLaBr3];
    bool fOCLCollimator_presence[numberOf_OCLLaBr3];
    G4double fOCLLaBr3_Distance[numberOf_OCLLaBr3];
    G4double OCLLaBr3_theta[numberOf_OCLLaBr3];
    G4double OCLLaBr3_phi[numberOf_OCLLaBr3];
    G4RotationMatrix  rotmOCLLaBr3[numberOf_OCLLaBr3];
    G4ThreeVector positionOCLLaBr3[numberOf_OCLLaBr3];
    G4ThreeVector positionCollimator[numberOf_OCLLaBr3];

    bool fUseCSGOldTargetChamber;
    bool fUseCSGOldTarget;
    bool fUseCSGRadSource;
    bool fUseCSGSiRi;
    bool fUseCSGNiff;



};

#endif

