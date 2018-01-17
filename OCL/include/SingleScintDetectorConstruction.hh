
#ifndef SingleScintDetectorConstruction_h
#define SingleScintDetectorConstruction_h 1

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
class OCLCollimator;

const G4int numberOf_OCLLaBr3=32;

class SingleScintDetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    SingleScintDetectorConstruction();
    ~SingleScintDetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();
    
    //virtual void ConstructSDandField();

  private:
    
    void SetPlacementParameters();

    G4Box* solidWorld;
    G4LogicalVolume* WorldLog;
    G4PVPlacement* physiWorld;


    OCLLaBr3* labr3[numberOf_OCLLaBr3];
    OCLCollimator* collimator[numberOf_OCLLaBr3];

    bool OCLLaBr3_presence[numberOf_OCLLaBr3];   
    bool OCLCollimator_presence[numberOf_OCLLaBr3];
    G4double OCLLaBr3_Distance[numberOf_OCLLaBr3];
    G4double OCLLaBr3_theta[numberOf_OCLLaBr3];
    G4double OCLLaBr3_phi[numberOf_OCLLaBr3];
    G4RotationMatrix  rotmOCLLaBr3[numberOf_OCLLaBr3];
    G4ThreeVector positionOCLLaBr3[numberOf_OCLLaBr3];
    G4ThreeVector positionCollimator[numberOf_OCLLaBr3];
    
};

#endif

