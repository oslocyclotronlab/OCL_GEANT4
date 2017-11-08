#ifndef TChamberAuspuff_h
#define TChamberAuspuff_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include "G4RunManager.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Polycone.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4SubtractionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

//#include "G4PhysicalConstants.hh"
#include "CLHEP/Units/PhysicalConstants.h"

class G4VPhysicalVolume;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.

class TChamberAuspuff
{
  public:
    TChamberAuspuff();
    ~TChamberAuspuff();

    // virtual G4VPhysicalVolume* Construct();
    // virtual void ConstructSDandField();

  	// void SetPosition( G4ThreeVector );
  	// void SetRotation( G4RotationMatrix );
    void CreateSolids();
  	void Placement(G4int, G4VPhysicalVolume*, G4bool);
  	// void CreatePlacementParameters();

  protected:
    G4LogicalVolume*  fScoringVolume;

  private:
  	// G4ThreeVector        translatePos;
   //  G4RotationMatrix     rotation;

    //
    // Materials
    //
    G4Material* aluminium;
    G4Material* plexiGlass;

    //
    // Geometry
    //

    G4double rInnerTube           ;
    G4double thicknessMiddleTube  ;
    G4double thicknessOuterTube   ;
    G4double halfLengthTube       ;
    G4double thicknessEndCaps     ;
    G4double halfLengthEndCaps    ;


	//
	// Geometries
	// 

    G4Polycone* solidChamberTube;
    G4LogicalVolume* logicChamberTube;
    G4VPhysicalVolume* physiChamberTube;


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

