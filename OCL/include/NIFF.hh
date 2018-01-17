#ifndef NIFF_h
#define NIFF_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include "G4RunManager.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
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

class NIFF : public G4VUserDetectorConstruction
{
  public:
    NIFF();
    virtual ~NIFF();

    virtual G4VPhysicalVolume* Construct();
    // virtual void ConstructSDandField();

  	// void SetPosition( G4ThreeVector );
  	// void SetRotation( G4RotationMatrix );
    void CreateSolids();
  	void Placement(G4int, G4VPhysicalVolume*, G4bool);
  	// void CreatePlacementParameters();
	
    
    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

  protected:
    G4LogicalVolume*  fScoringVolume;

  private:
  	G4ThreeVector        translatePos;
    G4RotationMatrix     rotation;

    //
    // Materials
    //
    G4Material* aluminium;
    G4Material* isobutane;
    G4Material* isobutane_ppac;
    G4Material* Mylar;
    G4Material* vacuum;


	//
	// Geometries
	// 
  	G4ThreeVector        pos1;
    G4RotationMatrix     rMat;



	G4LogicalVolume* trapLog;    
	G4LogicalVolume* trapAlLog;
	G4LogicalVolume* trapMylarLog;


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

