
#ifndef SiRi_h
#define SiRi_h 1

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
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4PhysicalConstants.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TouchableHistory.hh"


#include "CLHEP/Units/PhysicalConstants.h"

/// Detector construction class to define materials and geometry.
class G4LogicalVolume;
class G4VPhysicalVolume;

class SiRi : public G4VUserDetectorConstruction
{
  public:

    SiRi();
    ~SiRi();

  public:
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

	void SetPosition( G4ThreeVector );
	void SetRotation( G4RotationMatrix );
	void Placement(G4int, G4VPhysicalVolume*, G4bool);
	void CreatePlacementParameters();

  private:
    G4ThreeVector        translatePos;
    G4RotationMatrix     rotation;

    //
    // Parameters/Constants
    //

    G4double dist;
    G4double a;
    G4double h;
    G4double alfa;
    G4double b;
    G4double d_deltaE;
    G4double d_E;

    G4double h_thin;


    //
    // Materials
    //

    G4Material* silicon;
    G4Material* vacuum;


    // G4LogicalVolume*   deltaE_Log;
    G4LogicalVolume*   E_Log;
    G4LogicalVolume*   DE_Log1;
    G4LogicalVolume*   DE_Log2;
    G4LogicalVolume*   DE_Log3;
    G4LogicalVolume*   DE_Log4;
    G4LogicalVolume*   DE_Log5;
    G4LogicalVolume*   DE_Log6;
    G4LogicalVolume*   DE_Log7;
    G4LogicalVolume*   DE_Log8;
    
    G4VPhysicalVolume*  E_phys;
    G4VPhysicalVolume*  DE_phys1;
    G4VPhysicalVolume*  DE_phys2;
    G4VPhysicalVolume*  DE_phys3;
    G4VPhysicalVolume*  DE_phys4;
    G4VPhysicalVolume*  DE_phys5;
    G4VPhysicalVolume*  DE_phys6;
    G4VPhysicalVolume*  DE_phys7;
    G4VPhysicalVolume*  DE_phys8;


  private:
  void CreateSolids();
};

#endif

