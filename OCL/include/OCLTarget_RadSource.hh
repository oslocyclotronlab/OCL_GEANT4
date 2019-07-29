
#ifndef OCLTarget_RadSource_hh
#define OCLTarget_RadSource_hh 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include "G4RunManager.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Transform3D.hh"

#include "G4PhysicalConstants.hh"

#include <stdlib.h>     /* exit, EXIT_FAILURE */

/// Detector construction class to define materials and geometry.
class G4LogicalVolume;
class G4VPhysicalVolume;


class OCLTarget_RadSource
{
  public:

    OCLTarget_RadSource();
    ~OCLTarget_RadSource();

  public:
    // virtual G4VPhysicalVolume* Construct();
    // virtual void ConstructSDandField();

	// void SetPosition( G4ThreeVector );
	// void SetRotation( G4RotationMatrix );
    // void SetAngle(G4double);
    void CreateSolids();

	void Placement(G4int, G4VPhysicalVolume*, G4bool);

  private:
    //
    // Parameters/Constants
    //

     G4double pX2Target;
     G4double pY2Target;
     G4double pZ2Target;

     // G4double pX2TargetHole;
     // G4double pY2TargetHole;
     // G4double pZ2TargetHole;

     // G4double pX2TargetHolder;
     // G4double pY2TargetHolder;
     // G4double pZ2TargetHolder;

    //
    // Materials
    //

    G4Material* silicon;
    G4Material* vacuum;
    G4Material* copper;
    G4Material* aluminum;
    G4Material* pyrex;


    // Shapes


     G4Box* solidTarget;
     // G4Box* solidTargetHole;
     // G4Box* solidTargetHolder;

     G4LogicalVolume* logTarget;
     // G4LogicalVolume* logTargetHole;
     // G4LogicalVolume*  logTargetHolder;

     G4VPhysicalVolume* physTarget;
     // G4VPhysicalVolume* physTargetHole;
     // G4VPhysicalVolume* physTargetHolder;

  // private:
  // void CreateSolids();
};

#endif // OCLTarget_RadSource_hh

