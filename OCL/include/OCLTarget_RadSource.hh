
#ifndef OCLTarget_RadSource_hh
#define OCLTarget_RadSource_hh 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include "G4RunManager.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Orb.hh"
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

    G4double amberliteRadius;
    //
    // Materials
    //
    G4Material* plexiglass;
    G4Material* amberlite;

    // Shapes

    G4Box* solidTarget;
    G4LogicalVolume* logTarget;
    G4VPhysicalVolume* physTarget;

    G4Orb* solidAmberlite;
    G4LogicalVolume* logAmberlite;
    G4VPhysicalVolume* physAmberlite;
};

#endif // OCLTarget_RadSource_hh

