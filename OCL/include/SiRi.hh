
#ifndef SiRi_h
#define SiRi_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include "G4RunManager.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
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

#include <stdlib.h>     /* exit, EXIT_FAILURE */

#include "CLHEP/Units/PhysicalConstants.h"

#include "G4Sphere.hh"

/// Detector construction class to define materials and geometry.
class G4LogicalVolume;
class G4VPhysicalVolume;

const G4int nFront = 8;
const G4int nBack = 8;

class SiRi : public G4VUserDetectorConstruction
{
  public:

    SiRi();
    ~SiRi();

  public:
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

	// void SetPosition( G4ThreeVector );
	// void SetRotation( G4RotationMatrix );
    void SetAngle(G4double);
	void Placement(G4int, G4VPhysicalVolume*, G4bool);
	void CreatePlacementParameters();

  private:
    G4ThreeVector        translatePos;
    // G4RotationMatrix     rotation;
    G4double angle;     // another way to give the angle wrt beam
    G4double theta;     // another way to give the angle wrt beam
    G4double anotherAngle;     // another way to give the angle wrt beam
    G4double pmone;

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

    G4double r_delta;

    G4double h_thin;

    G4double rInSiriHolder;
    G4double thicknessSiriHolder;
    G4double halfLengthSiriHolder;

    G4double rCableCu;
    G4double rCableAl;
    G4double halfLengthCable;
    G4double zTransCable;

    //
    // Materials
    //

    G4Material* aluminum;
    G4Material* copper;
    G4Material* silicon;
    G4Material* vacuum;


    // G4LogicalVolume*   deltaE_Log;
    G4Trd* DE_shape[nFront];
    G4Trd* dE_shape;
    G4Trd* pad_shape;

    G4LogicalVolume*   E_Log;
    G4LogicalVolume*   dE_Log;
    G4LogicalVolume*   DE_Log[nFront];
    G4LogicalVolume*   pad_Log;

    G4VPhysicalVolume*  E_phys;
    G4VPhysicalVolume*  dE_phys_strip;
    G4VPhysicalVolume*  dE_phys;
    G4VPhysicalVolume*  pad_phys;

    G4Tubs*  solidSiriHolder;
    G4LogicalVolume*  logSiriHolder;
    G4VPhysicalVolume*  physSiriHolder;

    G4Tubs* solidCableAl;
    G4Tubs*  solidCableCu;
    G4LogicalVolume*  logCableCu;
    G4LogicalVolume*  logCable;
    G4VPhysicalVolume*  physCableCu;
    G4VPhysicalVolume*  physCable;


  private:
  void CreateSolids();
};

#endif

