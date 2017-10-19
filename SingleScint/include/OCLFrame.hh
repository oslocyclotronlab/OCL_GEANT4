#ifndef OCLFRAME_H
#define OCLFRAME_H 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Polyhedra.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"

// #include "G4Colour.hh"


/// Detector construction class to define materials and geometry.

// class G4LogicalVolume;
// class G4VPhysicalVolume;
// class G4PVPlacement;
// class G4Box;

//////////////////////////////////////////////////////////
//                  Frame SETUP                         //
//////////////////////////////////////////////////////////

  const G4int numberOf_Pentagons = 12;
  const G4int numberOf_Hexagons  = 20;

class OCLFrame {

public:
  OCLFrame();
  ~OCLFrame();

public:
  void SetPosition( G4ThreeVector );
  void SetRotation( G4RotationMatrix );
  void Placement(G4int, G4VPhysicalVolume*, G4bool);
  void CreatePlacementParameters();

private:
  G4ThreeVector        translatePos;
  G4RotationMatrix     rotation;

  //
  //
  //
  // int test[] = {2,3}; 


  //
  // Rotations
  //

  G4double frameHexagon_theta[numberOf_Hexagons];
  G4double frameHexagon_phi[numberOf_Hexagons];
  G4double frameHexagon_psi[numberOf_Hexagons];
  G4RotationMatrix rotmHexagon[numberOf_Hexagons];
  G4ThreeVector positionHexagon[numberOf_Hexagons];

  G4double framePentagon_theta[numberOf_Pentagons];
  G4double framePentagon_phi[numberOf_Pentagons];
  G4double framePentagon_psi[numberOf_Pentagons];
  G4RotationMatrix rotmPentagon[numberOf_Pentagons];
  G4ThreeVector positionPentagon[numberOf_Pentagons];

  //
  // Elements & Materials
  //

  G4Material* Aluminium;

  //
  // Frame (Ball) Solids & Volumes
  //

  G4Tubs*             solidFrameHoles;

  G4Polyhedra*        solidPentagon;
  G4SubtractionSolid* subtractFramePentagon;
  G4LogicalVolume*    logicFramePentagon;
  G4VPhysicalVolume*  physiFramePentagon;

  G4Polyhedra*        solidHexagon;
  G4SubtractionSolid* subtractFrameHexagon;
  G4LogicalVolume*    logicFrameHexagon;
  G4VPhysicalVolume*  physiFrameHexagon;


private:
  void CreateSolids();
  // void MakeMaterials();

};

#endif
