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

#include "G4VisAttributes.hh"
#include "G4Colour.hh"


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
  // void SetPosition( G4ThreeVector );
  // void SetRotation( G4RotationMatrix );
  void Placement(G4int, G4VPhysicalVolume*, G4bool);

private:
  // G4ThreeVector        translatePos;
  // G4RotationMatrix     rotation;


  //
  // Parameters
  //

  G4double dx_from_center; // Distance of global model "centre" to centre of ball/OSCAR

  G4double dx_FrameRing;  // arb number at the moment
  G4double dz_FrameRing; // arb number at the moment
  G4double dRot;

  G4double dx_FrameTopBase; // arb number at the moment
  G4double dy_FrameTopBase; // arb number at the moment
  G4double dRotFrameTopBase;

  //
  // Rotations and Transofrmations
  //
  G4ThreeVector offsetFrameBall;
  G4RotationMatrix rotmFrameBallCAD[2];
  G4ThreeVector translatFrameBallCAD[2];
  G4Transform3D transformFrameBallCAD[2];

  G4ThreeVector offsetFrameRing;
  G4RotationMatrix rotmFrameRingCAD[2];
  G4ThreeVector translatFrameRingCAD[2];
  G4Transform3D transformFrameRingCAD[2];

  G4ThreeVector offsetFrameTopBase;
  G4RotationMatrix rotmFrameTopBaseCAD[2];
  G4ThreeVector translatFrameTopBaseCAD[2];
  G4Transform3D transformFrameTopBaseCAD[2];

  //
  // Elements & Materials
  //

  G4Material* G4_Al_Material;

  //
  // Frame (Ball) Solids & Volumes
  //


  G4VSolid * FrameBallCADSolid;
  G4LogicalVolume * FrameBallCADlog;
  G4VPhysicalVolume * FrameBallCADphys[2];


  G4VSolid * FrameRingCADSolid;
  G4LogicalVolume * FrameRingCADlog;
  G4VPhysicalVolume * FrameRingCADphys[2];

  G4VSolid * FrameTopBaseCADSolid;
  G4LogicalVolume * FrameTopBaseCADlog;
  G4VPhysicalVolume * FrameTopBaseCADphys[2];

  //------------------------------------------------------
  // visualization attributes
  //------------------------------------------------------

  G4VisAttributes* FrameBallCAD_VisAtt;
  G4VisAttributes* FrameRingCAD_VisAtt;
  G4VisAttributes* FrameTopBaseCAD_VisAtt;

private:
  void CreateSolids();
  // void MakeMaterials();

};

#endif
