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

  G4double dx_FrameBase; // arb number at the moment
  G4double dy_FrameBase; // arb number at the moment
  G4double dRotFrameBase;

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

  G4ThreeVector offsetFrameBase;
  G4RotationMatrix rotmFrameBaseCAD[2];
  G4ThreeVector translatFrameBaseCAD[2];
  G4Transform3D transformFrameBaseCAD[2];

  //
  // Elements & Materials
  //

  G4Element* C ;
  G4Element* Mn;
  G4Element* Si;
  G4Element* P ;
  G4Element* S ;
  G4Element* N ;
  G4Element* Cu;
  G4Element* Cr;
  G4Element* Ni;
  G4Element* Fe;

  G4Material* SteelAISI304; // SteelAISI304 (Stainless Steel)
  G4Material* SteelS235JR; // SteelS235JR (Carbon Steel)
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

  G4VSolid * FrameBaseCADSolid;
  G4LogicalVolume * FrameBaseCADlog;
  G4VPhysicalVolume * FrameBaseCADphys[2];

  //------------------------------------------------------
  // visualization attributes
  //------------------------------------------------------

  G4VisAttributes* FrameBallCAD_VisAtt;
  G4VisAttributes* FrameRingCAD_VisAtt;
  G4VisAttributes* FrameBaseCAD_VisAtt;

private:
  void CreateSolids();
  // void MakeMaterials();

};

#endif
