#ifndef OCLLABR3_H
#define OCLLABR3_H 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4UnionSolid.hh"

#include "G4VisAttributes.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Colour.hh"


/// Detector construction class to define materials and geometry.

// class G4LogicalVolume;
// class G4VPhysicalVolume;
// class G4PVPlacement;
// class G4Box;


class OCLLaBr3 {

public:
  OCLLaBr3();
  ~OCLLaBr3();

public:
  void SetPosition( G4ThreeVector );
  void SetRotation( G4RotationMatrix );
  void Placement(G4int, G4VPhysicalVolume*, G4bool);

private:
  G4ThreeVector        translatePos;
  G4RotationMatrix     rotation;

  //
  // Elements & Materials
  //

  G4Material* Aluminium;
  G4Material* PlexiGlass;

  G4Material* Sili;
  G4Material* SiO2;
  G4Material* Na2O ;
  G4Material* K2O ;
  G4Material* Al2O3;
  G4Material* B2O3;

  G4Material* LaBr3;
  G4Material* CeBr3;
  G4Material* MgO;
  G4Material* LaBr3_Ce;
  G4Material* vacuum;

  G4Material* Borosilicate;
  G4Material* Bialkali;

  //
  // Material Properties (for scintillation physics, if enabled)
  //

  G4MaterialPropertiesTable* MgOMPT;
  G4MaterialPropertiesTable* LaBr3MPT;
  G4MaterialPropertiesTable* PlexiGlasMPT;
  G4MaterialPropertiesTable* QuartzMPT;
  G4MaterialPropertiesTable* K2CsSbMPT;
  G4MaterialPropertiesTable* vacMPT;


  //
  // Detector Solids & Volumes
  //

  G4Tubs* solidOCLDetector;
  G4LogicalVolume*  logicOCLDetector;
  G4VPhysicalVolume* physiOCLDetector;

  G4Tubs* solidShieldingMain;
  G4Cons* solidShieldingConical;
  G4Tubs* solidShieldingLid ;

  G4ThreeVector translationUnion1;
  G4UnionSolid* unionShielding1  ;
  G4ThreeVector translationUnion2;
  G4UnionSolid* unionShielding2;
  G4LogicalVolume* logicShielding;
  G4ThreeVector positionShielding;
  G4VPhysicalVolume* physiShield;

  G4ThreeVector positionCoating;
  G4Tubs* solidCoating;
  G4LogicalVolume* logicCoating;
  G4VPhysicalVolume* physiCoating;

  G4ThreeVector positionCoatingPlexi;
  G4Tubs* solidCoatingPlexi;
  G4LogicalVolume* logicCoatingPlexi;
  G4VPhysicalVolume* physiCoatingPlexi;

  G4ThreeVector positionReflector;
  G4Tubs* solidReflector;
  G4LogicalVolume* logicReflector;
  G4VPhysicalVolume* physiReflector;

  G4ThreeVector positionCrystal;
  G4Tubs* solidCrystal;
  G4LogicalVolume* logicCrystal;
  G4VPhysicalVolume* physiCrystal;


  G4Tubs* solidPlexiWindow;
  G4LogicalVolume* logicPlexiWindow;
  G4ThreeVector positionPlexiWindow;
  G4VPhysicalVolume* physiPlexiWindow;

  G4ThreeVector positionPMTWindow;
  G4Tubs* solidPMTWindow;
  G4LogicalVolume* logicPMTWindow;
  G4VPhysicalVolume* physiPMTWindow;

  G4ThreeVector positionCathode;
  G4Tubs* solidCathode;
  G4LogicalVolume* logicCathode;
  G4VPhysicalVolume* physiCathode;

  //------------------------------------------------------
  // Surfaces and boundary processes
  //------------------------------------------------------

  // Reflector - scintillator surface

  G4OpticalSurface* OpCryRefSurface;
  G4LogicalBorderSurface* CryRefSurface;
  G4OpticalSurface* OpCryPlexiSurface;
  G4LogicalBorderSurface* CryPlexiSurface;
  G4OpticalSurface* OpPlexiPMTWinSurface;
  G4LogicalBorderSurface* CryPlexiPMTWinSurface;
  G4OpticalSurface* OpPMTWinCathSurface;
  G4LogicalBorderSurface* PMTWinCathSurface;

  //------------------------------------------------------
  // visualization attributes
  //------------------------------------------------------

  // White color for Detector
  G4VisAttributes* VisAtt1;
  G4VisAttributes* VisAtt2;
  G4VisAttributes* VisAtt3;
  G4VisAttributes* VisAtt4;
  G4VisAttributes* VisAtt5;
  G4VisAttributes* VisAtt6;
  G4VisAttributes* VisAtt7;

private:
  void CreateSolids();
  // void MakeMaterials();

};

#endif
