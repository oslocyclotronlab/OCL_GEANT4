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
#include "G4Polycone.hh"
#include "G4MultiUnion.hh"

#include "G4VisAttributes.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Colour.hh"


/// Detector construction class to define materials and geometry.

// class G4LogicalVolume;
// class G4VPhysicalVolume;
// class G4PVPlacement;
// class G4Box;

class OCLLaBr3Messenger;

class OCLLaBr3 {

public:
  OCLLaBr3(OCLLaBr3Messenger*);
  ~OCLLaBr3();

public:
  void SetPosition( G4ThreeVector );
  void SetRotation( G4RotationMatrix );
  void Placement(G4int, G4VPhysicalVolume*, G4bool);

// public:
  G4double GetCoatingPlasticThickness() { return coatingPlasticThickness; };
  // void SetCoatingAlThicknessFront(G4double thick);
  // void SetCoatingAlThickness(G4double thick);
  // void SetShieldingHalfThicknessLid(G4double thick);

private:
  G4ThreeVector        fTranslatePos;
  G4RotationMatrix     fRotation;

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
  G4Material* Air;

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

  G4Cons* solidShieldingConical;

  G4Tubs* solidShieldTubs[4];

  G4MultiUnion* unionShielding;

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

  G4Tubs* solidPMTandAir;
  G4LogicalVolume* logicPMTandAir;
  G4ThreeVector positionPMTandAir;
  G4VPhysicalVolume* physiPMTandAir;

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

  G4ThreeVector positionPMT;
  G4Polycone* solidPMT;
  G4LogicalVolume* logicPMT;
  G4VPhysicalVolume* physiPMT;

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


  //
  // Dimensions
  //

  G4double startPhi;
  G4double deltaPhi;

  G4double crystalOuterR;
  G4double crystalInnerR;
  G4double crystalHalfLength;

  // Some reflector like MgO
  G4double reflectorThickness;
  G4double reflectorHalfLength;
  G4double reflectorInnerR;
  G4double reflectorOuterR;

  // Alumium(?) coating
  G4double coatingAlThickness;
  G4double coatingAlThicknessFront;
  G4double coatingOuterR;

  // in between reflector and coating, there will be some plastic
  G4double coatingPlasticThickness;
  G4double coatingHalfLength;

  // G4double shieldingThickness1;
  // G4double shieldingThickness2;
  // G4double shieldingThickness3;

  G4double dShieldTubs[3];
  G4double rOuters[3];

  G4double shieldingHalfThicknessLid;
  G4double shieldingInnerR;
  G4double shieldingOuterR;

  //in the front, the shielding tube diameter is reduces. It's later modeled by a conical section
  G4double shieldingConeHalfLength;
  G4double shieldingConeOuterRFront;
  G4double shieldingConeOuterRBack;
  G4double shieldingHalfLength;

  G4double plexiGlasWindowOuterR;
  G4double plexiGlasWindowHalfLength;

  // PMT
  G4double PMTWindowHalfLength;
  G4double PMTWindowRadius;
  G4double cathodeHalfLength;
  G4double cathodeRadius;

  G4double PMTStartRadius;
  G4double PMTMidtRadius;
  G4double PMTEndRadius;
  G4double PMTMidtZ;
  G4double PMTEndZ;
  G4double PMTHalfLength;

  // Whole detector incl. PMT (-> Logical unit)
  G4double detectorHalfinclPMT;
  G4double PMTandAirHalfLength;

private:
  void GetMaterials();
  void CreateSolids();
  void SetOpticalProperties();
  void CalculateGeometryParameters();

  // OCLLaBr3Messenger* fMessenger;
  // void MakeMaterials();

};

#endif
