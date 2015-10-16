//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// -------------------------------------------------------------------
// $Id: DetectorConstruction.cc,v 1.5 2010-10-06 14:39:41 sincerti Exp $
// -------------------------------------------------------------------

#include "DetectorConstruction.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include "G4Transform3D.hh"
#include "G4PhysicalConstants.hh"

//#include "G4MultiFunctionalDetector.hh"
//#include "G4VPrimitiveScorer.hh"
//#include "G4PSEnergyDeposit.hh"
//#include "G4TransportationManager.hh"
//#include "G4SDManager.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::DetectorConstruction()
:
     solidWorld(0), WorldLog(0), WorldPhys(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::~DetectorConstruction()
{}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* DetectorConstruction::Construct()
{

    //----------------------------------------------------
	// Material definitions
	//----------------------------------------------------

  G4double a, z;                    //a=mass of a mole;
  G4double density;                 //z=mean number of protons;

  G4int ncomponents, natoms;
  G4double abundance, fractionmass;

  // load NIST material database manager
  G4NistManager * man = G4NistManager::Instance();

  //
  // Define Elements
  //

  // there were here already
    G4Element* H =  new G4Element("Hydrogen",    "H",   z=1.,   a=1.01*g/mole);
    G4Element* C =  new G4Element("Carbon",      "C",   z=6.,   a=12.01*g/mole);
    G4Element* F =  new G4Element("Fluorine",    "F",   z=9.,   a=18.9984*g/mole);
    G4Element* Br = new G4Element("Bromium",    "Br",   z=35.,  a=79.904*g/mole);
    G4Element* La = new G4Element("Lanthanum",  "La",   z=57.,  a=138.90547*g/mole); // TODO Test for Todo
    G4Element* Ce = new G4Element("Cerium",     "Cl",   z=58.,  a=140.116*g/mole);
    G4Element* Tl = new G4Element("Thallium",   "Tl",   z=81.,  a=204.383*g/mole);

  // add more elements from NIST database
    G4Element* B  = man->FindOrBuildElement("B");
    G4Element* O  = man->FindOrBuildElement("O");
    G4Element* Na = man->FindOrBuildElement("Na");
    G4Element* Al = man->FindOrBuildElement("Al");
    G4Element* Si = man->FindOrBuildElement("Si");
    G4Element* K  = man->FindOrBuildElement("K");
    G4Element* Sc = man->FindOrBuildElement("Sc");
    G4Element* Sb = man->FindOrBuildElement("Sb");
    G4Element* Cs = man->FindOrBuildElement("Cs");
    G4Element* Mg = man->FindOrBuildElement("Mg");

  //
  // define materials from elements.
  //

  // Sili
  G4Material* Sili =    new G4Material("Silicon", density= 2.330*g/cm3, ncomponents=1);
  Sili->AddElement(Si, natoms=1);

  // SiO2 = Quartz
  G4Material* SiO2 =    new G4Material("quartz", density= 2.200*g/cm3, ncomponents=2);
  SiO2->AddElement(Si, natoms=1);
  SiO2->AddElement(O , natoms=2);

  // Na2O
  G4Material* Na2O =    man->FindOrBuildMaterial("G4_SODIUM_MONOXIDE");

  // K2O
  G4Material* K2O =     man->FindOrBuildMaterial("G4_POTASSIUM_OXIDE");

  // Al2O3
  G4Material* Al2O3 =   man->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");

  //B2O3
  G4Material* B2O3 =    man->FindOrBuildMaterial("G4_BORON_OXIDE");

  //Lead
  G4Material* lead =    man->FindOrBuildMaterial("G4_Pb");

  //LaBr3
  G4Material* LaBr3 =   new G4Material("LaBr3", density = 5.07*g/cm3, ncomponents=2);
  LaBr3->AddElement(La, natoms=1);
  LaBr3->AddElement(Br, natoms=3);

  // MgO reflector
  density = 2.0*g/cm3;
  G4Material* MgO = new G4Material("MgO", density, ncomponents=2);
  MgO->AddElement(Mg, natoms=1);
  MgO->AddElement(O, natoms=1);

  //LaBr3_Ce
  //with 5% dopping, see technical note "BrilLanCe Scintillators Performance Summary"
  G4Material* LaBr3_Ce = new G4Material("LaBr3_Ce", density = 5.08*g/cm3, ncomponents=2);
  LaBr3_Ce->AddMaterial(LaBr3,  fractionmass=95*perCent);
  LaBr3_Ce->AddElement(Ce,      fractionmass=5*perCent);


   // vacuum (non-STP)

    G4Material* vacuum = new G4Material("Vacuum",       //name as String
							1,		                    //atomic number (use 1 for Hydrogen)
                    		1.008*g/mole, 	            //molar mass (use 1.008*g/mole for Hydoren)
							1.e-25*g/cm3,  	            //density
							kStateGas,		            //kStateGas - the material is gas (see G4State)
                    		2.73*kelvin,	            //Temperature
							1.e-25*g/cm3);	            //pressure


	// Steel as non-NIST material
    G4Element* elFe = G4NistManager::Instance()->FindOrBuildElement("Fe");
    G4Element* elNi = G4NistManager::Instance()->FindOrBuildElement("Ni");
    G4Element* elCr = G4NistManager::Instance()->FindOrBuildElement("Cr");
    G4Material* iron = new G4Material("StainlessSteel", 7.80 * g/cm3, 3 /* components */);
    iron -> AddElement(elFe, 70 * perCent);
    iron -> AddElement(elCr, 18 * perCent);
    iron -> AddElement(elNi, 12 * perCent);


    // PMT-materials

    // Borosilicate
    G4Material* Borosilicate = new G4Material("Borosilicate glass", density= 2.23*g/cm3, ncomponents=5);
    Borosilicate->AddMaterial(SiO2,   fractionmass=80.6 * perCent);
    Borosilicate->AddMaterial(B2O3,  fractionmass=13.0 * perCent);
    Borosilicate->AddMaterial(Na2O,  fractionmass=2.   * perCent); // 1/2 of wt% for (Na20+K20)
    Borosilicate->AddMaterial(K2O,   fractionmass=2.   * perCent); // 1/2 of wt% for (Na20+K20)
    Borosilicate->AddMaterial(Al2O3, fractionmass=2.31  * perCent);

    // Bialkali
    // (Bialkali KCsSb,  Density=?, Thickness=?)?
    G4Material* Bialkali = new G4Material("Bialkali", density= 2*g/cm3, ncomponents=3);
	Bialkali->AddElement(K,  natoms=2);
	Bialkali->AddElement(Cs, natoms=1);
	Bialkali->AddElement(Sb, natoms=1);


	//------------------------------------------------------
	// Optical properties
	//------------------------------------------------------


    // at the moment taken from the Scintiallator example
  	const G4int nEntries = 2;

	G4double PhotonEnergy[nEntries] = {1.0*eV,7.0*eV};

	// MgO reflector

	G4double MgORefractionIndex[nEntries] = {1.0,1.0};

	G4double MgOAbsorptionLength[nEntries] = {1.0E-9*m,1.0E-9*m};

	G4MaterialPropertiesTable* MgOMPT = new G4MaterialPropertiesTable();

	MgOMPT->AddProperty("RINDEX",PhotonEnergy,MgORefractionIndex,
						nEntries);
	MgOMPT->AddProperty("ABSLENGTH",PhotonEnergy,MgOAbsorptionLength,
						nEntries);

	MgO->SetMaterialPropertiesTable(MgOMPT);


  	// LaBr3

	G4double LaBr3RefractionIndex[nEntries] = {1.9,1.9};

	G4double LaBr3AbsorptionLength[nEntries] = {50.*cm,50.*cm};

	G4MaterialPropertiesTable* LaBr3MPT = new G4MaterialPropertiesTable();

	LaBr3MPT->AddProperty("RINDEX",PhotonEnergy,LaBr3RefractionIndex,
						  nEntries);
	LaBr3MPT->AddProperty("ABSLENGTH",PhotonEnergy,LaBr3AbsorptionLength,
						  nEntries);

	G4double ScintEnergy[nEntries] = {3.26*eV,3.44*eV};
	G4double ScintFast[nEntries] = {1.0,1.0};

	LaBr3MPT->AddProperty("FASTCOMPONENT",ScintEnergy,ScintFast,nEntries);

	LaBr3MPT->AddConstProperty("RESOLUTIONSCALE",1.);
	LaBr3MPT->AddConstProperty("FASTTIMECONSTANT",20.*ns);
	LaBr3MPT->AddConstProperty("YIELDRATIO",1.);
	LaBr3MPT->AddConstProperty("SCINTILLATIONYIELD",63./keV);     // manifact. info

	LaBr3_Ce->SetMaterialPropertiesTable(LaBr3MPT);

	// Quartz

	G4double QuartzRefractionIndex[nEntries] = {1.47,1.47};

	G4double QuartzAbsorptionLength[nEntries] = {3.0*cm,3.0*cm};

	G4MaterialPropertiesTable* QuartzMPT = new G4MaterialPropertiesTable();

	QuartzMPT->AddProperty("RINDEX",PhotonEnergy,QuartzRefractionIndex,
						   nEntries);
	QuartzMPT->AddProperty("ABSLENGTH",PhotonEnergy,QuartzAbsorptionLength,
						   nEntries);

	SiO2->SetMaterialPropertiesTable(QuartzMPT);

	// K2CsSb (Bialcali Photocathode)

	G4double K2CsSbRefractionIndex[nEntries] = {1.47,1.47};

	G4double K2CsSbAbsorptionLength[nEntries] = {1.0E-9*m,1.0E-9*m};

	G4MaterialPropertiesTable* K2CsSbMPT = new G4MaterialPropertiesTable();

	K2CsSbMPT->AddProperty("RINDEX",PhotonEnergy,K2CsSbRefractionIndex,
						   nEntries);
	K2CsSbMPT->AddProperty("ABSLENGTH",PhotonEnergy,K2CsSbAbsorptionLength,
						   nEntries);

	Bialkali->SetMaterialPropertiesTable(K2CsSbMPT);

	// Vacuum

	G4double vacRefractionIndex[nEntries] = {1.0,1.0};


	G4MaterialPropertiesTable* vacMPT = new G4MaterialPropertiesTable();
	vacMPT->AddProperty("RINDEX",PhotonEnergy,vacRefractionIndex,
						nEntries);


	vacuum->SetMaterialPropertiesTable(vacMPT);






  	//------------------------------------------------------
	// Detector geometry
	//------------------------------------------------------

	//
	// World
	//

  G4double world_sizeXYZ = 35*cm;

  G4Box* solidWorld = new G4Box("World",
                         world_sizeXYZ/2, world_sizeXYZ/2, world_sizeXYZ/2);     				//size (defined through half-sizes)

  G4LogicalVolume* WorldLog =  new G4LogicalVolume(solidWorld,        		//solid defining the World
                        		  vacuum,           	//material of the World
                        		  "World");         	//name

  G4VPhysicalVolume* WorldPhys = new G4PVPlacement(0,       //specifies rotation: 0 = no rotation
                      			    G4ThreeVector(),     	//at (0,0,0)
                      				WorldLog,            	//logical volume
									"World",               	//name
									0,                     	//mother  volume
									false,                 	//no boolean operation
									0);                     //copy number




  	//
	// Detector
	//


	// Reflector
	// Assumption: We (currently) don't know whether we really have a Reflector/Material properiets
	//             the specifications are taken from the Scintillation example

	G4double crystalOuterR = 50.*mm;
  	G4double crystalInnerR = 0.0*mm;
  	G4double crystalHalfLength = 213.5*0.5*mm;
  	G4double startAngleOfTheCrystal = 0.*deg;
  	G4double spanningAngleOfTheCrystal = 360.*deg;

  	G4double ReflectorThickness = 0.5*mm;
	G4double ReflectorHalfLength = crystalHalfLength + ReflectorThickness;
	//G4double ReflectorInnerR = crystalOuterR;
	G4double ReflectorInnerR = 0;
	G4double ReflectorOuterR = crystalOuterR + ReflectorThickness;
	G4double StartPhi = 0.*deg;
	G4double DeltaPhi = 360.*deg;

	G4Tubs* solidReflector = new G4Tubs("Reflector",
										ReflectorInnerR, ReflectorOuterR, ReflectorHalfLength,
										StartPhi,DeltaPhi);

	G4LogicalVolume* logicReflector = new G4LogicalVolume(solidReflector,MgO,
														  "Reflector");

	G4ThreeVector positionReflector = G4ThreeVector(0.*cm,0.*cm,0.*cm);

	G4VPhysicalVolume* physiReflector = new G4PVPlacement(0,positionReflector,
														  "Reflector",logicReflector,
														  WorldPhys,false,0);

	// Crystal

  	G4Tubs* detectorSolid = new G4Tubs("Crystal",
				 			crystalInnerR, crystalOuterR, crystalHalfLength,
                                 		startAngleOfTheCrystal,spanningAngleOfTheCrystal);

    G4LogicalVolume* detectorLog = new G4LogicalVolume(detectorSolid, LaBr3_Ce, "Crystal");

    G4ThreeVector positionCrystal = G4ThreeVector(0.*cm,0.*cm,ReflectorThickness);

	G4VPhysicalVolume* physiCrystal = new G4PVPlacement(0, positionCrystal, "Crystal",
	                                                    detectorLog, physiReflector, false, 0);



    // Shielding

    G4double shieldOuterR = 60.*mm + ReflectorThickness;
  	G4double shieldInnerR = 50.*mm + ReflectorThickness;
  	G4double shieldHeight = 213.5*0.5*mm;
  	G4double startAngleOfTheShield = 0.*deg;
  	G4double spanningAngleOfTheShield = 360.*deg;


  	G4Tubs* shieldSolid = new G4Tubs("Shield",
				 			shieldInnerR, shieldOuterR, shieldHeight,
                                 		startAngleOfTheShield,spanningAngleOfTheShield);

    G4LogicalVolume* shieldLog = new G4LogicalVolume(shieldSolid, iron, "Shield");

	G4VPhysicalVolume* physiShield = new G4PVPlacement(0, G4ThreeVector(), "Shield", shieldLog, WorldPhys, false, 0);

   // define lid - addition of two disks

	G4double lidOuterR = 60.*mm;
	G4double lidFullInnerR = 0.*mm;
	G4double lidDiskInnerR = 90.5*0.5*mm;
	G4double lidHeight = 2.5*mm;
	G4double startAngle = 0.*deg;
  	G4double spanningAngle = 360.*deg;

  	G4ThreeVector translate = G4ThreeVector(0, 0, -5*mm);
  	G4RotationMatrix rMat;									// rotation matrix around X axis
  		rMat.rotateX(0.*deg);
  	G4Transform3D transform(rMat,translate);

  	G4Tubs* lidFullSolid = new G4Tubs("LidFull", lidFullInnerR, lidOuterR, lidHeight, startAngle, spanningAngle);
  	G4Tubs* lidDiskSolid = new G4Tubs("LidDisk", lidDiskInnerR, lidOuterR, lidHeight, startAngle, spanningAngle);
  	G4UnionSolid* lidSolid = new G4UnionSolid("Lid", lidFullSolid, lidDiskSolid, transform);

  	G4LogicalVolume* lidLog = new G4LogicalVolume( lidSolid, iron, "Lid");

  	G4ThreeVector positionLid = G4ThreeVector(0.,0., -(crystalHalfLength + ReflectorThickness)  );

	G4VPhysicalVolume* physiLid = new G4PVPlacement(0, positionLid, "Lid", lidLog, WorldPhys, false, 0);


	//
	// PMT
	//

	G4double PMTWindowHalfLength = 1.0*mm;
	G4double PMTWindowRadius = 85*0.5*mm;

	G4double CathodeHalfLength = 0.005*mm;
	G4double CathodeRadius =85*0.5*mm;


    // PMT window

	G4Tubs* solidPMTWindow = new G4Tubs("PMTWindow",0.*cm,PMTWindowRadius,
										PMTWindowHalfLength,StartPhi,DeltaPhi);

	G4LogicalVolume* logicPMTWindow = new G4LogicalVolume(solidPMTWindow,
														  SiO2,"PMTWindow");

	G4ThreeVector positionPMTWindow = G4ThreeVector(0.*cm,0.*cm,
													crystalHalfLength+PMTWindowHalfLength);

	G4VPhysicalVolume* physiPMTWindow = new G4PVPlacement(0,positionPMTWindow,
														  "PMTWindow",logicPMTWindow,
														  WorldPhys,false,0);

	// Photocathode

	G4Tubs* solidCathode = new G4Tubs("Cathode",0.*cm,CathodeRadius,
									  CathodeHalfLength,StartPhi,DeltaPhi);

	G4LogicalVolume* logicCathode = new G4LogicalVolume(solidCathode,
														Bialkali,"Cathode");

	G4ThreeVector positionCathode = G4ThreeVector(0.*cm,0.*cm,
												  crystalHalfLength+2.*PMTWindowHalfLength
												  +CathodeHalfLength);

	G4VPhysicalVolume* physiCathode = new G4PVPlacement(0,positionCathode,
														"Cathode",logicCathode,
														WorldPhys,false,0);


	//------------------------------------------------------
	// Surfaces and boundary processes
	//------------------------------------------------------


	// Reflector - sintillator surface

			//some sort of reflector!
		// Reflector - sintillator surface
		// HELP HELP HELP!! This needs to be reviewed!

	G4OpticalSurface* OpRefCrySurface =
	new G4OpticalSurface("RefCrySurface");

	OpRefCrySurface->SetType(dielectric_metal);
	OpRefCrySurface->SetModel(glisur);
	OpRefCrySurface->SetFinish(polished);

	G4LogicalBorderSurface* RefCrySurface =
    new G4LogicalBorderSurface("RefCrySurface",physiCrystal,
							   physiReflector,OpRefCrySurface);



//	G4OpticalSurface* OpRefCrySurface =
//	new G4OpticalSurface("RefCrySurface");
//
//	OpRefCrySurface->SetType(dielectric_metal);
//	OpRefCrySurface->SetModel(glisur);
//	OpRefCrySurface->SetFinish(polished);
//
//	G4LogicalBorderSurface* RefCrySurface =
//    new G4LogicalBorderSurface("RefCrySurface",physiCrystal,
//							   physiShield,OpRefCrySurface);


	// Scintillator - PMT window surface

	G4OpticalSurface* OpCryPMTWinSurface =
	new G4OpticalSurface("CryPMTWinSurface");

	OpCryPMTWinSurface->SetType(dielectric_dielectric);
	OpCryPMTWinSurface->SetModel(glisur);
	OpCryPMTWinSurface->SetFinish(polished);

	G4LogicalBorderSurface* CryPMTWinSurface =
    new G4LogicalBorderSurface("CryPMTWinSurface",physiCrystal,physiPMTWindow,
							   OpCryPMTWinSurface);

	// PMT window - photocathode surface

	G4OpticalSurface* OpPMTWinCathSurface = new G4OpticalSurface("PMTWinCathSurface");

	OpPMTWinCathSurface->SetType(dielectric_dielectric);
	OpPMTWinCathSurface->SetModel(glisur);
	OpPMTWinCathSurface->SetFinish(polished);

	G4LogicalBorderSurface* PMTWinCathSurface =
    new G4LogicalBorderSurface("CathodeSurface",physiPMTWindow,physiCathode,
							   OpPMTWinCathSurface);






	//------------------------------------------------------
	// visualization attributes
	//------------------------------------------------------

//  detectorLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  G4VisAttributes* worldVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0)); //White
  worldVisAtt->SetVisibility(true);
  detectorLog->SetVisAttributes(worldVisAtt);

  G4VisAttributes* worldVisAtt1 = new G4VisAttributes(G4Colour(1.0,0.0,0.0)); //red
  worldVisAtt1->SetVisibility(true);
  shieldLog->SetVisAttributes(worldVisAtt1);
//    shieldLog->SetVisAttributes(G4VisAttributes::GetInvisible());

  G4VisAttributes* worldVisAtt2 = new G4VisAttributes(G4Colour(0.0,1.0,0.0)); //red
  worldVisAtt2->SetVisibility(true);
  lidLog->SetVisAttributes(worldVisAtt2);

	//Blue color for PMT window
	G4VisAttributes* Att3= new G4VisAttributes(G4Colour(0.0,0.0,1.0));
	logicPMTWindow->SetVisAttributes(Att3);

	//White color for the absorber photocathode
	G4VisAttributes* Att4= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
	logicCathode->SetVisAttributes(Att4);

		//Yellow color for reflector
	G4VisAttributes* Att2= new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	logicReflector->SetVisAttributes(Att2);






	//
	// always return the physical World
	//

  return WorldPhys;
}

//void DetectorConstruction::ConstructSDandField()
//{
//  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

//// declare trackers as a MultiFunctionalDetector scorer
//  //
//  G4MultiFunctionalDetector* ScintDet = new G4MultiFunctionalDetector("LaBrScint");
//  G4VPrimitiveScorer* primitiv1 = new G4PSEnergyDeposit("edep1");
//  ScintDet->RegisterPrimitive(primitiv1);
//  SetSensitiveDetector("Crystal",ScintDet);
//   //

//}
