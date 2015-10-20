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


#include "DetectorConstruction.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4UnionSolid.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include "G4Transform3D.hh"
#include "G4PhysicalConstants.hh"

#include "G4Colour.hh"

//#include "G4MultiFunctionalDetector.hh"
//#include "G4VPrimitiveScorer.hh"
//#include "G4PSEnergyDeposit.hh"
//#include "G4TransportationManager.hh"
//#include "G4SDManager.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::DetectorConstruction()
//:
//     solidWorld(0), WorldLog(0), WorldPhys(0)
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

	// Aluminium
	//a = 26.98*g/mole;
	G4Material* Aluminium = new G4Material("Aluminum", density = 2.700*g/cm3, ncomponents=1);
	Aluminium->AddElement(Al, natoms=1);

	// Plexiglas
	G4Material* PlexiGlass = man->FindOrBuildMaterial("G4_PLEXIGLASS");

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

	G4double PhotonEnergy[nEntries] = {1.0*eV,7.0*eV}; // 1eV -> 1.2399 µm; 7eV -> 0.1771µm // TODO more detailed; adopt all of them

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



	// PlexiGlas

	G4double PhotonEnergy_PlexiGlass[nEntries] = {3.2626*eV, 3.4439*eV}; // 3.4439eV <- 360nm; 3.2626eV <- 380nm

	//info from: http://refractiveindex.info/?shelf=organic&book=poly%28methyl_methacrylate%29&page=Szczurowski
	// Szczurowski 2013 - n 0.4047-1.083 µm; extrapolated with
	// n=\sqrt( 1+\frac{0.99654λ^2}{λ^2-0.00787}+\frac{0.18964λ^2}{λ^2-0.02191}
	//          +\frac{0.00411λ^2}{λ^2-3.85727} ), where λ is in µm
	G4double PlexiGlasRefractionIndex[nEntries] = {1.47996,1.47996};

	// values from
	// Polycast Acrylic Sheets. Davis Earle, Ron Deal and Earl Gaudette. 	SNO-STR-93-042	revised and expanded Jan 24, 1994
	G4double PlexiGlasAbsorptionLength[nEntries] = {1./(0.04e-2)*m,1./(0.02e-2)*m};

	G4MaterialPropertiesTable* PlexiGlasMPT = new G4MaterialPropertiesTable();

	MgOMPT->AddProperty("RINDEX",PhotonEnergy,MgORefractionIndex,
						nEntries);
	MgOMPT->AddProperty("ABSLENGTH",PhotonEnergy,MgOAbsorptionLength,
						nEntries);

	PlexiGlass->SetMaterialPropertiesTable(PlexiGlasMPT);

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

  G4double world_sizeXYZ = 50*cm;

  G4Box* solidWorld = new G4Box("World",
                         world_sizeXYZ/2, world_sizeXYZ/2, world_sizeXYZ/2);     				//size (defined through half-sizes)

  G4LogicalVolume* WorldLog =  new G4LogicalVolume(solidWorld,        		//solid defining the World
                        		  vacuum,           	//material of the World
                        		  "World");         	//name

  G4VPhysicalVolume* physiWorld = new G4PVPlacement(0,       //specifies rotation: 0 = no rotation
                      			    G4ThreeVector(),     	//at (0,0,0)
                      				WorldLog,            	//logical volume
									"World",               	//name
									0,                     	//mother  volume
									false,                 	//no boolean operation
									0);                     //copy number




    // parameters

	G4double crystalOuterR = 8.89*cm/2.; // 3.5 in in cm
  	G4double crystalInnerR = 0.0*mm;
  	G4double crystalHalfLength = 203.2*0.5*mm; // 8 in in cm
  	G4double startPhi = 0.*deg;
  	G4double deltaPhi = 360.*deg;

  	G4double reflectorThickness = 1.*mm; // assumption: 1 mm thick reflector on the front side
	G4double reflectorHalfLength = crystalHalfLength + 0.5 * reflectorThickness; // assumption: backside doesn't have a reflector
	//G4double ReflectorInnerR = crystalOuterR;
	G4double reflectorInnerR = 0.*mm;
	G4double reflectorOuterR = crystalOuterR + reflectorThickness;

  	G4double coatingThickness = 2.*mm; // thickness as in the radius part
  	G4double coatingThicknessFront = 1.*mm; // we assume a smaller thickness at the front of the detector
  	G4double coatingOuterR = 100.*mm/2. ;
	// in between reflector and coating, there will be some plastic
	G4double coatingPlasticThickness = coatingOuterR - coatingThickness -reflectorThickness - crystalOuterR; // assumption: 2.55 mm plexiglas coating around the reflector before the aluminium
  	G4double coatingHalfLength = reflectorHalfLength + 0.5 * coatingThicknessFront + 0.5 * coatingPlasticThickness; // backside doesn't have an (Aluminium) coating

  	G4double shieldingThickness = 5.*mm; 		// thickness of the tube
  	G4double shieldingHalfThicknessLid = 2.*mm/2.;
  	G4double shieldingInnerR = 0*mm; 			// as we use it as a mother volume
  	G4double shieldingOuterR = coatingOuterR + shieldingThickness;

  	//in the front, the shielding tube diameter is reduces. It's later modeled by a conical section
  	G4double shieldingConeHalfLength = 10.*mm;// in the front, the tube
  	//G4double shieldingConeInnerRFront = coatingOuterR;
  	G4double shieldingConeOuterRFront = coatingOuterR + 2.*mm;
  	//G4double shieldingConeInnerRBack = shieldingConeInnerRFront;
  	G4double shieldingConeOuterRBack = coatingOuterR + 5.*mm;

  	G4double shieldingHalfLength = coatingHalfLength - shieldingConeHalfLength; // without conical Section and Lid
  																						 //  we assume no coating at the back side

  	G4double plexiGlasWindowOuterR = shieldingOuterR; // currently we just assume a flat window on the top.
  	G4double plexiGlasWindowHalfLength= 0.5 * 1.*mm;



  	//
	// Detector Geometry
	//

  	//
    // Shielding
  	//

  	// Main tube

  	G4Tubs* solidShieldingMain = new G4Tubs("ShieldingMainTube",
											shieldingInnerR,
											shieldingOuterR,
											shieldingHalfLength,
											startPhi,
											deltaPhi);

	//G4LogicalVolume* logicShieldingMain = new G4LogicalVolume(solidShieldingMain, Aluminium, "ShieldingMainTube");

	// concial section

	G4Cons* solidShieldingConical = new G4Cons("ShieldingConnical",
										shieldingInnerR,   // inner radius = 0 because used as mother volume
										shieldingConeOuterRFront,
										shieldingInnerR,   // inner radius = 0 because used as mother volume
										shieldingConeOuterRBack,
										shieldingConeHalfLength,
										startPhi,
										deltaPhi);

	//G4LogicalVolume* logicShieldConical = new G4LogicalVolume(solidShieldingConical, Aluminium, "ShieldingConical");

	// lid

  	G4Tubs* solidShieldingLid = new G4Tubs("ShieldingLid",
  									shieldingInnerR,
  									shieldingConeOuterRFront, // we assume that there is no gap between the coating/crystal and shielding
  									shieldingHalfThicknessLid,
  									startPhi,
  									deltaPhi);

	//G4LogicalVolume* logicShieldingLid = new G4LogicalVolume(solidShieldingLid, Aluminium, "ShieldingLid");

	// add the shielding parts together

  	// The origin and the coordinates of the combined solid are the same as those of
  	// the first solid.

	G4ThreeVector translationUnion1 = G4ThreeVector(0, 0, - (shieldingConeHalfLength+shieldingHalfThicknessLid) );
    G4UnionSolid* unionShielding1 = new  G4UnionSolid("unionShielding1",
									solidShieldingConical,  // 1st object
									solidShieldingLid,	   // 2nd object
                    				0,				// no Rotation
                    				translationUnion1);   // translation of the second object

	G4ThreeVector translationUnion2 = G4ThreeVector(0, 0, -(shieldingHalfLength+shieldingConeHalfLength));
    G4UnionSolid* unionShielding2 = new  G4UnionSolid ("unionShielding2",
									solidShieldingMain,		     // 1nd object
									unionShielding1,			// 2st object
                    				0,				// no Rotation
                    				translationUnion2);   // translation of the second object

    G4LogicalVolume* logicShielding = new G4LogicalVolume(unionShielding2, Aluminium, "ShieldingLid");

    G4ThreeVector positionShielding = G4ThreeVector(0.*cm,
    											    0.*cm,
    											    -(coatingHalfLength-shieldingConeHalfLength-crystalHalfLength) // because of the shift in the coordinate system of the coating
    												);  														   // (center != origin)

	G4VPhysicalVolume* physiShield = new G4PVPlacement(0, positionShielding, "Shielding", logicShielding, physiWorld, false, 0);


	//
	// Detector Crystal (incl. Reflector & Coating)
	//

	// Coating
	// Coating: Aluminum part

	G4Tubs* solidCoating = new G4Tubs("Coating",
	  									0. * mm, // inner radius = 0 because used as mother volume
	  									coatingOuterR,
	  									coatingHalfLength,
	  									startPhi,
	  									deltaPhi);

	G4LogicalVolume* logicCoating = new G4LogicalVolume(solidCoating, Aluminium, "Coating");

	G4ThreeVector positionCoating = G4ThreeVector(0.*cm,0.*cm,-shieldingConeHalfLength); // because of the shift in the coordinate system of the shielding
																						 // (center != origin)

	G4VPhysicalVolume* physiCoating = new G4PVPlacement(0,positionCoating,
															  "Coating",logicCoating,
															  physiShield,false,0);
	// Coating: PlexiGlass part
	G4Tubs* solidCoatingPlexi = new G4Tubs("CoatingPlexiGlas",
	  									0. * mm, // inner radius = 0 because used as mother volume
	  									coatingOuterR-coatingThickness,
	  									coatingHalfLength - 0.5 * coatingThicknessFront, // in order to get the coatingPlasticHalfLength
	  									startPhi,
	  									deltaPhi);

	G4LogicalVolume* logicCoatingPlexi = new G4LogicalVolume(solidCoatingPlexi, PlexiGlass, "CoatingPlexiGlas");

	G4ThreeVector positionCoatingPlexi = G4ThreeVector(0.*cm,
													   0.*cm,
													   0.5 * coatingThicknessFront); // because of the shift in the coordinate system of the coating
																					 // (center != origin)

	G4VPhysicalVolume* physiCoatingPlexi = new G4PVPlacement(0,positionCoatingPlexi,
															  "CoatingPlexiGlas",logicCoatingPlexi,
															  physiCoating,false,0);

	// Reflector
	// Assumption: We (currently) don't know whether we really have a Reflector/Material properties
	//             the specifications are taken from the Scintillation example

	G4Tubs* solidReflector = new G4Tubs("Reflector",
										reflectorInnerR, reflectorOuterR, reflectorHalfLength,
										startPhi,deltaPhi);

	G4LogicalVolume* logicReflector = new G4LogicalVolume(solidReflector,MgO,
														  "Reflector");

	G4ThreeVector positionReflector = G4ThreeVector(0.*cm,
													0.*cm,
													0.5 * coatingPlasticThickness
													);

	G4VPhysicalVolume* physiReflector = new G4PVPlacement(0,positionReflector,
														  "Reflector",logicReflector,
														  physiCoatingPlexi,false,0);

	// Crystal

  	G4Tubs* solidCrystal = new G4Tubs("Crystal",
				 			crystalInnerR, crystalOuterR, crystalHalfLength,
                                 		startPhi,deltaPhi);

    G4LogicalVolume* logicCrystal = new G4LogicalVolume(solidCrystal, LaBr3_Ce, "Crystal");

    G4ThreeVector positionCrystal = G4ThreeVector(0.*cm,0.*cm, 0.5*reflectorThickness );

	G4VPhysicalVolume* physiCrystal = new G4PVPlacement(0, positionCrystal, "Crystal",
	                                                    logicCrystal, physiReflector, false, 0);

	// Plexiglas Window on Detector
	G4Tubs* solidPlexiWindow = new G4Tubs("PlexiGlasWindow",
	  									0. * mm, // inner radius = 0
	  									plexiGlasWindowOuterR,
	  									plexiGlasWindowHalfLength, // in order to get the coatingPlasticHalfLength
	  									startPhi,
	  									deltaPhi);

	G4LogicalVolume* logicPlexiWindow = new G4LogicalVolume(solidPlexiWindow, PlexiGlass, "PlexiGlasWindow");

	G4ThreeVector positionPlexiWindow = G4ThreeVector(0.*cm,
													  0.*cm,
													  crystalHalfLength+plexiGlasWindowHalfLength
													  ); // because of the shift in the coordinate system of the coating
																					 // (center != origin)

	G4VPhysicalVolume* physiPlexiWindow = new G4PVPlacement(0,positionPlexiWindow,
			"PlexiGlasWindow",logicPlexiWindow,physiWorld,false,0);

	//
	// PMT
	//

	G4double PMTWindowHalfLength = 1.0*mm;
	G4double PMTWindowRadius = 85*0.5*mm;

	G4double cathodeHalfLength = 0.005*mm;
	G4double cathodeRadius =85*0.5*mm;


    // PMT window

	G4Tubs* solidPMTWindow = new G4Tubs("PMTWindow",0.*cm,PMTWindowRadius,
										PMTWindowHalfLength,startPhi,deltaPhi);

	G4LogicalVolume* logicPMTWindow = new G4LogicalVolume(solidPMTWindow,
														  SiO2,"PMTWindow");

	G4ThreeVector positionPMTWindow = G4ThreeVector(0.*cm,0.*cm,
													crystalHalfLength+2.*plexiGlasWindowHalfLength+PMTWindowHalfLength);

	G4VPhysicalVolume* physiPMTWindow = new G4PVPlacement(0,positionPMTWindow,
														  "PMTWindow",logicPMTWindow,
														  physiWorld,false,0);

	// Photocathode

	G4Tubs* solidCathode = new G4Tubs("Cathode",0.*cm,cathodeRadius,
									  cathodeHalfLength,startPhi,deltaPhi);

	G4LogicalVolume* logicCathode = new G4LogicalVolume(solidCathode,
														Bialkali,"Cathode");

	G4ThreeVector positionCathode = G4ThreeVector(0.*cm,0.*cm,
												  crystalHalfLength+2.*plexiGlasWindowHalfLength+2.*PMTWindowHalfLength
												  +cathodeHalfLength);

	G4VPhysicalVolume* physiCathode = new G4PVPlacement(0,positionCathode,
														"Cathode",logicCathode,
														physiWorld,false,0);


	//------------------------------------------------------
	// Surfaces and boundary processes
	//------------------------------------------------------


	// Reflector - scintillator surface

	G4OpticalSurface* OpCryRefSurface =
	new G4OpticalSurface("CryRefSurface");

	OpCryRefSurface->SetType(dielectric_metal);
	OpCryRefSurface->SetModel(glisur);
	OpCryRefSurface->SetFinish(polished);

	G4LogicalBorderSurface* CryRefSurface =
    new G4LogicalBorderSurface("CryRefSurface",physiCrystal,
							   physiReflector,OpCryRefSurface);

	// scintillator - scintillatorPlexiGlass surface

	G4OpticalSurface* OpCryPlexiSurface =
	new G4OpticalSurface("CryPlexiSurface");

	OpCryPlexiSurface->SetType(dielectric_dielectric);
	OpCryPlexiSurface->SetModel(glisur);
	OpCryPlexiSurface->SetFinish(polished);

	G4LogicalBorderSurface* CryPlexiSurface =
    new G4LogicalBorderSurface("CryPlexiSurface",physiCrystal,
							   physiCoatingPlexi,OpCryPlexiSurface);

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


	// scintillatorPlexiGlass - PMT window surface

	G4OpticalSurface* OpPlexiPMTWinSurface =
	new G4OpticalSurface("CryPlexiPMTWinSurface");

	OpPlexiPMTWinSurface->SetType(dielectric_dielectric);
	OpPlexiPMTWinSurface->SetModel(glisur);
	OpPlexiPMTWinSurface->SetFinish(polished);

	G4LogicalBorderSurface* CryPlexiPMTWinSurface =
    new G4LogicalBorderSurface("CryPlexiPMTWinSurface",physiCoatingPlexi,physiPMTWindow,
							   OpPlexiPMTWinSurface);

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

	// White color for Detector
	G4VisAttributes* VisAtt1= new G4VisAttributes(G4Colour(1.0,1.0,1.0)); //White
	logicCrystal->SetVisAttributes(VisAtt1);

	// Red color for Shielding
	G4VisAttributes* VisAtt2 = new G4VisAttributes(G4Colour(1.0,0.0,0.0)); //red
	logicShielding->SetVisAttributes(VisAtt2);

	// Gray color for Coating (Aluminium)
	G4VisAttributes* VisAtt3= new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
	logicCoating->SetVisAttributes(VisAtt3);

	// Brown color for CoatingPlexiGlass
	G4VisAttributes* VisAtt4= new G4VisAttributes(G4Colour(0.45,0.25,0.0));
	logicCoatingPlexi->SetVisAttributes(VisAtt4);

	// Yellow color for reflector
	G4VisAttributes* VisAtt7= new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	logicReflector->SetVisAttributes(VisAtt7);

	// Blue color for PMT window
	G4VisAttributes* VisAtt5= new G4VisAttributes(G4Colour(0.0,0.0,1.0));
	logicPMTWindow->SetVisAttributes(VisAtt5);

	// White color for the absorber photocathode
	G4VisAttributes* VisAtt6= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
	logicCathode->SetVisAttributes(VisAtt6);








	//
	// always return the physical World
	//

  return physiWorld;
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
