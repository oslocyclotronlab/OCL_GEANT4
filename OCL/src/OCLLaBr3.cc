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


#include "OCLLaBr3.hh"
#include "OCLParameters.hh"

#include "G4NistManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

OCLLaBr3::OCLLaBr3()
{
	//some default Clover detector parameters
	// --> Parameter file


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

	Br = new G4Element("Bromium",    "Br",   z=35.,  a=79.904*g/mole);
	La = new G4Element("Lanthanum",  "La",   z=57.,  a=138.90547*g/mole);
	Ce = new G4Element("Cerium",     "Cl",   z=58.,  a=140.116*g/mole);

	// add more elements from NIST database
	O  = man->FindOrBuildElement("O");
	K  = man->FindOrBuildElement("K");
	Sb = man->FindOrBuildElement("Sb");
	Cs = man->FindOrBuildElement("Cs");
	Mg = man->FindOrBuildElement("Mg");

	//
	// define materials from elements.
	//

	Aluminium  = man->FindOrBuildMaterial("G4_Al");
	PlexiGlass = man->FindOrBuildMaterial("G4_PLEXIGLASS");
	SiO2       = man->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
	Na2O       = man->FindOrBuildMaterial("G4_SODIUM_MONOXIDE");
	K2O        = man->FindOrBuildMaterial("G4_POTASSIUM_OXIDE");
	Al2O3      = man->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
	B2O3       = man->FindOrBuildMaterial("G4_BORON_OXIDE");

	//LaBr3
	LaBr3 =   new G4Material("LaBr3", density = 5.07*g/cm3, ncomponents=2);
	LaBr3->AddElement(La, natoms=1);
	LaBr3->AddElement(Br, natoms=3);

	//LaBr3_Ce
	//with 5% dopping, see technical note "BrilLanCe Scintillators Performance Summary"
	LaBr3_Ce = new G4Material("LaBr3_Ce", density = 5.08*g/cm3, ncomponents=2);
	LaBr3_Ce->AddMaterial(LaBr3,  fractionmass=95.*perCent);
	LaBr3_Ce->AddElement(Ce,      fractionmass=5.*perCent);

	// MgO reflector
	density = 2.0*g/cm3;
	MgO = new G4Material("MgO", density, ncomponents=2);
	MgO->AddElement(Mg, natoms=1);
	MgO->AddElement(O, natoms=1);

	// vacuum (non-STP)

	vacuum = new G4Material("Vacuum",       //name as String
							1,		                    //atomic number (use 1 for Hydrogen)
	                		1.008*g/mole, 	            //molar mass (use 1.008*g/mole for Hydoren)
							1.e-25*g/cm3,  	            //density
							kStateGas,		            //kStateGas - the material is gas (see G4State)
	                		2.73*kelvin,	            //Temperature
							1.e-25*g/cm3);	            //pressure


	// Steel as non-NIST material
	// elFe = G4NistManager::Instance()->FindOrBuildElement("Fe");
	// elNi = G4NistManager::Instance()->FindOrBuildElement("Ni");
	// elCr = G4NistManager::Instance()->FindOrBuildElement("Cr");
	// iron = new G4Material("StainlessSteel", 7.80 * g/cm3, 3 /* components */);
	// iron -> AddElement(elFe, 70 * perCent);
	// iron -> AddElement(elCr, 18 * perCent);
	// iron -> AddElement(elNi, 12 * perCent);

	// PMT-materials

	// Borosilicate
	Borosilicate = new G4Material("Borosilicate glass", density= 2.23*g/cm3, ncomponents=5);
	Borosilicate->AddMaterial(SiO2,   fractionmass=80.6 * perCent);
	Borosilicate->AddMaterial(B2O3,  fractionmass=13.0 * perCent);
	Borosilicate->AddMaterial(Na2O,  fractionmass=2.   * perCent); // 1/2 of wt% for (Na20+K20)
	Borosilicate->AddMaterial(K2O,   fractionmass=2.   * perCent); // 1/2 of wt% for (Na20+K20)
	Borosilicate->AddMaterial(Al2O3, fractionmass=2.31  * perCent);

	// Bialkali
	// (Bialkali KCsSb,  Density=?, Thickness=?)?
	Bialkali = new G4Material("Bialkali", density= 2.*g/cm3, ncomponents=3);
	Bialkali->AddElement(K,  natoms=2);
	Bialkali->AddElement(Cs, natoms=1);
	Bialkali->AddElement(Sb, natoms=1);


	//------------------------------------------------------
	// Optical properties
	//------------------------------------------------------


	// at the moment taken from the Scintiallator example
		const G4int nEntries = 2;

	// 1eV -> 1.2399 µm; 7eV -> 0.1771µm // TODO more detailed; adopt all of them
	G4double PhotonEnergy[nEntries] = {1.0*eV,7.0*eV}; 

	// MgO reflector

	G4double MgORefractionIndex[nEntries] = {1.0,1.0};
	G4double MgOAbsorptionLength[nEntries] = {1.0E-9*m,1.0E-9*m};
	MgOMPT = new G4MaterialPropertiesTable();
	MgOMPT->AddProperty("RINDEX",
						PhotonEnergy,
						MgORefractionIndex,
						nEntries);
	MgOMPT->AddProperty("ABSLENGTH",
						PhotonEnergy,
						MgOAbsorptionLength,
						nEntries);
	MgO->SetMaterialPropertiesTable(MgOMPT);


	// LaBr3

	G4double LaBr3RefractionIndex[nEntries] = {1.9,1.9};
	G4double LaBr3AbsorptionLength[nEntries] = {50.*cm,50.*cm};

	LaBr3MPT = new G4MaterialPropertiesTable();
	LaBr3MPT->AddProperty("RINDEX",
						  PhotonEnergy,
						  LaBr3RefractionIndex,
						  nEntries);
	LaBr3MPT->AddProperty("ABSLENGTH",
						  PhotonEnergy,
						  LaBr3AbsorptionLength,
						  nEntries);

	G4double ScintEnergy[nEntries] = {3.26*eV,3.44*eV};
	G4double ScintFast[nEntries] = {1.0,1.0};

	LaBr3MPT->AddProperty("FASTCOMPONENT",
						  ScintEnergy,
						  ScintFast,
						  nEntries);

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

	PlexiGlasMPT = new G4MaterialPropertiesTable();
	MgOMPT->AddProperty("RINDEX",
						PhotonEnergy,
						MgORefractionIndex,
						nEntries);
	MgOMPT->AddProperty("ABSLENGTH",
						PhotonEnergy,
						MgOAbsorptionLength,
						nEntries);

	PlexiGlass->SetMaterialPropertiesTable(PlexiGlasMPT);

	// Quartz

	G4double QuartzRefractionIndex[nEntries] = {1.47,1.47};
	G4double QuartzAbsorptionLength[nEntries] = {3.0*cm,3.0*cm};
	QuartzMPT = new G4MaterialPropertiesTable();
	QuartzMPT->AddProperty("RINDEX",
						   PhotonEnergy,
						   QuartzRefractionIndex,
						   nEntries);
	QuartzMPT->AddProperty("ABSLENGTH",
						   PhotonEnergy,
						   QuartzAbsorptionLength,
						   nEntries);
	SiO2->SetMaterialPropertiesTable(QuartzMPT);

	// K2CsSb (Bialcali Photocathode)

	G4double K2CsSbRefractionIndex[nEntries] = {1.47,1.47};
	G4double K2CsSbAbsorptionLength[nEntries] = {1.0E-9*m,1.0E-9*m};
	K2CsSbMPT = new G4MaterialPropertiesTable();
	K2CsSbMPT->AddProperty("RINDEX",
						   PhotonEnergy,
						   K2CsSbRefractionIndex,
						   nEntries);
	K2CsSbMPT->AddProperty("ABSLENGTH",
		                   PhotonEnergy,
		                   K2CsSbAbsorptionLength,
						   nEntries);
	Bialkali->SetMaterialPropertiesTable(K2CsSbMPT);

	// Vacuum

	G4double vacRefractionIndex[nEntries] = {1.0,1.0};
	vacMPT = new G4MaterialPropertiesTable();
	vacMPT->AddProperty("RINDEX",
						PhotonEnergy,
						vacRefractionIndex,
						nEntries);
	vacuum->SetMaterialPropertiesTable(vacMPT);
	// END of OPTICAL PROPERTIES
	////////////////////////////////////////////////////////

	//
	// Create the solids.....
	//
	CreateSolids();

}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Destructor
OCLLaBr3::~OCLLaBr3() {}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void OCLLaBr3::SetPosition(G4ThreeVector thisPos) {
  translatePos = thisPos*mm;
  G4cout << " ----> A OCLLaBr3 will be placed at distance: " << translatePos/mm << " mm" << G4endl;
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void OCLLaBr3::SetRotation(G4RotationMatrix thisRot) { 
	rotation = thisRot; 
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void  OCLLaBr3::CreateSolids()
{
	//
	// Detector Geometry
	//

	// Create a Logical Volume that contains the whole detector unit
	// -> need to specify placement only one logic volume in the world volume later
	solidOCLDetector = new G4Tubs("OCLDetector",
	  								0. * mm, // inner radius = 0 because used as mother volume
	  								shieldingOuterR,
	  								detectorHalfinclPMT,
	  								startPhi,
	  								deltaPhi);

    logicOCLDetector = new G4LogicalVolume(solidOCLDetector, vacuum, "OCLDetector");

  	//
    // Shielding
  	//

  	// Main tube

  	solidShieldingMain = new G4Tubs("ShieldingMainTube",
									 shieldingInnerR,
									 shieldingOuterR,
									 shieldingHalfLength,
									 startPhi,
									 deltaPhi);

	// concial section

	solidShieldingConical = new G4Cons("ShieldingConnical",
										shieldingInnerR,   // inner radius = 0 because used as mother volume
										shieldingConeOuterRFront,
										shieldingInnerR,   // inner radius = 0 because used as mother volume
										shieldingConeOuterRBack,
										shieldingConeHalfLength,
										startPhi,
										deltaPhi);

	// lid

  	solidShieldingLid = new G4Tubs("ShieldingLid",
  									shieldingInnerR,
  									shieldingConeOuterRFront, // we assume that there is no gap between the coating/crystal and shielding
  									shieldingHalfThicknessLid,
  									startPhi,
  									deltaPhi);

  	//
	// add the shielding parts together

  	// The origin and the coordinates of the combined solid are the same as those of
  	// the first solid.

	translationUnion1 = G4ThreeVector(0, 0, - (shieldingConeHalfLength + shieldingHalfThicknessLid) );
    unionShielding1 = 
    new  G4UnionSolid("unionShielding1",
						solidShieldingConical,  // 1st object
						solidShieldingLid,	   	// 2nd object
                    	0,					  	// no Rotation
                    	translationUnion1);   	// translation of the 2nd object

	translationUnion2 = G4ThreeVector(0, 0, -(shieldingHalfLength + shieldingConeHalfLength));
    unionShielding2 = 
    	new  G4UnionSolid ("unionShielding2",
							solidShieldingMain,	// 1nd object
							unionShielding1,	// 2st object
                    		0,					// no Rotation
                    		translationUnion2); // translation of the 2nd object

    logicShielding = new G4LogicalVolume(unionShielding2, Aluminium, "ShieldingLid");

	//
	// Detector Crystal (incl. Reflector & Coating)
	//

	// Coating

	// Coating: Aluminum part

	solidCoating = 
		new G4Tubs("Coating",
	  				0. * mm, // inner radius = 0 because used as mother volume
	  				coatingOuterR,
	  				coatingHalfLength,
	  				startPhi,
	  				deltaPhi);

	logicCoating = new G4LogicalVolume(solidCoating, Aluminium, "Coating");

	// Coating: PlexiGlass part
	solidCoatingPlexi = 
		new G4Tubs("CoatingPlexiGlas",
	  				0. * mm, // inner radius = 0 because used as mother volume
	  				coatingOuterR-coatingThickness,
	  				coatingHalfLength - 0.5 * coatingThicknessFront, // in order to get the coatingPlasticHalfLength
	  				startPhi,
	  				deltaPhi);

	logicCoatingPlexi = new G4LogicalVolume(solidCoatingPlexi, PlexiGlass, "CoatingPlexiGlas");

	// Reflector

	// Assumption: We (currently) don't know whether we really have a Reflector/Material properties
	//             the specifications are taken from the Scintillation example

	solidReflector = 
		new G4Tubs("Reflector",
					reflectorInnerR, 
					reflectorOuterR, 
					reflectorHalfLength,
					startPhi,
					deltaPhi);

	logicReflector = new G4LogicalVolume(solidReflector, MgO, "Reflector");
	
	// Crystal

  	solidCrystal = 
  		new G4Tubs("Crystal",
  					crystalInnerR,
  					crystalOuterR,
  					crystalHalfLength,
  					startPhi,
  					deltaPhi);

    logicCrystal = new G4LogicalVolume(solidCrystal, LaBr3_Ce, "Crystal");


	// Plexiglas Window on Detector
	solidPlexiWindow = 
		new G4Tubs("PlexiGlasWindow",
	  				0. * mm, // inner radius = 0
	  				plexiGlasWindowOuterR,
	  				plexiGlasWindowHalfLength, // in order to get the coatingPlasticHalfLength
	  				startPhi,
	  				deltaPhi);

	logicPlexiWindow = new G4LogicalVolume(solidPlexiWindow, PlexiGlass, "PlexiGlasWindow");

	//
	// PMT
	//

    // PMT window

	solidPMTWindow = 
		new G4Tubs("PMTWindow",
					0.*cm,
					PMTWindowRadius,
					PMTWindowHalfLength,
					startPhi,
					deltaPhi);

	logicPMTWindow = new G4LogicalVolume(solidPMTWindow, SiO2, "PMTWindow");



	// Photocathode

	solidCathode = 
		new G4Tubs("Cathode",
					0.*cm,
					cathodeRadius,
					cathodeHalfLength,
					startPhi,
					deltaPhi);

	logicCathode = new G4LogicalVolume(solidCathode, Bialkali, "Cathode");


}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void OCLLaBr3::Placement(G4int copyNo, G4VPhysicalVolume* physiMother, G4bool checkOverlaps)
{

	// Create the transformation vector to mother volume
	G4Transform3D transDetector = G4Transform3D(rotation,translatePos);

	G4int copyNoSub = 0; // copy number for the sub elements (could also be copyNo)

 	//
	// Detector Geometry
	//

	physiOCLDetector = 
		new G4PVPlacement(transDetector, 		// Transformation (Rot&Transl)
							"OCLDetector",		// its name
							logicOCLDetector, 	// its logical volume
							physiMother, 		// its physical mother volume
							false, 				// unknown "pMany"; def: false
							copyNo , 			// copy number
							checkOverlaps);		// checkOverlaps


  	//
    // Shielding

	// lid
    positionShielding = G4ThreeVector(0.*cm,
    								  0.*cm,  
    								  (2.*shieldingHalfThicknessLid 
    								  	+ 2.*shieldingConeHalfLength 
    								  	+ shieldingHalfLength)
			                          - detectorHalfinclPMT
    												);  	
   																								   
	physiShield = 
		new G4PVPlacement(0, 					// Rotation
							positionShielding, 	// Transformation (Rot&Transl)
							"Shielding", 		// its logical volume
							logicShielding, 	// its name
							physiOCLDetector, 	// its physical mother volume
							false, 				// unknown "pMany"; def: false
							copyNoSub, 			// copy number
							checkOverlaps);		// checkOverlaps


	//
	// Detector Crystal (incl. Reflector & Coating)
	//

	// Coating
	// Coating: Aluminum part
	positionCoating = G4ThreeVector(0.*cm,0.*cm,-shieldingConeHalfLength); // because of the shift in the coordinate system of the shielding
																		   // (center != origin)
	physiCoating = 
		new G4PVPlacement(0,					// Rotation
							positionCoating,	// Transformation (Rot&Transl)
							"Coating",			// its logical volume
							logicCoating,		// its name
							physiShield,		// its physical mother volume
							false,				// unknown "pMany"; def: false
							copyNoSub,			// copy number
							checkOverlaps);		// checkOverlaps
	
	// Coating: PlexiGlass part
	positionCoatingPlexi = G4ThreeVector(0.*cm, 0.*cm, 0.5*coatingThicknessFront); // because of the shift in the coordinate system of the coating

	physiCoatingPlexi = 
		new G4PVPlacement(0,					// Rotation
						  positionCoatingPlexi,	// Transformation (Rot&Transl)
						  "CoatingPlexiGlas",	// its logical volume		
						  logicCoatingPlexi,	// its name
						  physiCoating,			// its physical mother volume
						  false,				// unknown "pMany"; def: false
						  copyNoSub,			// copy number
						  checkOverlaps);		// checkOverlaps

	// Reflector

	positionReflector = G4ThreeVector(0.*cm,
									  0.*cm,
									  0.5 * coatingPlasticThickness);

	physiReflector = 
		new G4PVPlacement(0,					// Rotation
						  positionReflector,	// Transformation (Rot&Transl)
						  "Reflector",			// its logical volume
						  logicReflector,		// its name
						  physiCoatingPlexi,	// its physical mother volume
						  false,				// unknown "pMany"; def: false
						  copyNoSub,			// copy number
						  checkOverlaps);		// checkOverlaps

	// Crystal

    positionCrystal = G4ThreeVector(0.*cm, 0.*cm, 0.5*reflectorThickness);

	physiCrystal = 
		new G4PVPlacement(0,					// Rotation
						  positionCrystal,	// Transformation (Rot&Transl)
						  "Crystal",			// its logical volume
						  logicCrystal,		// its name
						  physiReflector,	// its physical mother volume
						  false,				// unknown "pMany"; def: false
						  copyNoSub,			// copy number
						  checkOverlaps);	// checkOverlaps

	// Plexiglas Window on Detector

	positionPlexiWindow = positionShielding 
	                      + G4ThreeVector(0.*cm,
										  0.*cm,
										  shieldingHalfLength 
										  + plexiGlasWindowHalfLength);

	physiPlexiWindow = 
	 new G4PVPlacement(0,						// Rotation
					   positionPlexiWindow,		// Transformation (Rot&Transl)
					   "PlexiGlasWindow",		// its logical volume
					   logicPlexiWindow,		// its name
					   physiOCLDetector,		// its physical mother volume
					   false,					// unknown "pMany"; def: false
					   copyNoSub,				// copy number
					   checkOverlaps);			// checkOverlaps

	//
	// PMT
	//


    // PMT window

	positionPMTWindow = positionPlexiWindow 
						+ G4ThreeVector(0.*cm,
							        	0.*cm,
							        	plexiGlasWindowHalfLength 
							        	+ PMTWindowHalfLength);

	physiPMTWindow = 
		new G4PVPlacement(0,					// Rotation
						  positionPMTWindow,	// Transformation (Rot&Transl)
						  "PMTWindow",			// its logical volume
						  logicPMTWindow,		// its name
						  physiOCLDetector,		// its physical mother volume
						  false,				// unknown "pMany"; def: false
						  copyNoSub,			// copy number
						  checkOverlaps);		// checkOverlaps


	// Photocathode

	positionCathode = positionPMTWindow 
					  + G4ThreeVector(0.*cm,
					  				  0.*cm,
					  	              PMTWindowHalfLength + cathodeHalfLength);

	physiCathode = 
		new G4PVPlacement(0,					// Rotation
						  positionCathode,		// Transformation (Rot&Transl)
						  "Cathode",			// its logical volume
						  logicCathode,			// its name
						  physiOCLDetector,		// its physical mother volume
						  false,				// unknown "pMany"; def: false
						  copyNoSub,			// copy number
						  checkOverlaps);		// checkOverlaps
      


	//------------------------------------------------------
	// Surfaces and boundary processes
	//------------------------------------------------------


	// Reflector - scintillator surface

	OpCryRefSurface = 
		new G4OpticalSurface("CryRefSurface");

	OpCryRefSurface->SetType(dielectric_metal);
	OpCryRefSurface->SetModel(glisur);
	OpCryRefSurface->SetFinish(polished);

	CryRefSurface =  
		new G4LogicalBorderSurface("CryRefSurface",
									physiCrystal,
							   		physiReflector,
							   		OpCryRefSurface);

	// scintillator - scintillatorPlexiGlass surface

	OpCryPlexiSurface =	
		new G4OpticalSurface("CryPlexiSurface");

	OpCryPlexiSurface->SetType(dielectric_dielectric);
	OpCryPlexiSurface->SetModel(glisur);
	OpCryPlexiSurface->SetFinish(polished);

	CryPlexiSurface = 
		new G4LogicalBorderSurface("CryPlexiSurface",
									physiCrystal,
							   		physiCoatingPlexi,
							   		OpCryPlexiSurface);

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

	OpPlexiPMTWinSurface = new G4OpticalSurface("CryPlexiPMTWinSurface");

	OpPlexiPMTWinSurface->SetType(dielectric_dielectric);
	OpPlexiPMTWinSurface->SetModel(glisur);
	OpPlexiPMTWinSurface->SetFinish(polished);

	CryPlexiPMTWinSurface = 
		new G4LogicalBorderSurface("CryPlexiPMTWinSurface",
									physiCoatingPlexi,
									physiPMTWindow,
							   		OpPlexiPMTWinSurface);

	// PMT window - photocathode surface

	OpPMTWinCathSurface = new G4OpticalSurface("PMTWinCathSurface");

	OpPMTWinCathSurface->SetType(dielectric_dielectric);
	OpPMTWinCathSurface->SetModel(glisur);
	OpPMTWinCathSurface->SetFinish(polished);

	PMTWinCathSurface = 
		new G4LogicalBorderSurface("CathodeSurface",
		                            physiPMTWindow,
		                            physiCathode,
							   	    OpPMTWinCathSurface);







	//------------------------------------------------------
	// visualization attributes
	//------------------------------------------------------

	// White color for Detector
	VisAtt1= new G4VisAttributes(G4Colour(1.0,1.0,1.0)); //White
	logicCrystal->SetVisAttributes(VisAtt1);

	// Red color for Shielding
	VisAtt2 = new G4VisAttributes(G4Colour(1.0,0.0,0.0)); //red
	logicShielding->SetVisAttributes(VisAtt2);

	// Gray color for Coating (Aluminium)
	VisAtt3= new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
	logicCoating->SetVisAttributes(VisAtt3);

	// Brown color for CoatingPlexiGlass
	VisAtt4= new G4VisAttributes(G4Colour(0.45,0.25,0.0));
	logicCoatingPlexi->SetVisAttributes(VisAtt4);

	// Yellow color for reflector
	VisAtt7= new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	logicReflector->SetVisAttributes(VisAtt7);

	// Blue color for PMT window
	VisAtt5= new G4VisAttributes(G4Colour(0.0,0.0,1.0));
	logicPMTWindow->SetVisAttributes(VisAtt5);

	// White color for the absorber photocathode
	VisAtt6= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
	logicCathode->SetVisAttributes(VisAtt6);


}
