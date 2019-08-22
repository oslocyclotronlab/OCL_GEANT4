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
#include "OCLLaBr3Messenger.hh"
#include "OCLMaterials.hh"
#include "OCLParameters.hh"

#include "G4ExceptionSeverity.hh"
#include "G4NistManager.hh"
#include "G4RunManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

OCLLaBr3::OCLLaBr3(OCLLaBr3Messenger* fMessenger)
{
	GetMaterials();

	reflectorThickness = fMessenger->GetReflectorThickness();
	coatingAlThickness = fMessenger->GetCoatingAlThickness();
	coatingAlThicknessFront = fMessenger->GetCoatingAlThicknessFront();
	shieldingHalfThicknessLid = fMessenger->GetShieldingHalfThicknessLid();

	CalculateGeometryParameters();
	CreateSolids();
	SetOpticalProperties();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Destructor
OCLLaBr3::~OCLLaBr3() {
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void OCLLaBr3::SetPosition(G4ThreeVector thisPos) {
  fTranslatePos = thisPos*mm;
  G4cout << " ----> A OCLLaBr3 will be placed at distance: " << fTranslatePos.getR()/cm << " cm" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void OCLLaBr3::SetRotation(G4RotationMatrix thisRot) {
	fRotation = thisRot;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

	//some default Clover detector parameters
	// --> Parameter file

void OCLLaBr3::GetMaterials() {
	//----------------------------------------------------
	// Material definitions
	//----------------------------------------------------

	G4double a, z;                    //a=mass of a mole;
	G4double density;                 //z=mean number of protons;

	G4int ncomponents, natoms;
	G4double abundance, fractionmass;

  // // Get materials
  OCLMaterials* fMat = OCLMaterials::GetInstance();

  // //Lead
	Al2O3      = fMat->GetMaterial("G4_ALUMINUM_OXIDE");
	Aluminium  = fMat->GetMaterial("G4_Al");
	B2O3       = fMat->GetMaterial("G4_BORON_OXIDE");
	K2O        = fMat->GetMaterial("G4_POTASSIUM_OXIDE");
	Na2O       = fMat->GetMaterial("G4_SODIUM_MONOXIDE");
	PlexiGlass = fMat->GetMaterial("G4_PLEXIGLASS");
	SiO2       = fMat->GetMaterial("G4_SILICON_DIOXIDE");

	LaBr3_Ce = fMat->GetMaterial("LaBr3_Ce");

	// MgO reflector
	MgO = fMat->GetMaterial("MgO");

	// Air (non-STP)
	Air = fMat->GetMaterial("Air");

	// Borosilicate
	Borosilicate = fMat->GetMaterial("Borosilicate glass");
	Bialkali = fMat->GetMaterial("Bialkali");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void  OCLLaBr3::CreateSolids()
{
	//
	// Detector Geometry
	//

    solidPMTandAir = new G4Tubs("PMTandAir",
                                    0. * mm,
                                    coatingOuterR,
                                    PMTandAirHalfLength,
                                    startPhi,
                                    deltaPhi);

    logicPMTandAir = new G4LogicalVolume(solidPMTandAir, Air, "PMTandAir");

  	//
    // Shielding
  	//

	// concial section

	solidShieldingConical = new G4Cons("ShieldingConnical",
										shieldingInnerR,   // inner radius = 0 because used as mother volume
										shieldingConeOuterRFront,
										shieldingInnerR,   // inner radius = 0 because used as mother volume
										shieldingConeOuterRBack,
										shieldingConeHalfLength,
										startPhi,
										deltaPhi);
    for(int i=0; i<3; i++){
        solidShieldTubs[i] = new G4Tubs(("ShieldingTube_"+std::to_string(i)).c_str(),
                                    shieldingInnerR,
                                    rOuters[i], // we assume that there is no gap between the coating/crystal and shielding
                                    dShieldTubs[i]/2.,
                                    startPhi,
                                    deltaPhi);
    }

    unionShielding = new G4MultiUnion("Shield_Union");

    G4RotationMatrix rotm = G4RotationMatrix();
    G4ThreeVector positionC = G4ThreeVector(0.,0.,shieldingConeHalfLength                                          - detectorHalfinclPMT);
    G4ThreeVector position0 = positionC + G4ThreeVector(0.,0.,dShieldTubs[0]/2.+shieldingConeHalfLength);
    G4ThreeVector position1 = position0 + G4ThreeVector(0.,0.,dShieldTubs[1]/2.+dShieldTubs[0]/2.);
    G4ThreeVector position2 = position1 + G4ThreeVector(0.,0.,dShieldTubs[2]/2.+dShieldTubs[1]/2.);
    G4Transform3D trC = G4Transform3D(rotm,positionC);
    G4Transform3D tr0 = G4Transform3D(rotm,position0);
    G4Transform3D tr1 = G4Transform3D(rotm,position1);
    G4Transform3D tr2 = G4Transform3D(rotm,position2);
    unionShielding->AddNode(*solidShieldingConical,trC);
    unionShielding->AddNode(*solidShieldTubs[0],tr0);
    unionShielding->AddNode(*solidShieldTubs[1],tr1);
    unionShielding->AddNode(*solidShieldTubs[2],tr2);

    // Finally close the structure
    //
    unionShielding->Voxelize();

    logicShielding = new G4LogicalVolume(unionShielding, Aluminium, "Shielding_LV");

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
	  				coatingOuterR-coatingAlThickness,
	  				coatingHalfLength - 0.5 * coatingAlThicknessFront, // in order to get the coatingPlasticHalfLength
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


    // PMT (simplified)
    G4double zPlane[3] = {0-PMTHalfLength, PMTMidtZ - PMTHalfLength, PMTEndZ - PMTHalfLength};
    G4double rInners[3] = {0};
    G4double rOuters[3] = {PMTStartRadius, PMTMidtRadius, PMTEndRadius};
    solidPMT =
        new G4Polycone("PMT",
                       startPhi,
                       deltaPhi,
                       3, // nZPlanes
                       zPlane,
                       rInners,
                       rOuters);

    logicPMT = new G4LogicalVolume(solidPMT, Aluminium, "PMT");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void OCLLaBr3::SetOpticalProperties() {
	// at the moment taken from the Scintiallator example
	G4int nEntries = 2;

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

	// Air

	G4double vacRefractionIndex[nEntries] = {1.0,1.0};
	vacMPT = new G4MaterialPropertiesTable();
	vacMPT->AddProperty("RINDEX",
						PhotonEnergy,
						vacRefractionIndex,
						nEntries);
	Air->SetMaterialPropertiesTable(vacMPT);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void OCLLaBr3::Placement(G4int copyNo, G4VPhysicalVolume* physiMother, G4bool checkOverlaps)
{

	// Create the transformation vector to mother volume
	// tranlation: distance middle + halflength = distance to face
	fTranslatePos.setR(fTranslatePos.getR()
	                   + detectorHalfinclPMT);
	G4Transform3D transDetector = G4Transform3D(fRotation, fTranslatePos);

	G4int copyNoSub = 0; // copy number for the sub elements (could also be copyNo)

 	//
	// Detector Geometry
	//

  	//
    // Shielding

	physiShield =
		new G4PVPlacement(transDetector, 	// Transformation (Rot&Transl)
						  "OCLDetector",    // its name
						  logicShielding, 	// its logical volume
						  physiMother, 	// its physical mother volume
						  false, 				// unknown "pMany"; def: false
						  copyNo, 			// copy number
						  checkOverlaps);		// checkOverlaps


	//
	// Detector Crystal (incl. Reflector & Coating)
	//

	// Coating
	// Coating: Aluminum part
    positionCoating = G4ThreeVector(0.*cm,0.*cm,
                                    2*shieldingHalfThicknessLid + coatingHalfLength
                                    - detectorHalfinclPMT); // because of the shift in the coordinate system of the shielding

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
	positionCoatingPlexi = G4ThreeVector(0.*cm, 0.*cm, 0.5*coatingAlThicknessFront); // because of the shift in the coordinate system of the coating

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

    positionPlexiWindow = positionCoating
                              + G4ThreeVector(0.*cm,
                                              0.*cm,
                                              coatingHalfLength
                                              + plexiGlasWindowHalfLength);


	physiPlexiWindow =
	 new G4PVPlacement(0,						// Rotation
					   positionPlexiWindow,		// Transformation (Rot&Transl)
					   "PlexiGlasWindow_on_Crystal",		// its logical volume
					   logicPlexiWindow,		// its name
					   physiShield,		        // its physical mother volume
					   false,					// unknown "pMany"; def: false
					   copyNoSub,				// copy number
					   checkOverlaps);			// checkOverlaps
    //
    // PMT and Air
    //
    positionPMTandAir = positionPlexiWindow + G4ThreeVector(0.*cm, 0.*cm,
                                                            plexiGlasWindowHalfLength
                                                            + PMTandAirHalfLength);

    physiPMTandAir =
        new G4PVPlacement(0,                    // Rotation
                          positionPMTandAir,    // Transformation (Rot&Transl)
                          "PMTandAir",    // its name
                          logicPMTandAir,   // its logical volume
                          physiShield,  // its physical mother volume
                          false,                // unknown "pMany"; def: false
                          copyNoSub,           // copy number
                          checkOverlaps);       // checkOverlaps


	//
	// PMT
	//


    // PMT window

	positionPMTWindow = G4ThreeVector(0.*cm,
							          0.*cm,
							        	plexiGlasWindowHalfLength
							        	- PMTandAirHalfLength);

	physiPMTWindow =
		new G4PVPlacement(0,					// Rotation
						  positionPMTWindow,	// Transformation (Rot&Transl)
						  "PMTWindow",			// its logical volume
						  logicPMTWindow,		// its name
						  physiPMTandAir,		// its physical mother volume
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
						  physiPMTandAir,		// its physical mother volume
						  false,				// unknown "pMany"; def: false
						  copyNoSub,			// copy number
						  checkOverlaps);		// checkOverlaps

    // PMT (simplified)

    positionPMT = positionCathode
                      + G4ThreeVector(0.*cm,
                                      0.*cm,
                                      cathodeHalfLength + PMTHalfLength);

    physiPMT =
        new G4PVPlacement(0,                    // Rotation
                          positionPMT,      // Transformation (Rot&Transl)
                          "PMT",            // its logical volume
                          logicPMT,         // its name
                          physiPMTandAir,       // its physical mother volume
                          false,                // unknown "pMany"; def: false
                          copyNoSub,            // copy number
                          checkOverlaps);       // checkOverlaps



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

	// White color for Crystal
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

    // Yellow color for PMTandAir
    // VisAtt7= new G4VisAttributes(G4Colour(1.0,1.0,0.0));
    logicPMTandAir->SetVisAttributes(VisAtt7);

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void OCLLaBr3::CalculateGeometryParameters(){

	// +-----------------------------------+
	// |                Aluminum coating | |
	// | +---------------------------------+     XXXXX
	// | |              Plexi Glass coat.| |         X
	// | | +-----------------------------+ |         X
	// | | |            Reflector        | |         X
	// | | | +---------------------------+ |         X
	// | | | |          Crystal          | |         X
	// | | | |                           | |         XX = Detector
	// | | | +---------------------------+ |         X
	// | | |                             | |         X
	// | | +-----------------------------+ |         X
	// | |                               | |         X
	// | +---------------------------------+     XXXXX
	// |                                 | |
	// ++----------------------------------+
	//  ^                                 ^
	//  |                                 |
	//  +Coating Front                    Plexiglass-Window on crystal


	//      Lid
	//      +
	//      |                  ||
	//      v   +--------------++----------+
	//        +       Shielding            |
	//      +++------------+P+-	   	       |
	//      |||  Detector  |M    PMTandAir |
	//      +++------------+T+-            |
	//        +                            |
	//          +--------------+-----------+
	//         ^               ||
	//         |
	//         +
	//         Cone


	// For all circular objects
	startPhi = 0.*deg;
	deltaPhi = 360.*deg;

	// other defaults -> Now in Messenger class
	// reflectorThickness = 1.*mm; // assumption: 1 mm thick reflector on the front side
	// coatingAlThickness = 2.*mm; // thickness as in the radius part
	// coatingAlThicknessFront = 1.*mm; // we assume a smaller thickness at the front of the detector
	// shieldingHalfThicknessLid = 2.*mm/2.;

	crystalOuterR = 8.89*cm/2.; // 3.5 in in cm
	crystalInnerR = 0.0*mm;
	crystalHalfLength = 203.2*mm/2.; // 8 in in cm

	coatingOuterR = 100.*mm/2. ;
	plexiGlasWindowHalfLength = 1.*mm/2.;


	//
	// Detector & Shielding
	//
	reflectorHalfLength = crystalHalfLength + 0.5 * reflectorThickness; // assumption: backside doesn't have a reflector
	//ReflectorInnerR = crystalOuterR;
	reflectorInnerR = 0.*mm;
	reflectorOuterR = crystalOuterR + reflectorThickness;

	// in between reflector and coating, there will be some plastic
	// default assumption: 2.55 mm plexiglas coating around the reflector before the aluminium
	coatingPlasticThickness = coatingOuterR - coatingAlThickness -reflectorThickness - crystalOuterR;
	coatingHalfLength = reflectorHalfLength + 0.5 * coatingAlThicknessFront + 0.5 * coatingPlasticThickness; // backside doesn't have an (Aluminium) coating


    //in the front, the shielding tube diameter is reduces. It's later modeled by a conical section
    //shieldingConeInnerRFront = coatingOuterR;
    shieldingConeOuterRFront = coatingOuterR + 3.*mm;
    G4double deltaShieldingHalfThicknessLid = shieldingHalfThicknessLid - defaultshieldingHalfThicknessLid;
    shieldingConeHalfLength = 27./2.*mm + deltaShieldingHalfThicknessLid ;// in the front, the tube (subtracting default lid)


    G4double shieldingThickness1 = 5.*mm;         // thickness of the 1. tube
    G4double shieldingThickness2 = 30*mm;
    G4double shieldingThickness3 = 10*mm;

    G4double colimatorLength = 70.*mm;
    G4double shieldingLength1 = (110.74*mm + 114.26*mm) - colimatorLength + 2*deltaShieldingHalfThicknessLid;
    G4double shieldingLength2 = 49*mm;
    G4double shieldingLength3 = 415*mm + 2*deltaShieldingHalfThicknessLid
                       -( shieldingLength1 + shieldingLength2);

	shieldingInnerR = 0*mm; 			// as we use it as a mother volume
    G4double shieldingOuterR1 = coatingOuterR + shieldingThickness1;
    G4double shieldingOuterR2 = coatingOuterR + shieldingThickness2;
	G4double shieldingOuterR3 = coatingOuterR + shieldingThickness3;

    dShieldTubs[0] = shieldingLength1;
    dShieldTubs[1] = shieldingLength2;
    dShieldTubs[2] = shieldingLength3;

    rOuters[0] = shieldingOuterR1;
    rOuters[1] = shieldingOuterR2;
    rOuters[2] = shieldingOuterR3;

	//shieldingConeInnerRBack = shieldingConeInnerRFront;
	shieldingConeOuterRBack = coatingOuterR + shieldingThickness1;

	G4cout << "coatingHalfLength " << coatingHalfLength/mm << " mm" << G4endl;
	// without conical Section and Lid
	//  we assume no coating at the back side
	// shieldingHalfLength = coatingHalfLength - shieldingConeHalfLength; // - default value

	plexiGlasWindowOuterR = coatingOuterR; // currently we just assume a flat window on the top.


	//
	// PMT
	//

	PMTWindowHalfLength = 1.0*mm;
	PMTWindowRadius = 85*0.5*mm;

	cathodeHalfLength = 0.005*mm;
	cathodeRadius =85*0.5*mm;

    PMTStartRadius = cathodeRadius;
    PMTMidtRadius = 51.5*mm/2.;
    PMTEndRadius = 56.5*mm/2.;

    PMTMidtZ = 7.*cm;
    PMTEndZ = PMTMidtZ+3.*cm;
    PMTHalfLength = PMTEndZ/2.;

	//
	// Whole detector incl. PMT (-> Logical unit)
	//

    detectorHalfinclPMT = shieldingConeHalfLength
                          + shieldingLength1/2. +shieldingLength2/2. + shieldingLength3/2.;

    // PMTandAir

    PMTandAirHalfLength = detectorHalfinclPMT - shieldingHalfThicknessLid - coatingHalfLength
                          - plexiGlasWindowHalfLength;


  if (coatingPlasticThickness<0) {
  	G4Exception("OCLLaBr3::CalculateGeometryParameters","coatingPlasticThickness",
		            FatalErrorInArgument, "Value cannot be negative.\n"
		            "Check choice of geometry parameters via marco");}
  else { 	G4cout << "coatingPlasticThickness: " << coatingPlasticThickness/mm << " mm" << G4endl; }

}
