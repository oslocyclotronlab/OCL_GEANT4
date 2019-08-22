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


#include "OCLDetectorConstruction.hh"
#include "OCLParallelWorldFrameOuter.hh"
#include "OCLDetectorMessenger.hh"
#include "OCLMaterials.hh"
#include "OCLLaBr3.hh"
#include "OCLLaBr3Messenger.hh"
#include "OCLCollimator.hh"
#include "OCLFrame.hh"
#include "SiRi.hh"
#include "NIFF.hh"
#include "TChamberAuspuff.hh"
#include "OCLTarget.hh"
#include "OCLTarget_RadSource.hh"

#include "G4SystemOfUnits.hh"
#include "G4Transform3D.hh"
#include "G4PhysicalConstants.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4StateManager.hh"

#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4Colour.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

OCLDetectorConstruction::OCLDetectorConstruction()
:
     solidWorld(0), WorldLog(0), physiWorld(0),
     fworld_sizeXYZ(290*cm),
     fUseCSGOldTargetChamber(false),
     fUseCSGOldTarget(false),
     fUseCSGRadSource(false),
     fUseCSGSiRi(false),
     fUseCSGNiff(false)
{

	for(G4int i=0; i<numberOf_OCLLaBr3; i++){
    fOCLLaBr3_Distance[i] = 16.3*cm;

    if ((i!=0) and (i!=31)){ //beamline
      fOCLLaBr3_presence[i] = true;
    }
    else {fOCLLaBr3_presence[i] = false;}
  }

  fMessenger = new OCLDetectorMessenger(this);
  fMessenger_labr = new OCLLaBr3Messenger(); // to modify the static members
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

OCLDetectorConstruction::~OCLDetectorConstruction()
{
	delete fMessenger;
  delete fMessenger_labr;

	// for(G4int i=0; i<numberOf_OCLLaBr3; i++){
 //  	delete rotmOCLLaBr3[i];
	// }
}

G4VPhysicalVolume* OCLDetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* OCLDetectorConstruction::ConstructVolumes()
{
	// Cleanup old geometry
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // Get materials
  OCLMaterials* fMat = OCLMaterials::GetInstance();
  G4Material* air = fMat->GetMaterial("Air");

	//------------------------------------------------------
	// Detector geometry
	//------------------------------------------------------

	//
	// All parameters have been moved to Parameters.hh
	//

	//
	// World
	//

	solidWorld = 					 //size (defined through half-sizes)
		new G4Box("World",
	               fworld_sizeXYZ/2,
	               fworld_sizeXYZ/2,
	               fworld_sizeXYZ/2);

	WorldLog =
		new G4LogicalVolume(solidWorld,        	//solid defining the World
	                    	air,           	//material of the World
	                    	"World_mass_log");         	//name

	physiWorld =
		new G4PVPlacement(0,                    //specifies rotation: 0 = no rotation
	                  	G4ThreeVector(),     	//at (0,0,0)
	                  	WorldLog,            	//logical volume
										  "World_mass_phys",    //name
										  0,                    //mother  volume
										  false,                //no boolean operation
										  0);                   //copy number

	////////////////////////
	// Positinging
	////////////////////////

	SetPlacementParameters();

	//
	// possitioning
	//

	G4int copynumber;
	bool pSurfChk = false;

	//
	// LaBr3
	//

	for(G4int i=0; i<numberOf_OCLLaBr3; i++){
	if( fOCLLaBr3_presence[i])
		{
		copynumber=i;
		labr3[i] = new OCLLaBr3(fMessenger_labr);
		labr3[i]->SetRotation(rotmOCLLaBr3[i]);
		labr3[i]->SetPosition(positionOCLLaBr3[i]);
		labr3[i]->Placement(copynumber,  physiWorld, pSurfChk);
		}
	}


	///////////

	//
	// Collimator
	//

  // NOTE: This class has not been updated for a while
  //       As we don't use these (anymore), it may be ignored
	// for(G4int i=0; i<numberOf_OCLLaBr3; i++){
	// if( fOCLCollimator_presence[i])
	// 	{
	// 	copynumber=i;
	// 	collimator[i] = new OCLCollimator();
	// 	collimator[i]->SetRotation(rotmOCLLaBr3[i]); // same rotation as Detector
	// 	collimator[i]->SetPosition(positionCollimator[i]);
	// 	// collimator[i]->Placement(copynumber,  physiWorld, pSurfChk);
	// 	}
	// }


	//
	//	Frame (Ball)
	//
	G4RotationMatrix rotmFrame = G4RotationMatrix(); 	// Check implementation: Need to be empty?
	G4ThreeVector 	 positionFrame = G4ThreeVector(); 	//

	OCLFrame* frame;
	frame = new OCLFrame();
	frame->SetRotation(rotmFrame);
	frame->SetPosition(positionFrame);
	frame->Placement(0,  physiWorld, pSurfChk);


	//
	// Old Target Chamber
	//

  if (fUseCSGOldTargetChamber) {
  	TChamberAuspuff* chamber;
  	chamber = new TChamberAuspuff();
  	chamber->Placement(0,  physiWorld, pSurfChk);
  }

	//
	// (Old) Target & Target Holder
	//

  if (fUseCSGOldTarget) {
  	OCLTarget* target;
  	target = new OCLTarget();
    target->Placement(0,  physiWorld, pSurfChk);
  }

  //
  // Radioactive Calibration source (rectangular)
  //

  if (fUseCSGRadSource) {
    if (fUseCSGOldTarget) {
      G4Exception("OCLDetectorConstruction","Fatal error in Argument",
                  FatalErrorInArgument,
                  G4String("fUseCSGRadSource and fUseCSGOldTarget cannot "
                           "be true at the same time").c_str());
    }
  	OCLTarget_RadSource* radSource;
  	radSource = new OCLTarget_RadSource();
  	radSource->Placement(0,  physiWorld, pSurfChk);
  }

	//
	//	SiRi
	//
  if (fUseCSGSiRi) {
  	// G4RotationMatrix rotmSiRi = G4RotationMatrix(); 	// Check implementation: Need to be empty?
  	// G4ThreeVector 	 positionSiRi = G4ThreeVector();	//
  	SiRi* siri;
  	siri = new SiRi();
  	siri->SetAngle(137*deg); // backward
  	// siri->SetAngle(43*deg); // forward

  	// siri->SetRotation(rotmSiRi);
  	// siri->SetPosition(positionSiRi);
  	siri->Placement(0,  physiWorld, pSurfChk);
  }

	//
	//	NIFF
	//

  if (fUseCSGNiff) {
  	// G4RotationMatrix rotmNIFF = G4RotationMatrix(); 	// Check implementation: Need to be empty?
  	// G4ThreeVector 	 positionNIFF = G4ThreeVector();	//
  	NIFF* niff;
  	niff = new NIFF();
  	// niff->SetRotation(rotmNIFF);
  	// niff->SetPosition(positionNIFF);
  	niff->Placement(0,  physiWorld, pSurfChk);
  }


 	//
	// always return the physical World
	//

  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector SpherToCatG4three(G4double r,G4double theta,G4double phi){
	return r * G4ThreeVector( sin(theta) * cos(phi),
					   sin(theta) * sin(phi),
					   cos(theta));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void OCLDetectorConstruction::SetPlacementParameters()
{
  // G4double distColltoDet = 10*mm; // Distance between Collimator and Detector (surface to Surf)

  //Beamline
  // fOCLLaBr3_presence[0]		= false;
  // fOCLCollimator_presence[0]	= false;
  // fOCLLaBr3_Distance[0]		= 16.3*cm;
  OCLLaBr3_theta[00]			= 180.000000*deg;
  OCLLaBr3_phi[0]			= 180.000000*deg;

  //Det. number 1 at OCL:
  // fOCLLaBr3_presence[ 1]		= true;
  // fOCLCollimator_presence[ 1]	= true;
  // fOCLLaBr3_Distance[ 1]		= 16.3*cm;
  OCLLaBr3_theta[ 1]			= 142.622561*deg;
  OCLLaBr3_phi[ 1]			= 359.999847*deg;

  //Det. number 2 at OCL
  // fOCLLaBr3_presence[ 2]		= true;
  // fOCLCollimator_presence[ 2]	= true;
  // fOCLLaBr3_Distance[ 2]		= 16.3*cm;
  OCLLaBr3_theta[ 2]			= 142.622528*deg;
  OCLLaBr3_phi[ 2]			= 71.999866*deg;

  //Det. number 3 at OCL:
  // fOCLLaBr3_presence[ 3]		= true;
  // fOCLCollimator_presence[ 3]	= true;
  // fOCLLaBr3_Distance[ 3]		= 16.3*cm;
  OCLLaBr3_theta[ 3]			= 142.622528*deg;
  OCLLaBr3_phi[ 3]			= 144.000134*deg;

  //Det. number 4 at OCL
  // fOCLLaBr3_presence[ 4]		= true;
  // fOCLCollimator_presence[ 4]	= true;
  // fOCLLaBr3_Distance[ 4]		= 16.3*cm;
  OCLLaBr3_theta[ 4]			= 142.622561*deg;
  OCLLaBr3_phi[ 4]			= 216.000153*deg;

  //Det. number 5 at OCL
  // fOCLLaBr3_presence[ 5]		= true;
  // fOCLCollimator_presence[ 5]	= true;
  // fOCLLaBr3_Distance[ 5]		= 16.3*cm;
  OCLLaBr3_theta[ 5]			= 142.622535*deg;
  OCLLaBr3_phi[ 5]			= 288.000000*deg;

  //Det. number 6
  // fOCLLaBr3_presence[6]		= true;
  // fOCLCollimator_presence[6]	= true;
  // fOCLLaBr3_Distance[6]		= 16.3*cm;
  OCLLaBr3_theta[6]			= 116.564893*deg;
  OCLLaBr3_phi[6]				= 323.999989*deg;

  //Det. number 7
  // fOCLLaBr3_presence[ 7]		= true;
  // fOCLCollimator_presence[ 7]	= true;
  // fOCLLaBr3_Distance[ 7]		= 16.3*cm;
  OCLLaBr3_theta[ 7]			= 100.812208*deg;
  OCLLaBr3_phi[ 7]			= 0.000044*deg;

  //Det. number 8
  // fOCLLaBr3_presence[ 8]		= true;
  // fOCLCollimator_presence[ 8]	= true;
  // fOCLLaBr3_Distance[ 8]		= 16.3*cm;
  OCLLaBr3_theta[ 8]			= 116.564908*deg;
  OCLLaBr3_phi[ 8]			= 35.999955*deg;

  //Det. number 9 at OCL:
  // fOCLLaBr3_presence[ 9]		= true;
  // fOCLCollimator_presence[ 9]	= true;
  // fOCLLaBr3_Distance[ 9]		= 16.3*cm;
  OCLLaBr3_theta[ 9]			= 100.812191*deg;
  OCLLaBr3_phi[ 9]			= 71.999954*deg;

  //Det. number 10
  // fOCLLaBr3_presence[10]		= true;
  // fOCLCollimator_presence[10]	= true;
  // fOCLLaBr3_Distance[10]		= 16.3*cm;
  OCLLaBr3_theta[10]			= 116.564844*deg;
  OCLLaBr3_phi[10]			= 108.000000*deg;

  //Det. number 11
  // fOCLLaBr3_presence[11]		= true;
  // fOCLCollimator_presence[11]	= true;
  // fOCLLaBr3_Distance[11]		= 16.3*cm;
  OCLLaBr3_theta[11]			= 100.812191*deg;
  OCLLaBr3_phi[11]			= 144.000046*deg;

  //Det. number 12
  // fOCLLaBr3_presence[12]		= true;
  // fOCLCollimator_presence[12]	= true;
  // fOCLLaBr3_Distance[12]		= 16.3*cm;
  OCLLaBr3_theta[12]			= 116.564908*deg;
  OCLLaBr3_phi[12]			= 180.000045*deg;

  //Det. number 13:
  // fOCLLaBr3_presence[13]		= true;
  // fOCLCollimator_presence[13]	= true;
  // fOCLLaBr3_Distance[13]		= 16.3*cm;
  OCLLaBr3_theta[13]			= 100.812208*deg;
  OCLLaBr3_phi[13]			= 215.999956*deg;

  //Det. number 14
  // fOCLLaBr3_presence[14]		= true;
  // fOCLCollimator_presence[14]	= true;
  // fOCLLaBr3_Distance[14]		= 16.3*cm;
  OCLLaBr3_theta[14]			= 116.564893*deg;
  OCLLaBr3_phi[14]			= 252.000011*deg;

  //Det. number 15
  // fOCLLaBr3_presence[15]		= true;
  // fOCLCollimator_presence[15]	= true;
  // fOCLLaBr3_Distance[15]		= 16.3*cm;
  OCLLaBr3_theta[15]			= 100.812175*deg;
  OCLLaBr3_phi[15]			= 288.000000*deg;

  //Det. number 16
  // fOCLLaBr3_presence[16]		= true;
  // fOCLCollimator_presence[16]	= true;
  // fOCLLaBr3_Distance[16]		= 16.3*cm;
  OCLLaBr3_theta[16]			= 79.187575*deg;
  OCLLaBr3_phi[16]			= 324.000046*deg;

  //Det. number 17
  // fOCLLaBr3_presence[17]		= true;
  // fOCLCollimator_presence[17]	= true;
  // fOCLLaBr3_Distance[17]		= 16.3*cm;
  OCLLaBr3_theta[17]			= 63.434885*deg;
  OCLLaBr3_phi[17]			= 0.000045*deg;

  //Det. number 18
  // fOCLLaBr3_presence[18]		= true;
  // fOCLCollimator_presence[18]	= true;
  // fOCLLaBr3_Distance[18]		= 16.3*cm;
  OCLLaBr3_theta[18]			= 79.187559*deg;
  OCLLaBr3_phi[18]			= 35.999956*deg;

  //Det. number 19
  // fOCLLaBr3_presence[19]		= true;
  // fOCLCollimator_presence[19]	= true;
  // fOCLLaBr3_Distance[19]		= 16.3*cm;
  OCLLaBr3_theta[19]			= 63.434900*deg;
  OCLLaBr3_phi[19]			= 72.000011*deg;

  //Det. number 20
  // fOCLLaBr3_presence[20]		= true;
  // fOCLCollimator_presence[20]	= true;
  // fOCLLaBr3_Distance[20]		= 16.3*cm;
  OCLLaBr3_theta[20]			= 79.187591*deg;
  OCLLaBr3_phi[20]			= 108.000000*deg;

  //Det. number 21
  // fOCLLaBr3_presence[21]		= true;
  // fOCLCollimator_presence[21]	= true;
  // fOCLLaBr3_Distance[21]		= 16.3*cm;
  OCLLaBr3_theta[21]			= 63.434900*deg;
  OCLLaBr3_phi[21]			= 143.999989*deg;

  //Det. number 22
  // fOCLLaBr3_presence[22]		= true;
  // fOCLCollimator_presence[22]	= true;
  // fOCLLaBr3_Distance[22]		= 16.3*cm;
  OCLLaBr3_theta[22]			= 79.187559*deg;
  OCLLaBr3_phi[22]			= 180.000044*deg;

  //Det. number 23
  // fOCLLaBr3_presence[23]		= true;
  // fOCLCollimator_presence[23]	= true;
  // fOCLLaBr3_Distance[23]		= 16.3*cm;
  OCLLaBr3_theta[23]			= 63.434885*deg;
  OCLLaBr3_phi[23]			= 215.999955*deg;

  //Det. number 24
  // fOCLLaBr3_presence[24]		= true;
  // fOCLCollimator_presence[24]	= true;
  // fOCLLaBr3_Distance[24]		= 16.3*cm;
  OCLLaBr3_theta[24]			= 79.187575*deg;
  OCLLaBr3_phi[24]			= 251.999954*deg;

  //Det. number 25
  // fOCLLaBr3_presence[25]		= true;
  // fOCLCollimator_presence[25]	= true;
  // fOCLLaBr3_Distance[25]		= 16.3*cm;
  OCLLaBr3_theta[25]			= 63.434949*deg;
  OCLLaBr3_phi[25]			= 288.000000*deg;

  //Det. number 26
  // fOCLLaBr3_presence[26]		= true;
  // fOCLCollimator_presence[26]	= true;
  // fOCLLaBr3_Distance[26]		= 16.3*cm;
  OCLLaBr3_theta[26]			= 37.377328*deg;
  OCLLaBr3_phi[26]			= 324.000134*deg;

  //Det. number 27
  // fOCLLaBr3_presence[27]		= true;
  // fOCLCollimator_presence[27]	= true;
  // fOCLLaBr3_Distance[27]		= 16.3*cm;
  OCLLaBr3_theta[27]			= 37.377294*deg;
  OCLLaBr3_phi[27]			= 36.000153*deg;

  //Det. number 28
  // fOCLLaBr3_presence[28]		= true;
  // fOCLCollimator_presence[28]	= true;
  // fOCLLaBr3_Distance[28]		= 16.3*cm;
  OCLLaBr3_theta[28]			= 37.377321*deg;
  OCLLaBr3_phi[28]			= 108.000000*deg;

  //Det. number 29
  // fOCLLaBr3_presence[29]		= true;
  // fOCLCollimator_presence[29]	= true;
  // fOCLLaBr3_Distance[29]		= 16.3*cm;
  OCLLaBr3_theta[29]			= 37.377294*deg;
  OCLLaBr3_phi[29]			= 179.999847*deg;

  //Det. number 30
  // fOCLLaBr3_presence[30]		= true;
  // fOCLCollimator_presence[30]	= true;
  // fOCLLaBr3_Distance[30]		= 16.3*cm;
  OCLLaBr3_theta[30]			= 37.377328*deg;
  OCLLaBr3_phi[30]			= 251.999866*deg;

  //Beamline
  // fOCLLaBr3_presence[31]		= false;
  // fOCLCollimator_presence[31]	= false;
  // fOCLLaBr3_Distance[31]		= 16.3*cm;
  OCLLaBr3_theta[31]			= 0.000000*deg;
  OCLLaBr3_phi[31]			= 0.000000*deg;



  for(G4int i=0; i<numberOf_OCLLaBr3; i++){

  	G4double disttoLaBr3_face = fOCLLaBr3_Distance[i];
  	// disttoCollHalf =  fOCLLaBr3_Distance[i] - offsettoCollimator;

  	positionOCLLaBr3[i] = SpherToCatG4three(disttoLaBr3_face, OCLLaBr3_theta[i], OCLLaBr3_phi[i]);
  	// if (!rotmOCLLaBr3[i]){
  			// for(G4int i=0; i<numberOf_OCLLaBr3; i++){
    	rotmOCLLaBr3[i] = G4RotationMatrix();
  	// }
  		rotmOCLLaBr3[i].rotateY(OCLLaBr3_theta[i]);
  		rotmOCLLaBr3[i].rotateZ(OCLLaBr3_phi[i]);
    // }
  	// positionCollimator[i] = SpherToCatG4three(disttoCollHalf, OCLLaBr3_theta[i], OCLLaBr3_phi[i]);
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void OCLDetectorConstruction::SetLaBr3_Distance(G4int i, G4double dist){
  fOCLLaBr3_Distance[i] = dist;
  G4RunManager::GetRunManager()->ReinitializeGeometry(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void OCLDetectorConstruction::SetUseLaBr3(G4int i, G4bool use){
  fOCLLaBr3_presence[i] = use;
  G4RunManager::GetRunManager()->ReinitializeGeometry(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void OCLDetectorConstruction::SetUseCSGOldTargetChamber(G4bool use){
  fUseCSGOldTargetChamber = use;
  G4RunManager::GetRunManager()->ReinitializeGeometry(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void OCLDetectorConstruction::SetUseCSGOldTarget(G4bool use){
  fUseCSGOldTarget = use;
  G4RunManager::GetRunManager()->ReinitializeGeometry(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void OCLDetectorConstruction::SetUseCSGRadSource(G4bool use){
  fUseCSGRadSource = use;
  G4RunManager::GetRunManager()->ReinitializeGeometry(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void OCLDetectorConstruction::SetUseCSGSiRi(G4bool use){
  fUseCSGSiRi = use;
  G4RunManager::GetRunManager()->ReinitializeGeometry(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void OCLDetectorConstruction::SetUseCSGNiff(G4bool use){
  fUseCSGNiff = use;
  G4RunManager::GetRunManager()->ReinitializeGeometry(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
