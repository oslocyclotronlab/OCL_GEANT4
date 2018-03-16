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
#include "OCLParameters.hh"
#include "OCLLaBr3.hh"
#include "OCLCollimator.hh"
#include "OCLFrame.hh"
#include "SiRi.hh"
#include "NIFF.hh"
#include "TChamberAuspuff.hh"
#include "OCLTarget.hh"

#include "G4SystemOfUnits.hh"
#include "G4Transform3D.hh"
#include "G4PhysicalConstants.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
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
     solidWorld(0), WorldLog(0), physiWorld(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

OCLDetectorConstruction::~OCLDetectorConstruction()
{}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* OCLDetectorConstruction::Construct()
{

   // vacuum (non-STP)
	
    G4Material* vacuum = 
    	new G4Material("Vacuum",       				//name as String
						1,		                    //atomic number (use 1 for Hydrogen)
                    	1.008*g/mole, 	            //molar mass (use 1.008*g/mole for Hydoren)
						1.e-25*g/cm3,  	            //density
						kStateGas,		            //kStateGas - the material is gas (see G4State)
                    	2.73*kelvin,	            //Temperature
						1.e-25*g/cm3);	            //pressure


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
	               world_sizeXYZ/2, 
	               world_sizeXYZ/2, 
	               world_sizeXYZ/2); 

	WorldLog =  
		new G4LogicalVolume(solidWorld,        	//solid defining the World
	                    	vacuum,           	//material of the World
	                    	"World");         	//name

	physiWorld = 
		new G4PVPlacement(0,                    //specifies rotation: 0 = no rotation
	                  	  G4ThreeVector(),     	//at (0,0,0)
	                  	  WorldLog,            	//logical volume
						  "World",              //name
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
	if( OCLLaBr3_presence[i])
		{
		copynumber=i;
		labr3[i] = new OCLLaBr3();
		labr3[i]->SetRotation(rotmOCLLaBr3[i]);
		labr3[i]->SetPosition(positionOCLLaBr3[i]);
		labr3[i]->Placement(copynumber,  physiWorld, pSurfChk);
		}
	}
	 

	///////////

	//
	// Collimator
	//

	for(G4int i=0; i<numberOf_OCLLaBr3; i++){
	if( OCLCollimator_presence[i])
		{
		copynumber=i;
		collimator[i] = new OCLCollimator();
		collimator[i]->SetRotation(rotmOCLLaBr3[i]); // same rotation as Detector
		collimator[i]->SetPosition(positionCollimator[i]);
		// collimator[i]->Placement(copynumber,  physiWorld, pSurfChk);
		}
	}


	//
	//	Frame
	//
	G4RotationMatrix rotmFrame = G4RotationMatrix(); 	// Check implementation: Need to be empty?
	G4ThreeVector 	 positionFrame = G4ThreeVector(); 	//  

	OCLFrame* frame;
	frame = new OCLFrame();
	// frame->SetRotation(rotmFrame);
	// frame->SetPosition(positionFrame);
	frame->Placement(0,  physiWorld, pSurfChk);


	//
	// Target Chamber
	//

	TChamberAuspuff* chamber;
	chamber = new TChamberAuspuff();
	chamber->Placement(0,  physiWorld, pSurfChk);

	//
	// Target
	//

	OCLTarget* target;
	target = new OCLTarget();
	target->Placement(0,  physiWorld, pSurfChk);


	//
	//	SiRi
	//
	// G4RotationMatrix rotmSiRi = G4RotationMatrix(); 	// Check implementation: Need to be empty?
	// G4ThreeVector 	 positionSiRi = G4ThreeVector();	//  

	SiRi* siri;
	siri = new SiRi();
	siri->SetAngle(137*deg); // backward
	// siri->SetAngle(43*deg); // forward
	
	// siri->SetRotation(rotmSiRi);
	// siri->SetPosition(positionSiRi);
	siri->Placement(0,  physiWorld, pSurfChk);

	//
	//	NIFF
	//
	// G4RotationMatrix rotmNIFF = G4RotationMatrix(); 	// Check implementation: Need to be empty?
	// G4ThreeVector 	 positionNIFF = G4ThreeVector();	//  

	NIFF* niff;
	niff = new NIFF();
	// // niff->SetRotation(rotmNIFF);
	// // niff->SetPosition(positionNIFF);
	// niff->Placement(0,  physiWorld, pSurfChk);



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

G4ThreeVector SpherToCatG4three(G4double r,G4double theta,G4double phi){
	return r * G4ThreeVector( sin(theta) * cos(phi),
					   sin(theta) * sin(phi), 
					   cos(theta));
}

void OCLDetectorConstruction::SetPlacementParameters()
{

G4double distColltoDet = 10*mm; // Distance between Collimator and Detector (surface to Surf)

//Beamline            
OCLLaBr3_presence[0]		= false;             
OCLCollimator_presence[0]	= false;             
OCLLaBr3_Distance[0]		= 20*cm;             
OCLLaBr3_theta[00]			= 180.000000*deg;             
OCLLaBr3_phi[0]			= 180.000000*deg;

//Det. number 1 at OCL:            
OCLLaBr3_presence[ 1]		= true;             
OCLCollimator_presence[ 1]	= true;             
OCLLaBr3_Distance[ 1]		= 20*cm;             
OCLLaBr3_theta[ 1]			= 142.622561*deg;             
OCLLaBr3_phi[ 1]			= 359.999847*deg;

//Det. number 2 at OCL             
OCLLaBr3_presence[ 2]		= true;             
OCLCollimator_presence[ 2]	= true;             
OCLLaBr3_Distance[ 2]		= 20*cm;             
OCLLaBr3_theta[ 2]			= 142.622528*deg;             
OCLLaBr3_phi[ 2]			= 71.999866*deg;

//Det. number 3 at OCL:
OCLLaBr3_presence[ 3]		= true;             
OCLCollimator_presence[ 3]	= true;             
OCLLaBr3_Distance[ 3]		= 20*cm;             
OCLLaBr3_theta[ 3]			= 142.622528*deg;             
OCLLaBr3_phi[ 3]			= 144.000134*deg;

//Det. number 4 at OCL             
OCLLaBr3_presence[ 4]		= true;             
OCLCollimator_presence[ 4]	= true;             
OCLLaBr3_Distance[ 4]		= 20*cm;             
OCLLaBr3_theta[ 4]			= 142.622561*deg;             
OCLLaBr3_phi[ 4]			= 216.000153*deg;
             
//Det. number 5 at OCL  
OCLLaBr3_presence[ 5]		= true;             
OCLCollimator_presence[ 5]	= true;             
OCLLaBr3_Distance[ 5]		= 20*cm;             
OCLLaBr3_theta[ 5]			= 142.622535*deg;             
OCLLaBr3_phi[ 5]			= 288.000000*deg;

//Det. number 6           
OCLLaBr3_presence[6]		= true;             
OCLCollimator_presence[6]	= true;             
OCLLaBr3_Distance[6]		= 20*cm;             
OCLLaBr3_theta[6]			= 116.564893*deg;             
OCLLaBr3_phi[6]				= 323.999989*deg;

//Det. number 7
OCLLaBr3_presence[ 7]		= true;             
OCLCollimator_presence[ 7]	= true;             
OCLLaBr3_Distance[ 7]		= 20*cm;             
OCLLaBr3_theta[ 7]			= 100.812208*deg;             
OCLLaBr3_phi[ 7]			= 0.000044*deg;

//Det. number 8           
OCLLaBr3_presence[ 8]		= true;             
OCLCollimator_presence[ 8]	= true;             
OCLLaBr3_Distance[ 8]		= 20*cm;             
OCLLaBr3_theta[ 8]			= 116.564908*deg;             
OCLLaBr3_phi[ 8]			= 35.999955*deg;

//Det. number 9 at OCL:             
OCLLaBr3_presence[ 9]		= true;             
OCLCollimator_presence[ 9]	= true;             
OCLLaBr3_Distance[ 9]		= 20*cm;             
OCLLaBr3_theta[ 9]			= 100.812191*deg;             
OCLLaBr3_phi[ 9]			= 71.999954*deg;

//Det. number 10             
OCLLaBr3_presence[10]		= true;             
OCLCollimator_presence[10]	= true;             
OCLLaBr3_Distance[10]		= 20*cm;             
OCLLaBr3_theta[10]			= 116.564844*deg;             
OCLLaBr3_phi[10]			= 108.000000*deg;

//Det. number 11        
OCLLaBr3_presence[11]		= true;             
OCLCollimator_presence[11]	= true;             
OCLLaBr3_Distance[11]		= 20*cm;             
OCLLaBr3_theta[11]			= 100.812191*deg;             
OCLLaBr3_phi[11]			= 144.000046*deg;

//Det. number 12             
OCLLaBr3_presence[12]		= true;             
OCLCollimator_presence[12]	= true;             
OCLLaBr3_Distance[12]		= 20*cm;             
OCLLaBr3_theta[12]			= 116.564908*deg;             
OCLLaBr3_phi[12]			= 180.000045*deg;

//Det. number 13:             
OCLLaBr3_presence[13]		= true;             
OCLCollimator_presence[13]	= true;             
OCLLaBr3_Distance[13]		= 20*cm;             
OCLLaBr3_theta[13]			= 100.812208*deg;             
OCLLaBr3_phi[13]			= 215.999956*deg;

//Det. number 14             
OCLLaBr3_presence[14]		= true;             
OCLCollimator_presence[14]	= true;             
OCLLaBr3_Distance[14]		= 20*cm;             
OCLLaBr3_theta[14]			= 116.564893*deg;             
OCLLaBr3_phi[14]			= 252.000011*deg;

//Det. number 15             
OCLLaBr3_presence[15]		= true;             
OCLCollimator_presence[15]	= true;             
OCLLaBr3_Distance[15]		= 20*cm;             
OCLLaBr3_theta[15]			= 100.812175*deg;             
OCLLaBr3_phi[15]			= 288.000000*deg;

//Det. number 16           
OCLLaBr3_presence[16]		= true;             
OCLCollimator_presence[16]	= true;             
OCLLaBr3_Distance[16]		= 20*cm;             
OCLLaBr3_theta[16]			= 79.187575*deg;             
OCLLaBr3_phi[16]			= 324.000046*deg;

//Det. number 17             
OCLLaBr3_presence[17]		= true;             
OCLCollimator_presence[17]	= true;             
OCLLaBr3_Distance[17]		= 20*cm;             
OCLLaBr3_theta[17]			= 63.434885*deg;             
OCLLaBr3_phi[17]			= 0.000045*deg;

//Det. number 18      
OCLLaBr3_presence[18]		= true;             
OCLCollimator_presence[18]	= true;             
OCLLaBr3_Distance[18]		= 20*cm;             
OCLLaBr3_theta[18]			= 79.187559*deg;             
OCLLaBr3_phi[18]			= 35.999956*deg;

//Det. number 19             
OCLLaBr3_presence[19]		= true;             
OCLCollimator_presence[19]	= true;             
OCLLaBr3_Distance[19]		= 20*cm;             
OCLLaBr3_theta[19]			= 63.434900*deg;             
OCLLaBr3_phi[19]			= 72.000011*deg;

//Det. number 20             
OCLLaBr3_presence[20]		= true;             
OCLCollimator_presence[20]	= true;             
OCLLaBr3_Distance[20]		= 20*cm;             
OCLLaBr3_theta[20]			= 79.187591*deg;             
OCLLaBr3_phi[20]			= 108.000000*deg;

//Det. number 21             
OCLLaBr3_presence[21]		= true;             
OCLCollimator_presence[21]	= true;             
OCLLaBr3_Distance[21]		= 20*cm;             
OCLLaBr3_theta[21]			= 63.434900*deg;             
OCLLaBr3_phi[21]			= 143.999989*deg;

//Det. number 22             
OCLLaBr3_presence[22]		= true;             
OCLCollimator_presence[22]	= true;             
OCLLaBr3_Distance[22]		= 20*cm;             
OCLLaBr3_theta[22]			= 79.187559*deg;             
OCLLaBr3_phi[22]			= 180.000044*deg;

//Det. number 23             
OCLLaBr3_presence[23]		= true;             
OCLCollimator_presence[23]	= true;             
OCLLaBr3_Distance[23]		= 20*cm;             
OCLLaBr3_theta[23]			= 63.434885*deg;             
OCLLaBr3_phi[23]			= 215.999955*deg;

//Det. number 24             
OCLLaBr3_presence[24]		= true;             
OCLCollimator_presence[24]	= true;             
OCLLaBr3_Distance[24]		= 20*cm;             
OCLLaBr3_theta[24]			= 79.187575*deg;             
OCLLaBr3_phi[24]			= 251.999954*deg;

//Det. number 25
OCLLaBr3_presence[25]		= true;             
OCLCollimator_presence[25]	= true;             
OCLLaBr3_Distance[25]		= 20*cm;             
OCLLaBr3_theta[25]			= 63.434949*deg;             
OCLLaBr3_phi[25]			= 288.000000*deg;

//Det. number 26             
OCLLaBr3_presence[26]		= true;             
OCLCollimator_presence[26]	= true;             
OCLLaBr3_Distance[26]		= 20*cm;             
OCLLaBr3_theta[26]			= 37.377328*deg;             
OCLLaBr3_phi[26]			= 324.000134*deg;

//Det. number 27             
OCLLaBr3_presence[27]		= true;             
OCLCollimator_presence[27]	= true;             
OCLLaBr3_Distance[27]		= 20*cm;             
OCLLaBr3_theta[27]			= 37.377294*deg;             
OCLLaBr3_phi[27]			= 36.000153*deg;

//Det. number 28           
OCLLaBr3_presence[28]		= true;             
OCLCollimator_presence[28]	= true;             
OCLLaBr3_Distance[28]		= 20*cm;             
OCLLaBr3_theta[28]			= 37.377321*deg;             
OCLLaBr3_phi[28]			= 108.000000*deg;
 
//Det. number 29             
OCLLaBr3_presence[29]		= true;             
OCLCollimator_presence[29]	= true;             
OCLLaBr3_Distance[29]		= 20*cm;             
OCLLaBr3_theta[29]			= 37.377294*deg;             
OCLLaBr3_phi[29]			= 179.999847*deg;

//Det. number 30             
OCLLaBr3_presence[30]		= true;             
OCLCollimator_presence[30]	= true;             
OCLLaBr3_Distance[30]		= 20*cm;             
OCLLaBr3_theta[30]			= 37.377328*deg;             
OCLLaBr3_phi[30]			= 251.999866*deg;

//Beamline             
OCLLaBr3_presence[31]		= false;             
OCLCollimator_presence[31]	= false;             
OCLLaBr3_Distance[31]		= 20*cm;             
OCLLaBr3_theta[31]			= 0.000000*deg;             
OCLLaBr3_phi[31]			= 0.000000*deg;


G4double disttoLaBr3Half;
G4double disttoCollHalf;
G4double offsettoCollimator = distColltoDet + collimatorHalfLength;


for(G4int i=0; i<numberOf_OCLLaBr3; i++){
	

	disttoLaBr3Half = OCLLaBr3_Distance[i] 
					  + housingInnerHalfLength;
	disttoCollHalf =  OCLLaBr3_Distance[i] - offsettoCollimator;

	positionOCLLaBr3[i] = SpherToCatG4three(disttoLaBr3Half, OCLLaBr3_theta[i], OCLLaBr3_phi[i]); 
	rotmOCLLaBr3[i].rotateZ(30.*deg); 
	rotmOCLLaBr3[i].rotateY(OCLLaBr3_theta[i]); 
	rotmOCLLaBr3[i].rotateZ(OCLLaBr3_phi[i]);

	positionCollimator[i] = SpherToCatG4three(disttoCollHalf, OCLLaBr3_theta[i], OCLLaBr3_phi[i]); 
}

}