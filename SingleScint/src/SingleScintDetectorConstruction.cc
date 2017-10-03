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


#include "SingleScintDetectorConstruction.hh"
#include "SingleScintParameters.hh"
#include "OCLLaBr3.hh"
#include "OCLCollimator.hh"


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

SingleScintDetectorConstruction::SingleScintDetectorConstruction()
:
     solidWorld(0), WorldLog(0), WorldPhys(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SingleScintDetectorConstruction::~SingleScintDetectorConstruction()
{}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* SingleScintDetectorConstruction::Construct()
{

   // vacuum (non-STP)

    G4Material* vacuum = new G4Material("Vacuum",       //name as String
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


 ////////////////////////
 // Positinging
 ////////////////////////


G4double offsettoCollimator = 10*cm;          // Distance from source to Collimator (beginning);
G4double phi = 45*deg, 
G4double theta = 45*deg;

G4RotationMatrix rotm1 = G4RotationMatrix();
rotm1.rotateY(theta); 
rotm1.rotateZ(phi);
G4cout << "\n --> phi = " << phi/deg << " deg;  direct rotation matrix : ";
// rotm1.print(G4cout);   

//
// possitioning
//

G4int copynumber;
copynumber = 0;
bool pSurfChk = false;

//
// LaBr3
//

G4double disttoLaBr3Half = offsettoCollimator + 2*collimatorHalfLength +  detectorHalfinclPMT;

G4ThreeVector w = G4ThreeVector( std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));    
G4ThreeVector positionLaBr3 = disttoLaBr3Half*w;

OCLLaBr3* labr3;
labr3 = new OCLLaBr3();
labr3->SetRotation(rotm1);
labr3->SetPosition(positionLaBr3);
labr3->Placement(copynumber,  physiWorld, pSurfChk);

///////////

//
// Collimator
//
G4double disttoCol = offsettoCollimator + collimatorHalfLength;

// possitioning     
w = G4ThreeVector( std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));    
G4ThreeVector positionCol = disttoCol*w;

OCLCollimator* collimator;
collimator = new OCLCollimator();
collimator->SetRotation(rotm1);
collimator->SetPosition(positionCol);
collimator->Placement(copynumber,  physiWorld, pSurfChk);




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
