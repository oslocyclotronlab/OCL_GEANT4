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

#include "OCLCollimator.hh"
#include "OCLMaterials.hh"
#include "OCLParameters.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

OCLCollimator::OCLCollimator()
{
	//some default Clover detector parameters
	// --> Parameter file


	//----------------------------------------------------
	// Material definitions
	//----------------------------------------------------
	// Get materials
	OCLMaterials* fMat = OCLMaterials::GetInstance();

	//Lead
	lead = fMat->GetMaterial("G4_Pb");

	////////////////////////////////////////////////////////

	//create the solids.....
	CreateSolids();

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//Destructor
OCLCollimator::~OCLCollimator() { }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void OCLCollimator::SetPosition(G4ThreeVector thisPos) {
  translatePos = thisPos*mm;
  G4cout << " ----> A OCLCollimator will be placed at " << translatePos/mm << " mm" << G4endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void OCLCollimator::SetRotation(G4RotationMatrix thisRot) { rotation = thisRot; }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void  OCLCollimator::CreateSolids()
{

	//
	// Collimator
	//

    // Parameters are now in the header file "Parameters.hh"


	// Collimator geometry

	solidCollimator =
		new G4Cons("ShieldingConnical",
					colRmin1,   // inner radius = 0 because used as mother volume
					colRmax1,
					colRmin2,   // inner radius = 0 because used as mother volume
					colRmax2,
					collimatorHalfLength,
					startPhi,
					deltaPhi);

	logicCollimator = new G4LogicalVolume(solidCollimator, lead, "Collimator");

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void OCLCollimator::Placement(G4int copyNo, G4VPhysicalVolume* physiMother, G4bool checkOverlaps)
{

	//
	// Collimator
	//

	G4Transform3D transCollimator = G4Transform3D(rotation,translatePos);
	physiCollimator =
		new G4PVPlacement(transCollimator,		// Transformation (Rot&Transl)
							"Collimator",		// its name
							logicCollimator,	// its logical volume
							physiMother,		// its physical mother volume
							false,				// unknown "pMany"; def: false
							copyNo,				// copy number
							checkOverlaps);		// checkOverlaps

}
