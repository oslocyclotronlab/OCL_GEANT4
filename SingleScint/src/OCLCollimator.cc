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
#include "SingleScintParameters.hh"

// #include "G4VPhysicalVolume.hh"
// #include "G4LogicalVolume.hh"
// #include "G4Box.hh"
// #include "G4Tubs.hh"
// #include "G4Cons.hh"
// #include "G4UnionSolid.hh"
// #include "G4Material.hh"
#include "G4NistManager.hh"
// #include "G4PVPlacement.hh"
// #include "G4VisAttributes.hh"
// #include "G4SystemOfUnits.hh"
// #include "G4Transform3D.hh"
// #include "G4PhysicalConstants.hh"

// #include "G4Colour.hh"

// //#include "G4MultiFunctionalDetector.hh"
// //#include "G4VPrimitiveScorer.hh"
// //#include "G4PSEnergyDeposit.hh"
// //#include "G4TransportationManager.hh"
// //#include "G4SDManager.hh"

// #include "G4LogicalBorderSurface.hh"
// #include "G4OpticalSurface.hh"


OCLCollimator::OCLCollimator()
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

  //Lead
  lead =    man->FindOrBuildMaterial("G4_Pb");

  ////////////////////////////////////////////////////////

  //create the solids.....
  CreateSolids();

}

//Destructor
OCLCollimator::~OCLCollimator() { }

void OCLCollimator::SetPosition(G4ThreeVector thisPos) {
  translatePos = thisPos*mm;
  G4cout << " ----> A OCLCollimator will be placed at " << translatePos/mm << " mm" << G4endl;
}

void OCLCollimator::SetRotation(G4RotationMatrix thisRot) { rotation = thisRot; }

//---------------------------------------------------------------------
// Create the solids defining Phase-II Clovers
//---------------------------------------------------------------------
void  OCLCollimator::CreateSolids()
{

	//
	// Collimator
	//

    // Parameters are now in the header file "Parameters.hh"


	// Collimator geometry

	solidCollimator = new G4Cons("ShieldingConnical",
											colRmin1,   // inner radius = 0 because used as mother volume
											colRmax1,
											colRmin2,   // inner radius = 0 because used as mother volume
											colRmax2,
											collimatorHalfLength,
											startPhi,
											deltaPhi);

	logicCollimator = new G4LogicalVolume(solidCollimator, lead, "Collimator");

}


//------------------------------------------------------------------
void OCLCollimator::Placement(G4int copyNo, G4VPhysicalVolume* physiMother, G4bool checkOverlaps)
{

	//
	// Collimator
	//

	positionCollimator = G4ThreeVector(0.*cm,0.*cm,0.*cm);
	positionCollimator += translatePos;

	G4Transform3D transCollimator = G4Transform3D(rotation,translatePos);
	physiCollimator = new G4PVPlacement(transCollimator,
														   "Collimator",
														   logicCollimator,
														   physiMother,
														   false,copyNo,checkOverlaps);

}
