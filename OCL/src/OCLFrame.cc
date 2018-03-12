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


#include "OCLFrame.hh"

#include "G4PhysicalConstants.hh"
#include "G4NistManager.hh"

#include "CADMesh.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

OCLFrame::OCLFrame()
{
	//some default Clover detector parameters
	// --> Parameter file


	//----------------------------------------------------
	// Material definitions
	//----------------------------------------------------

	// G4double a, z;                    //a=mass of a mole;
	// G4double density;                 //z=mean number of protons;

	// G4int ncomponents, natoms;
	// G4double abundance, fractionmass;

	// load NIST material database manager
	G4NistManager * man = G4NistManager::Instance();

	//
	// Define Elements
	//

	G4_Al_Material  = man->FindOrBuildMaterial("G4_Al");

	//
	// Create the solids.....
	//

	CreateSolids();

}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Destructor
OCLFrame::~OCLFrame() {}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// void OCLFrame::SetPosition(G4ThreeVector thisPos) {
//   translatePos = thisPos*mm;
//   G4cout << " ----> A OCLFrame will be placed at distance: " << translatePos/mm << " mm" << G4endl;
// }




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// void OCLFrame::SetRotation(G4RotationMatrix thisRot) { 
// 	rotation = thisRot; 
// }





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void  OCLFrame::CreateSolids()
{
	///////////////////////////////////////////////////////////////////////////////
	// FrameBall as CAD
	///////////////////////////////////////////////////////////////////////////////

	// mm; Distance of global model "centre" to centre of ball/OSCAR
	dx_from_center = -267.6555*mm; // previously -267.655*mm, bu this lead to a small overlap.
	offsetFrameBall = G4ThreeVector(dx_from_center, 0, 0);

	CADMesh* FrameBallMesh = new CADMesh("../OCL/Mesh-Models/Structures/JMC_006_Assembly_BallHalf_Shell.stl", "STL", cm, offsetFrameBall, false);

	FrameBallCADSolid = FrameBallMesh->TessellatedMesh();
	FrameBallCADlog = new G4LogicalVolume(FrameBallCADSolid, G4_Al_Material, "FrameBallCADLog", 0, 0, 0);

    ///////////////////////////////////////////////////////////////////////////////
	// FrameRing as CAD
	///////////////////////////////////////////////////////////////////////////////

	// mm; Distance of global model "centre" to centre of ball/OSCAR
	offsetFrameRing = G4ThreeVector(0,0,0);

	CADMesh* FrameRingMesh = new CADMesh("../OCL/Mesh-Models/Structures/JMC_006_Frame_Ring.stl", "STL", cm, offsetFrameRing, false);

	FrameRingCADSolid = FrameRingMesh->TessellatedMesh();
	FrameRingCADlog = new G4LogicalVolume(FrameRingCADSolid, G4_Al_Material, "FrameRingCADLog", 0, 0, 0);

    ///////////////////////////////////////////////////////////////////////////////
	// FrameTopBase as CAD
	///////////////////////////////////////////////////////////////////////////////
	offsetFrameTopBase = G4ThreeVector(0,0,0);
	CADMesh* FrameTopBaseMesh = new CADMesh("../OCL/Mesh-Models/Structures/JMC_006_Frame_TopBase.stl", "STL", cm, offsetFrameTopBase, false);

	FrameTopBaseCADSolid = FrameTopBaseMesh->TessellatedMesh();
	FrameTopBaseCADlog = new G4LogicalVolume(FrameTopBaseCADSolid, G4_Al_Material, "FrameTopBaseCADLog", 0, 0, 0);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void OCLFrame::Placement(G4int copyNo, G4VPhysicalVolume* physiMother, G4bool checkOverlaps)
{

	///////////////////////////////////////////////////////////////////////////////
	// FrameBall as CAD
	///////////////////////////////////////////////////////////////////////////////
	rotmFrameBallCAD[0].rotateY(90.*deg); 
	rotmFrameBallCAD[0].rotateZ(-90.*deg);

	rotmFrameBallCAD[1].rotateY(-90.*deg); 
	rotmFrameBallCAD[1].rotateZ(90.*deg);
	

	for (int i=0; i<2; i++){
		translatFrameBallCAD[i] = G4ThreeVector();
		transformFrameBallCAD[i] = G4Transform3D(rotmFrameBallCAD[i],translatFrameBallCAD[i]);

		FrameBallCADphys[i] =
			    new G4PVPlacement(transformFrameBallCAD[i],
		                          "FrameBallCADphys",       // its name
		                          FrameBallCADlog,       // its logical volume
		                          physiMother,         // its mother  volume
		                          false,           // no boolean operations
		                          i,               // copy number
		                          true); // checking overlaps
	}
	//  Visualisation
    FrameBallCAD_VisAtt = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
    FrameBallCAD_VisAtt->SetForceSolid(true);
    FrameBallCADlog->SetVisAttributes(FrameBallCAD_VisAtt);
        

    ///////////////////////////////////////////////////////////////////////////////
	// FrameRing as CAD
	///////////////////////////////////////////////////////////////////////////////
	dx_FrameRing = 190.*mm; // arb number at the moment
	dz_FrameRing = -20*mm; // arb number at the moment

	dRot = -80*deg;

	translatFrameRingCAD[0] = G4ThreeVector(dx_FrameRing,0,dz_FrameRing);
	translatFrameRingCAD[1] = G4ThreeVector(-dx_FrameRing,0,-dz_FrameRing);

	rotmFrameRingCAD[0].rotateY(dRot); 
	rotmFrameRingCAD[1].rotateY(180.*deg + dRot);

	// rotmFrameRingCAD[1].rotateY(70.*deg); 
	// rotmFrameRingCAD[1].rotateZ(90.*deg);

	for (int i=0; i<2; i++){
		transformFrameRingCAD[i] = G4Transform3D(rotmFrameRingCAD[i],translatFrameRingCAD[i]);

		FrameRingCADphys[i] =
			    new G4PVPlacement(transformFrameRingCAD[i],
		                          "FrameRingCADphys",       // its name
		                          FrameRingCADlog,       // its logical volume
		                          physiMother,         // its mother  volume
		                          false,           // no boolean operations
		                          i,               // copy number
		                          true); // checking overlaps
	}
	//  Visualisation
    FrameRingCAD_VisAtt = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
    FrameRingCAD_VisAtt->SetForceSolid(true);
    FrameRingCADlog->SetVisAttributes(FrameRingCAD_VisAtt);
    ////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////
	// FrameTopBase as CAD
	///////////////////////////////////////////////////////////////////////////////
	
	dx_FrameTopBase =70*mm; // arb number at the moment
	dy_FrameTopBase = 120.*cm; // arb number at the moment

	dRotFrameTopBase = 90*deg;

	translatFrameTopBaseCAD[0] = G4ThreeVector(dx_FrameTopBase,-dy_FrameTopBase,0);
	translatFrameTopBaseCAD[1] = G4ThreeVector(-dx_FrameTopBase,-dy_FrameTopBase,0);

	rotmFrameTopBaseCAD[0].rotateY(dRotFrameTopBase); 
	rotmFrameTopBaseCAD[1].rotateY(180.*deg + dRotFrameTopBase);

	// rotmFrameTopBaseCAD[1].rotateY(70.*deg); 
	// rotmFrameTopBaseCAD[1].rotateZ(90.*deg);
	

	for (int i=0; i<2; i++){
		transformFrameTopBaseCAD[i] = G4Transform3D(rotmFrameTopBaseCAD[i],translatFrameTopBaseCAD[i]);

		FrameTopBaseCADphys[i] =
			    new G4PVPlacement(transformFrameTopBaseCAD[i],
		                          "FrameTopBaseCADphys",       // its name
		                          FrameTopBaseCADlog,       // its logical volume
		                          physiMother,         // its mother  volume
		                          false,           // no boolean operations
		                          i,               // copy number
		                          true); // checking overlaps
	}
	//  Visualisation
    FrameTopBaseCAD_VisAtt = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
    FrameTopBaseCAD_VisAtt->SetForceSolid(true);
    FrameTopBaseCADlog->SetVisAttributes(FrameTopBaseCAD_VisAtt);
    ////////////////////////////////////////////////////////////////

}