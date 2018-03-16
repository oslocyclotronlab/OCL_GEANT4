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
	G4double density;                 //z=mean number of protons;

	G4int ncomponents, natoms;
	G4double abundance, fractionmass;

	// load NIST material database manager
	G4NistManager * man = G4NistManager::Instance();

	//
	// Define Elements
	//

	// add more elements from NIST database
	C  = man->FindOrBuildElement("C");
	Mn = man->FindOrBuildElement("Mn");
	Si = man->FindOrBuildElement("Si");
	P  = man->FindOrBuildElement("P");
	S  = man->FindOrBuildElement("S");
	N = man->FindOrBuildElement("N");
	Cu = man->FindOrBuildElement("Cu");
	Cr = man->FindOrBuildElement("Cr");
	Ni = man->FindOrBuildElement("Ni");
	Fe = man->FindOrBuildElement("Fe");

	// SteelAISI304 (Stainless Steel)
	// source: approx. from https://www.azom.com/article.aspx?ArticleID=965
	SteelAISI304 =   new G4Material("SteelAISI304", density = 8.000*g/cm3, ncomponents=8);
	SteelAISI304->AddElement(C,  fractionmass=0.01*0.05); // approx avg. of 0.08 and 0.03
	SteelAISI304->AddElement(Mn, fractionmass=0.01*2.0);
	SteelAISI304->AddElement(Si, fractionmass=0.01*0.75);
	SteelAISI304->AddElement(P,  fractionmass=0.01*0.45);
	SteelAISI304->AddElement(S,  fractionmass=0.01*0.030);
	SteelAISI304->AddElement(Cr, fractionmass=0.01*19.); // avg. of 19 and 20;
	SteelAISI304->AddElement(Ni, fractionmass=0.01*9.);  // approx avg. of 8 and 10.5
	SteelAISI304->AddElement(Fe, fractionmass=0.01*68.72); // Rest

	// SteelS235JR (Carbon Steel)
	// source: approx. from http://www.b2bmetal.eu/en/pages/index/index/id/141/
	SteelS235JR =   new G4Material("SteelS235JR", density = 7.85*g/cm3, ncomponents=7);
	SteelS235JR->AddElement(C,  fractionmass=0.01*0.17);
	SteelS235JR->AddElement(Mn, fractionmass=0.01*1.4);
	SteelS235JR->AddElement(P,  fractionmass=0.01*0.040);
	SteelS235JR->AddElement(S,  fractionmass=0.01*0.040);
	SteelS235JR->AddElement(N,  fractionmass=0.01*0.012);
	SteelS235JR->AddElement(Cu, fractionmass=0.01*0.55); 
	SteelS235JR->AddElement(Fe, fractionmass=0.01*97.788); // Rest

	G4_Al_Material  = man->FindOrBuildMaterial("G4_Al");


// Stainless steel: aisi 304 (rings)
// Stainless steel: aisi 316 (stakene)

// Carbon steel: sj 235

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
	offsetFrameBall = G4ThreeVector(0, 0, 0);

	CADMesh* FrameBallMesh = new CADMesh("../OCL/Mesh-Models/Structures/Halv_Kaktus_overvie_Karbonstal_Copy1.stl", "STL", mm, offsetFrameBall, false);

	FrameBallCADSolid = FrameBallMesh->TessellatedMesh();
	FrameBallCADlog = new G4LogicalVolume(FrameBallCADSolid, SteelS235JR, "FrameBallCADLog", 0, 0, 0);

    ///////////////////////////////////////////////////////////////////////////////
	// FrameRing as CAD
	///////////////////////////////////////////////////////////////////////////////

	// mm; Distance of global model "centre" to centre of ball/OSCAR
	offsetFrameRing = G4ThreeVector(0,0,0);

	CADMesh* FrameRingMesh = new CADMesh("../OCL/Mesh-Models/Structures/halv_kaktus_overvie_Aluminium_lowres.stl", "STL", mm, offsetFrameRing, false);

	FrameRingCADSolid = FrameRingMesh->TessellatedMesh();
	// approximate(!) all material here as stainless steel
	FrameRingCADlog = new G4LogicalVolume(FrameRingCADSolid, G4_Al_Material, "SteelAISI304", 0, 0, 0);


    ///////////////////////////////////////////////////////////////////////////////
	// FrameBase as CAD
	///////////////////////////////////////////////////////////////////////////////
	offsetFrameBase = G4ThreeVector(0,0,0);
	CADMesh* FrameBaseMesh = new CADMesh("../OCL/Mesh-Models/Structures/JMC_006_00_Assembly_LowRes.stl", "STL", mm, offsetFrameBase, false);

	FrameBaseCADSolid = FrameBaseMesh->TessellatedMesh();
	FrameBaseCADlog = new G4LogicalVolume(FrameBaseCADSolid, SteelAISI304, "FrameBaseCADLog", 0, 0, 0);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void OCLFrame::Placement(G4int copyNo, G4VPhysicalVolume* physiMother, G4bool checkOverlaps)
{

	///////////////////////////////////////////////////////////////////////////////
	// FrameBall as CAD
	///////////////////////////////////////////////////////////////////////////////
	// rotmFrameBallCAD[0].rotateY(90.*deg); 
	// rotmFrameBallCAD[0].rotateZ(-90.*deg);

	// rotmFrameBallCAD[1].rotateY(-90.*deg); 
	// rotmFrameBallCAD[1].rotateZ(90.*deg);
	rotmFrameBallCAD[0] = G4RotationMatrix(0, 0, 0);
	rotmFrameBallCAD[1] = G4RotationMatrix(0, 0, 0);

	for (int i=0; i<1; i++){
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
	// dx_FrameRing = 190.*mm; // arb number at the moment
	// dz_FrameRing = -20*mm; // arb number at the moment

	// dRot = -80*deg;

	// translatFrameRingCAD[0] = G4ThreeVector(dx_FrameRing,0,dz_FrameRing);
	// translatFrameRingCAD[1] = G4ThreeVector(-dx_FrameRing,0,-dz_FrameRing);

	// rotmFrameRingCAD[0].rotateY(dRot); 
	// rotmFrameRingCAD[1].rotateY(180.*deg + dRot);

	// // rotmFrameRingCAD[1].rotateY(70.*deg); 
	// // rotmFrameRingCAD[1].rotateZ(90.*deg);
	rotmFrameRingCAD[0] = G4RotationMatrix(0, 0, 0);
	rotmFrameRingCAD[1] = G4RotationMatrix(0, 0, 0);
	translatFrameRingCAD[0] = G4ThreeVector(0,0,0);
	translatFrameRingCAD[1] = G4ThreeVector(0,0,0);

	for (int i=0; i<1; i++){
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

 //    /////////////////////////////////////////////////////////////////////////////
	// // FrameBase as CAD
	// ///////////////////////////////////////////////////////////////////////////////
	
	// dx_FrameBase =70*mm; // arb number at the moment
	// dy_FrameBase = 120.*cm; // arb number at the moment

	// dRotFrameBase = 90*deg;

	// translatFrameBaseCAD[0] = G4ThreeVector(dx_FrameBase,-dy_FrameBase,0);
	// translatFrameBaseCAD[1] = G4ThreeVector(-dx_FrameBase,-dy_FrameBase,0);

	// rotmFrameBaseCAD[0].rotateY(dRotFrameBase); 
	// rotmFrameBaseCAD[1].rotateY(180.*deg + dRotFrameBase);

	// // rotmFrameBaseCAD[1].rotateY(70.*deg); 
	// // rotmFrameBaseCAD[1].rotateZ(90.*deg);
	rotmFrameBaseCAD[0] = G4RotationMatrix(0, 0, 0);
	rotmFrameBaseCAD[1] = G4RotationMatrix(0, 0, 0);
	translatFrameBaseCAD[0] = G4ThreeVector(0,0,0);
	translatFrameBaseCAD[1] = G4ThreeVector(0,0,0);



	for (int i=0; i<1; i++){
		transformFrameBaseCAD[i] = G4Transform3D(rotmFrameBaseCAD[i],translatFrameBaseCAD[i]);

		FrameBaseCADphys[i] =
			    new G4PVPlacement(transformFrameBaseCAD[i],
		                          "FrameBaseCADphys",       // its name
		                          FrameBaseCADlog,       // its logical volume
		                          physiMother,         // its mother  volume
		                          false,           // no boolean operations
		                          i,               // copy number
		                          true); // checking overlaps
	}
	//  Visualisation
    FrameBaseCAD_VisAtt = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
    FrameBaseCAD_VisAtt->SetForceSolid(true);
    FrameBaseCADlog->SetVisAttributes(FrameBaseCAD_VisAtt);
 // // //    ////////////////////////////////////////////////////////////////

}