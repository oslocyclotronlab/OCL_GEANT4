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
#include "OCLMaterials.hh"

#include "G4PhysicalConstants.hh"
#include "G4NistManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

OCLFrame::OCLFrame()
{
	//some default Clover detector parameters
	// --> Parameter file


	//----------------------------------------------------
	// Material definitions
	//----------------------------------------------------

  OCLMaterials* fMat = OCLMaterials::GetInstance();

	//
	// Define Elements
	//

	Aluminium  = fMat->GetMaterial("G4_Al");

	//
	// Set all Placement Parameters
	//

	CreatePlacementParameters();

	//
	// Create the solids.....
	//

	CreateSolids();

}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Destructor
OCLFrame::~OCLFrame() {}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void OCLFrame::SetPosition(G4ThreeVector thisPos) {
  translatePos = thisPos*mm;
  G4cout << " ----> A OCLFrame will be placed at distance: " << translatePos/mm << " mm" << G4endl;
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void OCLFrame::SetRotation(G4RotationMatrix thisRot) {
	rotation = thisRot;
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void  OCLFrame::CreateSolids()
{
	//
	// Frame (Ball)
	//

	// Pentagons

	G4double startPhi = 0.*deg;
	G4double deltaPhi = twopi;

	G4double halfThicknessPentagon   = 10.*mm;
	G4double bendAnglePentagon   = 73.5*deg;
	G4double hightPentagon_z01 = 177.*mm;					// larger side (facing the detectors)
	G4double hightPentagon_z00 = hightPentagon_z01 - 2.*2.*halfThicknessPentagon / tan(bendAnglePentagon); //small sider (facing center)

	G4double sideLengthPentagon_z01  = 115.*mm;
    G4double sideLengthPentagon_z00  = hightPentagon_z00
    									* 2. / (sqrt(5.+2.*sqrt(5.)) );

	G4double rOuterPentagon_z01   = sideLengthPentagon_z01  // inradius of the larger side
									/ (2.*tan(pi/5.));
	G4double rOuterPentagon_z00   = sideLengthPentagon_z00  // inradius of the smaller side
									/ (2.*tan(pi/5.));

	G4double zPlanesPentagon[]  = {-halfThicknessPentagon, halfThicknessPentagon};
	G4double rInnerPentagon[]   = {0.*cm,0.*cm};
	G4double rOuterPentagon[]   = {rOuterPentagon_z00, rOuterPentagon_z01};

	G4double rFrameHole = 118.*mm/2.;

	// Hole & Pentagons

  	solidFrameHoles =
		new G4Tubs("FrameHoles",
				0.*cm,	// pRMin
				rFrameHole,	// pRMax
				halfThicknessPentagon, // pDz
				startPhi,
				deltaPhi);

	solidPentagon =
		new G4Polyhedra( "Pentagon",   // name
			             startPhi,   // initial phi starting angle
			             deltaPhi,        // total phi angle
			             5,			   // number sides
			             2,           // number of z planes
			             zPlanesPentagon,    // position of z planes
			       		 rInnerPentagon,    // tangent distance to inner surface
			       		 rOuterPentagon  ); // tangent distance to outer surface


  	subtractFramePentagon =
	    new  G4SubtractionSolid("Pentagon - Hole",
								solidPentagon,  		// 1st object
								solidFrameHoles);	   	// 2nd object


    logicFramePentagon = new G4LogicalVolume(subtractFramePentagon, Aluminium, "FramePentagon");

    // Hexagons

    G4double halfThicknessHexagon   = 10.*mm;
	G4double bendAngleHexagon   = 69.1*deg;
	G4double hightHexagon_z01 = 199.2*mm;					// larger side (facing the detectors)
	G4double hightHexagon_z00 = hightHexagon_z01 - 2.*2.*halfThicknessHexagon / tan(bendAngleHexagon); //small sider (facing center)

	G4double rOuterHexagon_z01   = hightHexagon_z01/2.;  // inradius of the larger side
	G4double rOuterHexagon_z00   = hightHexagon_z00/2.;  // inradius of the smaller side


	G4double zPlanesHexagon[]  = {-halfThicknessHexagon, halfThicknessHexagon};
	G4double rInnerHexagon[]   = {0.*cm,0.*cm};
	G4double rOuterHexagon[]   = {rOuterHexagon_z00, rOuterHexagon_z01};

	solidHexagon =
	new G4Polyhedra( "Hexagon",   // name
		             startPhi,   // initial phi starting angle
		             deltaPhi,        // total phi angle
		             6,			   // number sides
		             2,           // number of z planes
		             zPlanesHexagon,    // position of z planes
		       		 rInnerHexagon,    // tangent distance to inner surface
		       		 rOuterHexagon  ); // tangent distance to outer surface


  	subtractFrameHexagon =
	    new  G4SubtractionSolid("Hexagon - Hole",
								solidHexagon,  		// 1st object
								solidFrameHoles);	   	// 2nd object


    logicFrameHexagon = new G4LogicalVolume(subtractFrameHexagon, Aluminium, "FrameHexagon");


}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void OCLFrame::Placement(G4int copyNo, G4VPhysicalVolume* physiMother, G4bool checkOverlaps)
{

	// Create the transformation vector to mother volume
	// rotation.rotateZ(90.*deg);
	// G4Transform3D translationGlobal = G4Transform3D(rotation,translatePos);


 	//
	// Frame (Ball) Geometry
	//

	for(G4int i=0; i<numberOf_Pentagons; i++){
		G4Transform3D transformPentagon = G4Transform3D(rotmPentagon[i],positionPentagon[i]);
		physiFramePentagon =
			new G4PVPlacement(transformPentagon, 	// Transformation (Rot&Transl)
							  "FramePentagon",		// its name
							   logicFramePentagon, 	// its logical volume
							   physiMother, 		// its physical mother volume
							   false, 				// unknown "pMany"; def: false
							   i , 					// copy number
							   checkOverlaps);		// checkOverlaps
	}

	for(G4int i=0; i<numberOf_Hexagons; i++){
		G4Transform3D transformHexagon = G4Transform3D(rotmHexagon[i],positionHexagon[i]);
		physiFrameHexagon =
			new G4PVPlacement(transformHexagon, 	// Transformation (Rot&Transl)
							  "FrameHexagon",		// its name
							   logicFrameHexagon, 	// its logical volume
							   physiMother, 		// its physical mother volume
							   false, 				// unknown "pMany"; def: false
							   i, 					// copy number
							   checkOverlaps);		// checkOverlaps
	}





	// //------------------------------------------------------
	// // visualization attributes
	// //------------------------------------------------------

	// Frame (Ball)
	VisAtt1= new G4VisAttributes(G4Colour(0.5,0.5,0.5)); //White
	logicFramePentagon->SetVisAttributes(VisAtt1);
	logicFrameHexagon->SetVisAttributes(VisAtt1);

}

G4ThreeVector SpherToCatG4three(G4double r,G4double theta,G4double phi);

void OCLFrame::CreatePlacementParameters()
{

G4ThreeVector w;

G4double distToPentagonHalf = 240.738*mm + 10*mm; // dist. to inner side of Pentagon + halfwidth
G4double distToHexagonHalf  = 247.66 *mm + 10*mm; // dist. to inner side of Hexagon  + halfwidth

frameHexagon_theta[ 0]	= 142.622528*deg;
frameHexagon_phi[ 0]	= 144.000134*deg;

frameHexagon_theta[ 1]	= 142.622561*deg;
frameHexagon_phi[ 1]	= 216.000153*deg;

frameHexagon_theta[ 2]	= 142.622535*deg;
frameHexagon_phi[ 2]	= 288.000000*deg;

frameHexagon_theta[ 3]	= 142.622561*deg;
frameHexagon_phi[ 3]	= 359.999847*deg;

frameHexagon_theta[ 4]	= 142.622528*deg;
frameHexagon_phi[ 4]	= 71.999866*deg;

frameHexagon_theta[ 5]	= 100.812191*deg;
frameHexagon_phi[ 5]	= 144.000046*deg;

frameHexagon_theta[ 6]	= 79.187559*deg;
frameHexagon_phi[ 6]	= 180.000044*deg;

frameHexagon_theta[ 7]	= 100.812208*deg;
frameHexagon_phi[ 7]	= 215.999956*deg;

frameHexagon_theta[ 8]	= 79.187575*deg;
frameHexagon_phi[ 8]	= 251.999954*deg;

frameHexagon_theta[ 9]	= 100.812175*deg;
frameHexagon_phi[ 9]	= 288.000000*deg;

frameHexagon_theta[10]	= 79.187575*deg;
frameHexagon_phi[10]	= 324.000046*deg;

frameHexagon_theta[11]	= 100.812208*deg;
frameHexagon_phi[11]	= 0.000044*deg;

frameHexagon_theta[12]	= 79.187559*deg;
frameHexagon_phi[12]	= 35.999956*deg;

frameHexagon_theta[13]	= 100.812191*deg;
frameHexagon_phi[13]	= 71.999954*deg;

frameHexagon_theta[14]	= 79.187591*deg;
frameHexagon_phi[14]	= 108.000000*deg;

frameHexagon_theta[15]	= 37.377321*deg;
frameHexagon_phi[15]	= 108.000000*deg;

frameHexagon_theta[16]	= 37.377294*deg;
frameHexagon_phi[16]	= 36.000153*deg;

frameHexagon_theta[17]	= 37.377328*deg;
frameHexagon_phi[17]	= 324.000134*deg;

frameHexagon_theta[18]	= 37.377328*deg;
frameHexagon_phi[18]	= 251.999866*deg;

frameHexagon_theta[19]	= 37.377294*deg;
frameHexagon_phi[19]	= 179.999847*deg;

framePentagon_theta[ 0]	= 180.000000*deg;
framePentagon_phi[ 0]	= 180.000000*deg;

framePentagon_theta[ 1]	= 116.564844*deg;
framePentagon_phi[ 1]	= 108.000000*deg;

framePentagon_theta[ 2]	= 116.564908*deg;
framePentagon_phi[ 2]	= 180.000045*deg;

framePentagon_theta[ 3]	= 116.564893*deg;
framePentagon_phi[ 3]	= 252.000011*deg;

framePentagon_theta[ 4]	= 116.564893*deg;
framePentagon_phi[ 4]	= 323.999989*deg;

framePentagon_theta[ 5]	= 116.564908*deg;
framePentagon_phi[ 5]	= 35.999955*deg;

framePentagon_theta[ 6]	= 0.000000*deg;
framePentagon_phi[ 6]	= 0.000000*deg;

framePentagon_theta[ 7]	= 63.434900*deg;
framePentagon_phi[ 7]	= 72.000011*deg;

framePentagon_theta[ 8]	= 63.434885*deg;
framePentagon_phi[ 8]	= 0.000045*deg;

framePentagon_theta[ 9]	= 63.434949*deg;
framePentagon_phi[ 9]	= 288.000000*deg;

framePentagon_theta[10]	= 63.434885*deg;
framePentagon_phi[10]	= 215.999955*deg;

framePentagon_theta[11]	= 63.434900*deg;
framePentagon_phi[11]	= 143.999989*deg;

double rot1 = 36.*deg;

rotmPentagon[0].rotateZ(-180*deg);
rotmPentagon[7].rotateZ(rot1);
rotmPentagon[8].rotateZ(rot1);
rotmPentagon[9].rotateZ(rot1);
rotmPentagon[10].rotateZ(rot1);
rotmPentagon[11].rotateZ(rot1);


for(G4int i=0; i<numberOf_Pentagons; i++){

	positionPentagon[i] = SpherToCatG4three(distToPentagonHalf, framePentagon_theta[i], framePentagon_phi[i]);
	rotmPentagon[i].rotateY(framePentagon_theta[i]);
	rotmPentagon[i].rotateZ(framePentagon_phi[i]);
}

rotmHexagon[0].rotateZ(90*deg);
rotmHexagon[4].rotateZ(90*deg);

rotmHexagon[1].rotateZ(30*deg);
rotmHexagon[2].rotateZ(30*deg);
rotmHexagon[3].rotateZ(30*deg);

rotmHexagon[5].rotateZ(30*deg);
rotmHexagon[6].rotateZ(30*deg);
rotmHexagon[7].rotateZ(30*deg);
rotmHexagon[8].rotateZ(30*deg);
rotmHexagon[9].rotateZ(30*deg);
rotmHexagon[10].rotateZ(30*deg);
rotmHexagon[11].rotateZ(30*deg);
rotmHexagon[12].rotateZ(30*deg);
rotmHexagon[13].rotateZ(30*deg);
rotmHexagon[14].rotateZ(30*deg);
rotmHexagon[15].rotateZ(30*deg);
rotmHexagon[16].rotateZ(30*deg);
rotmHexagon[17].rotateZ(30*deg);
rotmHexagon[18].rotateZ(30*deg);
rotmHexagon[19].rotateZ(30*deg);

for(G4int i=0; i<numberOf_Hexagons; i++){
	positionHexagon[i] = SpherToCatG4three(distToHexagonHalf, frameHexagon_theta[i], frameHexagon_phi[i]);
	rotmHexagon[i].rotateY(frameHexagon_theta[i]);
	rotmHexagon[i].rotateZ(frameHexagon_phi[i]);
}



}
