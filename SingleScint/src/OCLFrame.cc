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

	Aluminium  = man->FindOrBuildMaterial("G4_Al");

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
    									/ (tan(pi/5.) * sqrt(5.+2.*sqrt(.5)) );

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


    logicFramePentagon = new G4LogicalVolume(solidPentagon, Aluminium, "FramePentagon");

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


    logicFrameHexagon = new G4LogicalVolume(solidHexagon, Aluminium, "FrameHexagon");


}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void OCLFrame::Placement(G4int copyNo, G4VPhysicalVolume* physiMother, G4bool checkOverlaps)
{

	// Create the transformation vector to mother volume
	// rotation.rotateZ(90.*deg);
	// G4Transform3D translationGlobal = G4Transform3D(rotation,translatePos);

	G4int copyNoSub = 0; // copy number for the sub elements (could also be copyNo)

 	//
	// Detector Geometry
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

	// // White color for Detector
	// VisAtt1= new G4VisAttributes(G4Colour(1.0,1.0,1.0)); //White
	// logicCrystal->SetVisAttributes(VisAtt1);

	// // Red color for Shielding
	// VisAtt2 = new G4VisAttributes(G4Colour(1.0,0.0,0.0)); //red
	// logicShielding->SetVisAttributes(VisAtt2);

}


G4ThreeVector SpherToCatG4three(G4double r,G4double theta,G4double phi){
	return r * G4ThreeVector( sin(theta) * cos(phi),
					   sin(theta) * sin(phi), 
					   cos(theta));
}


void OCLFrame::CreatePlacementParameters()
{
frameHexagon_theta[0]  = 69.094720*deg;      
frameHexagon_phi[ 0]  = 148.282473*deg;     
frameHexagon_psi[ 0]  = 0/twopi;
     
frameHexagon_theta[ 1]  = 54.735580*deg;      
frameHexagon_phi[ 1]  = 193.282664*deg;     
frameHexagon_psi[ 1]  = 0/twopi;
     
frameHexagon_theta[ 2]  = 90.000000*deg;      
frameHexagon_phi[ 2]  = 217.377465*deg;     
frameHexagon_psi[ 2]  = 0/twopi;
     
frameHexagon_theta[ 3]  = 125.264420*deg;     
frameHexagon_phi[ 3]  = 193.282664*deg;     
frameHexagon_psi[ 3]  = 0/twopi;
     
frameHexagon_theta[ 4]  = 110.905280*deg;     
frameHexagon_phi[ 4]  = 148.282473*deg;     
frameHexagon_psi[ 4]  = 0/twopi;
     
frameHexagon_theta[ 5]  = 54.735549*deg;      
frameHexagon_phi[ 5]  = 103.282380*deg;     
frameHexagon_psi[ 5]  = 0/twopi;
     
frameHexagon_theta[ 6]  = 20.905182*deg;      
frameHexagon_phi[ 6]  = 58.282163*deg;      
frameHexagon_psi[ 6]  = 0/twopi;
     
frameHexagon_theta[ 7]  = 20.905065*deg;      
frameHexagon_phi[ 7]  = 238.282731*deg;     
frameHexagon_psi[ 7]  = 0/twopi;
     
frameHexagon_theta[ 8]  = 54.735580*deg;      
frameHexagon_phi[ 8]  = 283.282664*deg;     
frameHexagon_psi[ 8]  = 0/twopi;
     
frameHexagon_theta[ 9]  = 90.000000*deg;      
frameHexagon_phi[ 9]  = 259.187825*deg;     
frameHexagon_psi[ 9]  = 0/twopi;
     
frameHexagon_theta[10]  = 125.264420*deg;     
frameHexagon_phi[10]  = 283.282664*deg;     
frameHexagon_psi[10]  = 0/twopi;
     
frameHexagon_theta[11]  = 159.094935*deg;     
frameHexagon_phi[11]  = 238.282731*deg;     
frameHexagon_psi[11]  = 0/twopi;
     
frameHexagon_theta[12]  = 159.094818*deg;     
frameHexagon_phi[12]  = 58.282163*deg;      
frameHexagon_psi[12]  = 0/twopi;
     
frameHexagon_theta[13]  = 125.264451*deg;     
frameHexagon_phi[13]  = 103.282380*deg;     
frameHexagon_psi[13]  = 0/twopi;
     
frameHexagon_theta[14]  = 90.000000*deg;      
frameHexagon_phi[14]  = 79.187591*deg;      
frameHexagon_psi[14]  = 0/twopi;
     
frameHexagon_theta[15]  = 90.000000*deg;      
frameHexagon_phi[15]  = 37.377321*deg;      
frameHexagon_psi[15]  = 0/twopi;
     
frameHexagon_theta[16]  = 125.264286*deg;     
frameHexagon_phi[16]  = 13.282597*deg;      
frameHexagon_psi[16]  = 0/twopi;
     
frameHexagon_theta[17]  = 110.905208*deg;     
frameHexagon_phi[17]  = 328.282607*deg;     
frameHexagon_psi[17]  = 0/twopi;
     
frameHexagon_theta[18]  = 69.094792*deg;      
frameHexagon_phi[18]  = 328.282607*deg;     
frameHexagon_psi[18]  = 0/twopi;
     
frameHexagon_theta[19]  = 54.735714*deg;      
frameHexagon_phi[19]  = 13.282597*deg;      
frameHexagon_psi[19]  = 0/twopi;
     
framePentagon_theta[ 0] = 90.000000*deg;      
framePentagon_phi[ 0] = 180.000000*deg;     
framePentagon_psi[ 0] = 0/twopi;
     
framePentagon_theta[ 1] = 90.000000*deg;      
framePentagon_phi[ 1] = 116.564844*deg;     
framePentagon_psi[ 1] = 0/twopi;
     
framePentagon_theta[ 2] = 31.717335*deg;      
framePentagon_phi[ 2] = 148.282427*deg;     
framePentagon_psi[ 2] = 0/twopi;
     
framePentagon_theta[ 3] = 58.282486*deg;      
framePentagon_phi[ 3] = 238.282706*deg;     
framePentagon_psi[ 3] = 0/twopi;
     
framePentagon_theta[ 4] = 121.717514*deg;     
framePentagon_phi[ 4] = 238.282706*deg;     
framePentagon_psi[ 4] = 0/twopi;
     
framePentagon_theta[ 5] = 148.282665*deg;     
framePentagon_phi[ 5] = 148.282427*deg;     
framePentagon_psi[ 5] = 0/twopi;
     
framePentagon_theta[ 6] = 90.000000*deg;      
framePentagon_phi[ 6] = 0.000000*deg;     
framePentagon_psi[ 6] = 0/twopi;
     
framePentagon_theta[ 7] = 121.717450*deg;     
framePentagon_phi[ 7] = 58.282475*deg;      
framePentagon_psi[ 7] = 0/twopi;
     
framePentagon_theta[ 8] = 148.282498*deg;     
framePentagon_phi[ 8] = 328.282658*deg;     
framePentagon_psi[ 8] = 0/twopi;
     
framePentagon_theta[ 9] = 90.000000*deg;      
framePentagon_phi[ 9] = 296.565051*deg;     
framePentagon_psi[ 9] = 0/twopi;
     
framePentagon_theta[10] = 31.717502*deg;      
framePentagon_phi[10] = 328.282658*deg;     
framePentagon_psi[10] = 0/twopi;
     
framePentagon_theta[11] = 58.282550*deg;      
framePentagon_phi[11] = 58.282475*deg;      
framePentagon_psi[11] = 0/twopi;

G4ThreeVector w;
G4double distToPentagonHalf = 240.9 *mm+10*mm; // TODO!!!
G4double distToHexagonHalf  = 247.66*mm    +10*mm; // TODO!!!

double arbrot1 = (90.-72.);

rotmPentagon[1].rotateZ(180*deg);
rotmPentagon[9].rotateZ(180*deg);

rotmPentagon[2].rotateZ(-arbrot1*deg);
rotmPentagon[3].rotateZ(arbrot1*deg);
rotmPentagon[4].rotateZ(-arbrot1*deg);
rotmPentagon[5].rotateZ(arbrot1*deg);

rotmPentagon[7].rotateZ(-arbrot1*deg);
rotmPentagon[8].rotateZ(arbrot1*deg);
rotmPentagon[10].rotateZ(-arbrot1*deg);
rotmPentagon[11].rotateZ(arbrot1*deg);



for(G4int i=0; i<numberOf_Pentagons; i++){
	
	positionPentagon[i] = SpherToCatG4three(distToPentagonHalf, framePentagon_theta[i], framePentagon_phi[i]); 

	rotmPentagon[i].rotateZ(-90*deg);
	rotmPentagon[i].rotateY(framePentagon_theta[i]); 
	rotmPentagon[i].rotateZ(framePentagon_phi[i]);
}

double arbrot = (90.-72.)/2.;

rotmHexagon[0].rotateZ(90*deg);
rotmHexagon[4].rotateZ(90*deg);

rotmHexagon[17].rotateZ(30*deg);
rotmHexagon[18].rotateZ(-30*deg);

rotmHexagon[6].rotateZ(30*deg);
rotmHexagon[7].rotateZ(30*deg);

rotmHexagon[11].rotateZ(30*deg);
rotmHexagon[12].rotateZ(30*deg);


rotmHexagon[1].rotateZ(arbrot*deg);
rotmHexagon[3].rotateZ(-arbrot*deg);

rotmHexagon[16].rotateZ(-arbrot*deg);
rotmHexagon[19].rotateZ(arbrot*deg);

rotmHexagon[5].rotateZ(-arbrot*deg);
rotmHexagon[8].rotateZ(-arbrot*deg);

rotmHexagon[10].rotateZ(arbrot*deg);
rotmHexagon[13].rotateZ(arbrot*deg);

for(G4int i=0; i<numberOf_Hexagons; i++){
	positionHexagon[i] = SpherToCatG4three(distToHexagonHalf, frameHexagon_theta[i], frameHexagon_phi[i]); 


	rotmHexagon[i].rotateY(frameHexagon_theta[i]); 
	rotmHexagon[i].rotateZ(frameHexagon_phi[i]);
}



}