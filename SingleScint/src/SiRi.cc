/*
	SiRi module
*/

#include "SiRi.hh"

SiRi::SiRi()
: G4VUserDetectorConstruction(),
	E_Log(0)
{  
  
  //----------------------------------------------------
  // Material definitions
  //----------------------------------------------------

  G4NistManager* nist = G4NistManager::Instance();
  silicon = nist->FindOrBuildMaterial("G4_Si");
  
  
   // vacuum (non-STP)
 
  vacuum = new G4Material("Vacuum",       //name as String
              1,          //atomic number (use 1 for Hydrogen)
                            1.008*g/mole,   //molar mass (use 1.008*g/mole for Hydoren) 
              1.e-25*g/cm3,   //density
              kStateGas,    //kStateGas - the material is gas (see G4State)
                            2.73*kelvin,  //Temperature 
              1.e-25*g/cm3);  //pressure
   

  //
  // Create the solids.....
  //

  CreateSolids();
}


SiRi::~SiRi()
{ }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Dummy class -- doesn't compile otherwise
G4VPhysicalVolume* SiRi::Construct(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SiRi::SetPosition(G4ThreeVector thisPos) {
  translatePos = thisPos*mm;
  G4cout << " ----> SiRi will be placed at distance: " << translatePos/mm << " mm" << G4endl;
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SiRi::SetRotation(G4RotationMatrix thisRot) { 
  rotation = thisRot; 
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void  SiRi::CreateSolids() {

  //
  // Detector Geometry
  //


	 // G4ThreeVector posE = G4ThreeVector(0, 0, 0);					// position of the E detector
	 // G4ThreeVector pos_deltaE = G4ThreeVector(0, -5.*mm, -5.*mm);			// position of the delta E detector
 	 
 	 G4RotationMatrix rMat;									// rotation matrix around X axis
   G4double phiX = -135 *deg;						// 90° (we want the detector up) + 47° (the angle)
 		// G4double phiZ = 45. *deg;
   rMat.rotateX(phiX);
 		// rMat.rotateZ(phiZ);
  

  	// dimensions 
  	a = 20*mm;
  	h = 13.6*mm;				//height
  	alfa = 3*pi/8;
  	b = a + 2*(h/tan(alfa));
  	d_deltaE = 0.13*mm;			// Delta E detector thickness: 130 micron
  	d_E = 1.55*mm;				// E detector thickness: 1550 micron
  	
  	dist = 5*cm;				// distance of the center of the thick detector from the source
  	
 	  // G4Trd* Delta_E_shape = new G4Trd("Delta_E", .5*a, .5*b, .5*d_deltaE, .5*d_deltaE, .5*h);
 	
 	  // G4Trd **DE_shape = new G4Trd*[8];
 	  // G4LogicalVolume **logVolume = new G4LogicalVolume*[8];
 	  // G4String nameShape[8] = {"shape1", "shape2", "shape3", "shape4", "shape5", "shape6", "shape7", "shape8"};
 	  // G4String nameLog[8] = {"pad40_LV", "pad42_LV", "pad44_LV", "pad46_LV", "pad48_LV", "pad50_LV", "pad52_LV", "pad54_LV"};
 	  // G4String placeLog[8] = {"Delta_E_1", "Delta_E_2", "Delta_E_3", "Delta_E_4", "Delta_E_5", "Delta_E_6", "Delta_E_7", "Delta_E_8"};
 	  
    //
    //  thick detector
    //
  	G4Trd* E_shape = new G4Trd("E_det", .5*a, .5*b, .5*d_E, .5*d_E, .5*h);
  	E_Log = new G4LogicalVolume(E_shape, silicon, "ThickLV");
 	  // new G4PVPlacement(G4Transform3D(rMat,posE), E_Log, "E_det", WorldLog, false, 0);
  	  
    //
    //	  thin detector 
    //
		/* We will define each of the 8 pads...by hand (yay!!)
		   If you know how to write this in a smarter way, go for it!
		   I did not get my string arrays to work.   */
	  
	  h_thin = 1.5*mm;
	  // G4double y0 = -(2.5 + h/2 - h_thin/2)*mm;
	  // G4double z0 = -1.5*mm;
	  G4double c[8];
	  // G4double y[8], z[8];
	  // G4ThreeVector posvector[8];
 
	  for(int i=0; i<8; i++){
	  	c[i] = a + (2+2*i)*(h_thin/tan(alfa));
	    	// y[i] = y0 + i*(1.2)*mm;					// step by 1.2 mm looks the best
	    	// z[i] = z0 -i*1.2*mm;					// same as above
	    	// posvector[i] = G4ThreeVector(0.,y[i],z[i]);
	  }// end of for loop

    //
    // now define shapes and logical volumes	  
	  //
    // 1
	  G4Trd* DE_shape1 = new G4Trd("Delta_E1", .5*a, .5*c[0], .5*d_deltaE, .5*d_deltaE, .5*h_thin);
	  DE_Log1 = new G4LogicalVolume(DE_shape1, silicon, "ThinLV1");
	  // 2
	  G4Trd* DE_shape2 = new G4Trd("Delta_E2", .5*c[0], .5*c[1], .5*d_deltaE, .5*d_deltaE, .5*h_thin);
	  DE_Log2 = new G4LogicalVolume(DE_shape2, silicon, "ThinLV2");
	  // 3
	  G4Trd* DE_shape3 = new G4Trd("Delta_E3", .5*c[1], .5*c[2], .5*d_deltaE, .5*d_deltaE, .5*h_thin);
	  DE_Log3 = new G4LogicalVolume(DE_shape3, silicon, "ThinLV3");
	  // 4
	  G4Trd* DE_shape4 = new G4Trd("Delta_E4", .5*c[2], .5*c[3], .5*d_deltaE, .5*d_deltaE, .5*h_thin);
	  DE_Log4 = new G4LogicalVolume(DE_shape4, silicon, "ThinLV4");
	  // 5
	  G4Trd* DE_shape5 = new G4Trd("Delta_E5", .5*c[3], .5*c[4], .5*d_deltaE, .5*d_deltaE, .5*h_thin);
	  DE_Log5 = new G4LogicalVolume(DE_shape5, silicon, "ThinLV5");
	  // 6
	  G4Trd* DE_shape6 = new G4Trd("Delta_E6", .5*c[4], .5*c[5], .5*d_deltaE, .5*d_deltaE, .5*h_thin);
	  DE_Log6 = new G4LogicalVolume(DE_shape6, silicon, "ThinLV6");
	  // 7
	  G4Trd* DE_shape7 = new G4Trd("Delta_E7", .5*c[5], .5*c[6], .5*d_deltaE, .5*d_deltaE, .5*h_thin);
	  DE_Log7 = new G4LogicalVolume(DE_shape7, silicon, "ThinLV7");
	  // 8
	  G4Trd* DE_shape8 = new G4Trd("Delta_E8", .5*c[6], .5*c[7], .5*d_deltaE, .5*d_deltaE, .5*h_thin);
	  DE_Log8 = new G4LogicalVolume(DE_shape8, silicon, "ThinLV8");
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void SiRi::Placement(G4int copyNo, G4VPhysicalVolume* physiMother, G4bool checkOverlaps)
{

  // place pads in a ring
  /* explain coordinate system */ 
  G4double nb_pads = 8;				// number of pads in each layer
  G4double dPhi = twopi/nb_pads;		// change in angle for each padS
  G4double theta = 47*deg;			// angle of the detector wrt to beam
  G4double r_delta = 5*mm;			// distance between thin and thick detector
  G4double R_E = dist*tan(theta);				// radius of the "circle with back pads", measured to the centre of the pads
  G4double y0 = -(2.5 + h/2 - h_thin/2)*mm;		// y-coordinate of the centre of the bottom strip, minus sign due to the definition of the coordinate system
  G4double z0 = -1.5*mm;					// distance between bottom strip and the back pad, minus sign due to the definition of the coordinate system
  	  
 	// G4ThreeVector uz_delta[8];
  G4ThreeVector pos_deltaE[8];
  G4RotationMatrix rotMat[8];
  	  

  for (G4int ipad = 0; ipad < nb_pads ; ipad++) {
    G4double phi = ipad*dPhi;
    G4double right_angle = 90*deg;
    G4double rot_angle = phi - right_angle;			// DO NOT CHANGE UNDER ANY CIRCUMSTANCES!
    G4RotationMatrix rotm  = G4RotationMatrix();		// global rotation matrix
    rotm.rotateX(-137*deg); 
    // rotm.rotateY(90*deg);
    rotm.rotateZ(rot_angle);

	for(int j=0; j<8; j++){						// strip loop
	
	  	G4double R_pad = R_E + y0 + j*(1.2)*mm;			// radius of each strip
	    	// y[i] = y0 + j*(1.2)*mm;					 
	    G4double z = z0 -j*1.2*mm;					 
	 
    	G4ThreeVector uz_delta = G4ThreeVector((1-r_delta/R_E)*std::cos(phi),(1-r_delta/R_E)*std::sin(phi),z/R_pad);    	// unitary vector
    	pos_deltaE[j] = R_pad*uz_delta;														// scale the unit vector to get the position of the strip
    	G4RotationMatrix rotmPad  = G4RotationMatrix();
    	G4double strip_angle = (40 + j*2)*deg;										// angle of each strip wrt beam axis
    	G4double X_angle = -(right_angle + strip_angle);
    	rotmPad.rotateX(X_angle);
    	rotmPad.rotateZ(rot_angle);
    	rotMat[j] = rotmPad;
        
	}	// j-loop
	
    G4ThreeVector uz = G4ThreeVector(std::cos(phi),  std::sin(phi),0.); 
    G4ThreeVector posE = R_E*uz;
    
    // thick detector                                
    E_phys = new G4PVPlacement(G4Transform3D(rotm,posE), "E_det", E_Log, physiMother, ipad, checkOverlaps);
    
    // thin strips (front detector)
    DE_phys1 = new G4PVPlacement(G4Transform3D(rotMat[0],pos_deltaE[0]), "Delta_E1", DE_Log1, physiMother, ipad, checkOverlaps);
    DE_phys2 = new G4PVPlacement(G4Transform3D(rotMat[1],pos_deltaE[1]), "Delta_E2", DE_Log2, physiMother, ipad, checkOverlaps);
    DE_phys3 = new G4PVPlacement(G4Transform3D(rotMat[2],pos_deltaE[2]), "Delta_E3", DE_Log3, physiMother, ipad, checkOverlaps);
    DE_phys4 = new G4PVPlacement(G4Transform3D(rotMat[3],pos_deltaE[3]), "Delta_E4", DE_Log4, physiMother, ipad, checkOverlaps);
    DE_phys5 = new G4PVPlacement(G4Transform3D(rotMat[4],pos_deltaE[4]), "Delta_E5", DE_Log5, physiMother, ipad, checkOverlaps);
    DE_phys6 = new G4PVPlacement(G4Transform3D(rotMat[5],pos_deltaE[5]), "Delta_E6", DE_Log6, physiMother, ipad, checkOverlaps);
    DE_phys7 = new G4PVPlacement(G4Transform3D(rotMat[6],pos_deltaE[6]), "Delta_E7", DE_Log7, physiMother, ipad, checkOverlaps);
    DE_phys8 = new G4PVPlacement(G4Transform3D(rotMat[7],pos_deltaE[7]), "Delta_E8", DE_Log8, physiMother, ipad, checkOverlaps);

  } //ipad-loop

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void SiRi::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

// declare trackers as a MultiFunctionalDetector scorer
  // 1 
  G4MultiFunctionalDetector* Si_deltaE1 = new G4MultiFunctionalDetector("Si_thin1");		// Delta E detector
  G4VPrimitiveScorer* primitiv1 = new G4PSEnergyDeposit("edep1_1");
  Si_deltaE1->RegisterPrimitive(primitiv1);
  SetSensitiveDetector("ThinLV1",Si_deltaE1);
  // 2
  G4MultiFunctionalDetector* Si_deltaE2 = new G4MultiFunctionalDetector("Si_thin2");		// Delta E detector
  G4VPrimitiveScorer* primitiv2 = new G4PSEnergyDeposit("edep1_2");
  Si_deltaE2->RegisterPrimitive(primitiv2);
  SetSensitiveDetector("ThinLV2",Si_deltaE2);
  // 3
  G4MultiFunctionalDetector* Si_deltaE3 = new G4MultiFunctionalDetector("Si_thin3");		// Delta E detector
  G4VPrimitiveScorer* primitiv3 = new G4PSEnergyDeposit("edep1_3");
  Si_deltaE3->RegisterPrimitive(primitiv3);
  SetSensitiveDetector("ThinLV3",Si_deltaE3);
  // 4
  G4MultiFunctionalDetector* Si_deltaE4 = new G4MultiFunctionalDetector("Si_thin4");		// Delta E detector
  G4VPrimitiveScorer* primitiv4 = new G4PSEnergyDeposit("edep1_4");
  Si_deltaE4->RegisterPrimitive(primitiv4);
  SetSensitiveDetector("ThinLV4",Si_deltaE4);
  // 5
  G4MultiFunctionalDetector* Si_deltaE5 = new G4MultiFunctionalDetector("Si_thin5");		// Delta E detector
  G4VPrimitiveScorer* primitiv5 = new G4PSEnergyDeposit("edep1_5");
  Si_deltaE5->RegisterPrimitive(primitiv5);
  SetSensitiveDetector("ThinLV5",Si_deltaE5);
  // 6
  G4MultiFunctionalDetector* Si_deltaE6 = new G4MultiFunctionalDetector("Si_thin6");		// Delta E detector
  G4VPrimitiveScorer* primitiv6 = new G4PSEnergyDeposit("edep1_6");
  Si_deltaE6->RegisterPrimitive(primitiv6);
  SetSensitiveDetector("ThinLV6",Si_deltaE6);
  // 7
  G4MultiFunctionalDetector* Si_deltaE7 = new G4MultiFunctionalDetector("Si_thin7");		// Delta E detector
  G4VPrimitiveScorer* primitiv7 = new G4PSEnergyDeposit("edep1_7");
  Si_deltaE7->RegisterPrimitive(primitiv7);
  SetSensitiveDetector("ThinLV7",Si_deltaE7);
  // 8 
  G4MultiFunctionalDetector* Si_deltaE8 = new G4MultiFunctionalDetector("Si_thin8");		// Delta E detector
  G4VPrimitiveScorer* primitiv8 = new G4PSEnergyDeposit("edep1_8");
  Si_deltaE8->RegisterPrimitive(primitiv8);
  SetSensitiveDetector("ThinLV8",Si_deltaE8);
  
  //.......................
  G4MultiFunctionalDetector* Si_E = new G4MultiFunctionalDetector("Si_thick");		// E detector
  G4VPrimitiveScorer* primitiv9 = new G4PSEnergyDeposit("edep2");
  Si_E->RegisterPrimitive(primitiv9);
  SetSensitiveDetector("ThickLV",Si_E);

}

  
                      

