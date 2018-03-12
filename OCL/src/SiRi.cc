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
  copper = nist->FindOrBuildMaterial("G4_Cu");
  aluminum = nist->FindOrBuildMaterial("G4_Al");
  
  
   // vacuum (non-STP)
 
  vacuum = new G4Material("Vacuum",       //name as String
              1,          //atomic number (use 1 for Hydrogen)
                            1.008*g/mole,   //molar mass (use 1.008*g/mole for Hydoren) 
              1.e-25*g/cm3,   //density
              kStateGas,    //kStateGas - the material is gas (see G4State)
                            2.73*kelvin,  //Temperature 
              1.e-25*g/cm3);  //pressure
   

  //
  // Geometry
  //

  dist = 5*cm;        // distance of the center of the thick detector from the source

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

// void SiRi::SetPosition(G4ThreeVector thisPos) {
//   translatePos = thisPos*mm;
//   G4cout << " ----> SiRi will be placed at distance: " << translatePos/mm << " mm" << G4endl;
// }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SiRi::SetAngle(G4double thisAngle) { 
  angle = thisAngle;
  if ( angle == 137*deg ) {    // backwards angles
    theta = 133*deg;     // angle of the detector wrt to beam
    translatePos = G4ThreeVector(0,0,dist*cos(theta));
    anotherAngle = 47*deg;
    pmone = -1;
    }
  else if ( angle == 43*deg ) {  // forward angles
    theta = 47*deg;      // angle of the detector wrt to beam
    translatePos = G4ThreeVector(0,0,dist*cos(theta));
    anotherAngle = -137*deg;
    pmone = +1;
    }
  else {
    G4cout << "SiRi implemented only for two agles: 43 & 135*deg;\n"
           <<  "You have choosen another angle" << G4endl;
    exit (EXIT_FAILURE);
    }
  G4cout << " ----> SiRi will be placed at this angle: " << thisAngle << G4endl;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void  SiRi::CreateSolids() {

  //
  // Detector Geometry
  //


  // dimensions 
  a = 20*mm;
  h = 13.6*mm;				//height
  alfa = 3*pi/8;
  b = a + 2*(h/tan(alfa));
  d_deltaE = 0.13*mm;			// Delta E detector thickness: 130 micron
  d_E = 1.55*mm;				// E detector thickness: 1550 micron

  r_delta = 0.01*mm;      // distance between thin and thick detector

  // Support structure that holds SiRi
  rInSiriHolder =  50.* mm; // approx.
  thicknessSiriHolder =  1.5* mm; // approx.
  halfLengthSiriHolder = 5.*cm/2.; // approx length
  
  //
  //  thick detector
  //
  G4Trd* E_shape = new G4Trd("E_det", .5*a, .5*b, .5*d_E, .5*d_E, .5*h);
  E_Log = new G4LogicalVolume(E_shape, silicon, "ThickLV");
    
  //
  //	  thin detector 
  //
  /* We will define each of the 8 pads...by hand (yay!!)
     If you know how to write this in a smarter way, go for it!
     I did not get my string arrays to work.   */

  h_thin = 1.5*mm;
  // G4double y0 = -(2.5 + h/2 - h_thin/2)*mm;
  // G4double z0 = -1.5*mm;
  G4double c[9];
  // G4double y[8], z[8];
  // G4ThreeVector posvector[8];

  for(int i=0; i<9; i++){
  	c[i] = a + (2*i)*(h_thin/tan(alfa));
    	// y[i] = y0 + i*(1.2)*mm;					// step by 1.2 mm looks the best
    	// z[i] = z0 -i*1.2*mm;					// same as above
    	// posvector[i] = G4ThreeVector(0.,y[i],z[i]);
  }// end of for loop

  //
  // now define shapes and logical volumes	  
  //
  for(int j=0; j<nFront; j++){ 
  DE_shape[j] = new G4Trd("Delta_E1", .5*c[j], .5*c[j+1], .5*d_deltaE, .5*d_deltaE, .5*h_thin);
  DE_Log[j] = new G4LogicalVolume(DE_shape[j], silicon, "ThinLV");
  }

  // logical shape to gather all dE's
  G4Trd* dE_shape = new G4Trd("dE_det", .5*a, .5*b, .5*(d_deltaE+r_delta), .5*(d_deltaE+r_delta), .5*h);
  dE_Log = new G4LogicalVolume(dE_shape, vacuum, "dE_det");

    // logical shape to gather a pad
  G4Trd* pad_shape = new G4Trd("dE_det", .5*a, .5*b, .5*(d_deltaE+d_E+r_delta), .5*(d_deltaE+d_E+r_delta), .5*h);
  pad_Log = new G4LogicalVolume(pad_shape, vacuum, "SiRi_pad");

  // Support structure that holds SiRi
  G4Tubs* solidSiriHolder = new G4Tubs("SiriHolder",      //pName,
                                     rInSiriHolder,            //pRMin,
                                     rInSiriHolder+thicknessSiriHolder,        //pRMax,
                                     halfLengthSiriHolder, //pDz,
                                     0*deg,           //pSPhi,
                                     360*deg     );    //pDPhi)
  logSiriHolder = new G4LogicalVolume(solidSiriHolder, aluminum, "SiriHolder");

  // Cabels
  rCableCu = 0.5*cm/2.;
  rCableAl = 1.*cm/2.;

  halfLengthCable = (45.5/2.*cm)/2;
  zTransCable = 0.8*cm; // arbitrary
  // halfLengthCable -= 2*cm; // arb red

  G4Tubs* solidCableCu = new G4Tubs("SiRiCableCu",      //pName,
                                     0*cm,            //pRMin,
                                     rCableCu,        //pRMax,
                                     halfLengthCable, //pDz,
                                     0*deg,           //pSPhi,
                                     360*deg     );    //pDPhi)
  logCableCu = new G4LogicalVolume(solidCableCu, copper, "SiRiCableCu");

  G4Tubs* solidCableAl = new G4Tubs("SiRiCableAl",      //pName,
                                     0*cm,            //pRMin,
                                     rCableAl,        //pRMax,
                                     halfLengthCable, //pDz,
                                     0*deg,           //pSPhi,
                                     360*deg     );    //pDPhi)
  logCable = new G4LogicalVolume(solidCableAl, aluminum, "SiRiCableAl");
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void SiRi::Placement(G4int copyNo, G4VPhysicalVolume* physiMother, G4bool checkOverlaps)
{

  // place pads in a ring
  /* explain coordinate system */ 
  G4double nb_pads = 8;				// number of pads in each layer
  G4double dPhi = twopi/nb_pads;		// change in angle for each padS
  G4double R_E = (dist + d_deltaE + r_delta)*sin(theta);				// radius of the "circle with back pads", measured to the centre of the pads
  G4double y0 = -(2.5 + h/2. - h_thin/2.)*mm;		// y-coordinate of the centre of the bottom strip, minus sign due to the definition of the coordinate system
  G4double z0 = -1.5*mm;					// distance between bottom strip and the back pad, minus sign due to the definition of the coordinate system

  G4cout << "R_E" <<  R_E << G4endl;

  G4ThreeVector pos_deltaE[8];
  G4RotationMatrix rotMat[8];
  	  

  // set the strips in the logical frame
  for(int ipad=0; ipad<8; ipad++){ 
      G4ThreeVector dEtranslate = G4ThreeVector(0,
                                                0,
                                                (-3.5+ipad) * h/8);
      
      new G4PVPlacement( 0,
                         dEtranslate,        //solidTransform,
                         DE_Log[ipad],       //pCurrentLogical,
                         "dE_strip_cp",      //pName,
                         dE_Log,             //pMotherLogical,
                         false,              //pMany (for future use),
                         ipad,               //pCopyNo,
                         checkOverlaps);     //pSurfChk=false )
   }

   dE_phys = new G4PVPlacement(0, 
                               G4ThreeVector(0,+d_E/2.,0), 
                               dE_Log, 
                               "SiRi_front", 
                               pad_Log, 
                               false, 
                               0, 
                               checkOverlaps);

   E_phys = new G4PVPlacement(0,  
                              G4ThreeVector(0,-(d_deltaE+r_delta)/2.,0),
                              E_Log, 
                              "SiRi_back", 
                              pad_Log,
                              false,
                              0, 
                              checkOverlaps);

   physSiriHolder = new G4PVPlacement(0,  
                              translatePos+G4ThreeVector(0,0,pmone * (halfLengthSiriHolder-h/2.)),
                              "SiRiHolder", 
                              logSiriHolder, 
                              physiMother,
                              false,
                              0, 
                              checkOverlaps);


   physCableCu = new G4PVPlacement(0,  
                              G4ThreeVector(0,0,0),
                              logCableCu, 
                              "SiRiCableCu", 
                              logCable,
                              false,
                              0, 
                              checkOverlaps);

  // place the pads
  for (G4int ipad = 0; ipad < nb_pads ; ipad++) {
    G4double phi = ipad*dPhi;
    G4double right_angle = 90*deg;
    G4double rot_angle = phi - pmone * right_angle;			// DO NOT CHANGE UNDER ANY CIRCUMSTANCES!

    G4RotationMatrix rotm  = G4RotationMatrix();		// global rotation matrix
    rotm.rotateX(anotherAngle); // rotation along theta
    rotm.rotateZ(rot_angle); // rotation along phi
	  
    G4ThreeVector uz = G4ThreeVector(std::cos(phi),
                                     std::sin(phi),
                                     0.);

    G4ThreeVector posE = R_E*uz;
    posE += translatePos;                      // translate the origin to have the source centred
    
                            
   new G4PVPlacement(G4Transform3D(rotm,posE), "SiRi_pad", pad_Log, physiMother, false, ipad, checkOverlaps);
   
   // Cables'
    uz = G4ThreeVector(std::cos(phi),
                       std::sin(phi),
                       0);
    posE = R_E*uz;
    posE += translatePos;
    posE += G4ThreeVector(std::cos(phi),
                       std::sin(phi),
                       pmone*(halfLengthCable+zTransCable));
    rotm  = G4RotationMatrix();    // global rotation matrix
    // rotm.rotateX(anotherAngle); // rotation along theta
    rotm.rotateZ(rot_angle); // rotation along phi
   new G4PVPlacement(G4Transform3D(rotm,posE), "SiRiCable", logCable, physiMother, false, ipad, checkOverlaps);
  } //ipad-loop

  

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void SiRi::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

// declare trackers as a MultiFunctionalDetector scorer
  // 1 
  // G4MultiFunctionalDetector* Si_deltaE1 = new G4MultiFunctionalDetector("Si_thin1");		// Delta E detector
  // G4VPrimitiveScorer* primitiv1 = new G4PSEnergyDeposit("edep1_1");
  // Si_deltaE1->RegisterPrimitive(primitiv1);
  // SetSensitiveDetector("ThinLV1",Si_deltaE1);
  // // 2
  // G4MultiFunctionalDetector* Si_deltaE2 = new G4MultiFunctionalDetector("Si_thin2");		// Delta E detector
  // G4VPrimitiveScorer* primitiv2 = new G4PSEnergyDeposit("edep1_2");
  // Si_deltaE2->RegisterPrimitive(primitiv2);
  // SetSensitiveDetector("ThinLV2",Si_deltaE2);
  // // 3
  // G4MultiFunctionalDetector* Si_deltaE3 = new G4MultiFunctionalDetector("Si_thin3");		// Delta E detector
  // G4VPrimitiveScorer* primitiv3 = new G4PSEnergyDeposit("edep1_3");
  // Si_deltaE3->RegisterPrimitive(primitiv3);
  // SetSensitiveDetector("ThinLV3",Si_deltaE3);
  // // 4
  // G4MultiFunctionalDetector* Si_deltaE4 = new G4MultiFunctionalDetector("Si_thin4");		// Delta E detector
  // G4VPrimitiveScorer* primitiv4 = new G4PSEnergyDeposit("edep1_4");
  // Si_deltaE4->RegisterPrimitive(primitiv4);
  // SetSensitiveDetector("ThinLV4",Si_deltaE4);
  // // 5
  // G4MultiFunctionalDetector* Si_deltaE5 = new G4MultiFunctionalDetector("Si_thin5");		// Delta E detector
  // G4VPrimitiveScorer* primitiv5 = new G4PSEnergyDeposit("edep1_5");
  // Si_deltaE5->RegisterPrimitive(primitiv5);
  // SetSensitiveDetector("ThinLV5",Si_deltaE5);
  // // 6
  // G4MultiFunctionalDetector* Si_deltaE6 = new G4MultiFunctionalDetector("Si_thin6");		// Delta E detector
  // G4VPrimitiveScorer* primitiv6 = new G4PSEnergyDeposit("edep1_6");
  // Si_deltaE6->RegisterPrimitive(primitiv6);
  // SetSensitiveDetector("ThinLV6",Si_deltaE6);
  // // 7
  // G4MultiFunctionalDetector* Si_deltaE7 = new G4MultiFunctionalDetector("Si_thin7");		// Delta E detector
  // G4VPrimitiveScorer* primitiv7 = new G4PSEnergyDeposit("edep1_7");
  // Si_deltaE7->RegisterPrimitive(primitiv7);
  // SetSensitiveDetector("ThinLV7",Si_deltaE7);
  // // 8 
  // G4MultiFunctionalDetector* Si_deltaE8 = new G4MultiFunctionalDetector("Si_thin8");		// Delta E detector
  // G4VPrimitiveScorer* primitiv8 = new G4PSEnergyDeposit("edep1_8");
  // Si_deltaE8->RegisterPrimitive(primitiv8);
  // SetSensitiveDetector("ThinLV8",Si_deltaE8);
  
  //.......................
  G4MultiFunctionalDetector* Si_E = new G4MultiFunctionalDetector("Si_thick");		// E detector
  G4VPrimitiveScorer* primitiv9 = new G4PSEnergyDeposit("edep2");
  Si_E->RegisterPrimitive(primitiv9);
  SetSensitiveDetector("ThickLV",Si_E);

}

  
                      

