/*

	NIFF module with simplified geometry.
	Material of the module is iso-butane, surrounded by vacuum.

*/

#include "NIFF.hh"
#include "OCLMaterials.hh"



NIFF::NIFF()
: G4VUserDetectorConstruction()

{

  //----------------------------------------------------
  // Material definitions
  //----------------------------------------------------

  OCLMaterials* fMat = OCLMaterials::GetInstance();
  aluminium = fMat->GetMaterial("G4_Al");
  isobutane = fMat->GetMaterial("isoC4H10");
  isobutane_ppac = fMat->GetMaterial("isoC4H10_PPAC");
  Mylar = fMat->GetMaterial("Mylar");
  vacuum = fMat->GetMaterial("Vacuum");

  /*

    The density of vacuum cannot be 0, hence we use some very small number.
    Temperature and pressure needs to be specified for correct dE/dx calculation in case of non-STP gases.

  */

  //
  // Create the solids.....
  //

  CreateSolids();


}


NIFF::~NIFF()
{ }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Dummy class -- doesn't compile otherwise
G4VPhysicalVolume* NIFF::Construct()
{
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// void NIFF::SetPosition(G4ThreeVector thisPos) {
//   translatePos = thisPos*mm;
//   G4cout << " ----> NIFF will be placed at distance: " << translatePos/mm << " mm" << G4endl;
// }




// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// void NIFF::SetRotation(G4RotationMatrix thisRot) {
//   rotation = thisRot;
// }



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void  NIFF::CreateSolids() {


  // Define the detector
  /*
  	We define the NIFF module as a hollow trapezoid (disjunction of 2 trapezoids).
  	Dimensions of the outer shell:
  						bottom base x1 = 77 mm
  						upper base  x2 = 20 mm
  						side length  b = 55 mm
  	Thickness of module:   5.5 mm

  	G4 defines trapezoid as: (bottom_lenght_X, upper_length_X, bottom_length_Y, upper_length_Y, height_Z).
  	Hence we need to first compute the height of our trapezoid.
  */

  pos1 = G4ThreeVector(0, 0, 0);					// position of the shape
  // G4RotationMatrix rMat;									// rotation matrix around X axis
  G4double phiX = 180. *deg;
  G4double phiZ = 45. *deg;
  //  	rMat.rotateX(phiX);
  rMat.rotateZ(phiZ);


  // set dimensions of outer trapezoid
  G4double outer_b = 60.*mm;
  G4double outer_x1 = 77.*mm;
  G4double outer_x2 = 20.*mm;
  G4double thickness = 5.5*mm;

  // compute dimensions of inner trapezoid
  G4double inner_x1 = outer_x1 - thickness;		// in mm
  G4double inner_x2 = outer_x2 - thickness;		// in mm

  // compute height
  G4double height = sqrt(pow(outer_b,2.) - 2.*pow((outer_x1 - outer_x2)/2.,2.));

  // define outer trapezoid, all dimensions have to be set as half-lengths!
  G4Trd* outer_trap = new G4Trd("Outer", .5*outer_x1, .5*outer_x2, .5*outer_x1, .5*outer_x2, .5*height);

  // define inner trapezoid, all dimensions have to be set as half-lengths!
  G4Trd* inner_trap = new G4Trd("Inner", .5*inner_x1, .5*inner_x2, .5*inner_x1, .5*inner_x2, .5*height);

  // subtract the two solids
  G4SubtractionSolid *trapezoid = new G4SubtractionSolid("Outer-Inner", outer_trap, inner_trap, 0, pos1);

  // add alufoil (10 micron thick), height is the same as for the gas part
  G4double out_alu_x1 = inner_x1;
  G4double out_alu_x2 = inner_x2;
  G4double thick_alu = 1.e-2*mm;

  G4double in_alu_x1 = out_alu_x1 - thick_alu;
  G4double in_alu_x2 = out_alu_x2 - thick_alu;

  G4Trd* out_trap_alu = new G4Trd("Out_Al", .5*out_alu_x1, .5*out_alu_x2, .5*out_alu_x1, .5*out_alu_x2, .5*height);
  G4Trd* in_trap_alu = new G4Trd("In_Al", .5*in_alu_x1, .5*in_alu_x2, .5*in_alu_x1, .5*in_alu_x2, .5*height);
  G4SubtractionSolid *trap_alu = new G4SubtractionSolid("Out-In_Al", out_trap_alu, in_trap_alu, 0, pos1);

  // add Mylar foil (1.5 micron), height is the same as for the gas part
  G4double out_Mylar_x1 = in_alu_x1;
  G4double out_Mylar_x2 = in_alu_x2;
  G4double thick_Mylar = 1.5e-3*mm;

  G4double in_Mylar_x1 = out_Mylar_x1 - thick_Mylar;
  G4double in_Mylar_x2 = out_Mylar_x2 - thick_Mylar;

  G4Trd* out_trap_Mylar = new G4Trd("Out_Mylar", .5*out_Mylar_x1, .5*out_Mylar_x2, .5*out_Mylar_x1, .5*out_Mylar_x2, .5*height);
  G4Trd* in_trap_Mylar = new G4Trd("In_Mylar", .5*in_Mylar_x1, .5*in_Mylar_x2, .5*in_Mylar_x1, .5*in_Mylar_x2, .5*height);
  G4SubtractionSolid *trap_Mylar = new G4SubtractionSolid("Out-In_Mylar", out_trap_Mylar, in_trap_Mylar, 0, pos1);

  // create the module
  trapLog = new G4LogicalVolume(trapezoid, isobutane, "NIFF");

  // Alu and Mylar layers
  trapAlLog = new G4LogicalVolume(trap_alu, aluminium, "NIFF_Al");
  trapMylarLog = new G4LogicalVolume(trap_Mylar, Mylar, "NIFF_Mylar");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void NIFF::Placement(G4int copyNo, G4VPhysicalVolume* physiMother, G4bool checkOverlaps)
{

  // // NIFF module
  new G4PVPlacement(G4Transform3D(rMat,pos1), "NIFF", trapLog, physiMother, copyNo, checkOverlaps);    // defines the position, rotation
  fScoringVolume = trapLog; // sets the module as scoring volume

  // add alu and Mylar layers
  new G4PVPlacement(G4Transform3D(rMat,pos1), "NIFF_Al",   trapAlLog,    physiMother, copyNo, checkOverlaps);
  new G4PVPlacement(G4Transform3D(rMat,pos1), "NIFF_Mylar",trapMylarLog, physiMother, copyNo, checkOverlaps);



  // //------------------------------------------------------
  // // visualization attributes
  // //------------------------------------------------------

  G4VisAttributes* worldVisAtt1 = new G4VisAttributes(G4Colour(1.0,0.0,0.0)); //red
  worldVisAtt1->SetVisibility(true);
  trapAlLog->SetVisAttributes(worldVisAtt1);

  G4VisAttributes* worldVisAtt2 = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  worldVisAtt2->SetVisibility(true);
  trapMylarLog->SetVisAttributes(worldVisAtt2);


}



