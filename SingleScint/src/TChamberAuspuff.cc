/*
	TChamberAuspuff module
*/

#include "TChamberAuspuff.hh"

TChamberAuspuff::TChamberAuspuff()
{  
  
  //----------------------------------------------------
  // Material definitions
  //----------------------------------------------------

  G4NistManager* man = G4NistManager::Instance();

  aluminium   = man->FindOrBuildMaterial("G4_Al");
  plexiGlass = man->FindOrBuildMaterial("G4_PLEXIGLASS");
  
  //
  // Create the solids.....
  //

  CreateSolids();
}


TChamberAuspuff::~TChamberAuspuff()
{ }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// void TChamberAuspuff::SetPosition(G4ThreeVector thisPos) {
//   translatePos = thisPos*mm;
//   G4cout << " ----> TChamberAuspuff will be placed at distance: " << translatePos/mm << " mm" << G4endl;
// }


// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// void TChamberAuspuff::SetRotation(G4RotationMatrix thisRot) { 
//   rotation = thisRot; 
// }



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void  TChamberAuspuff::CreateSolids() {

  //
  // Constants
  //

  rInnerTube           = 6.08*cm;
  thicknessMiddleTube  = 3.*mm;
  thicknessOuterTube   = 2.*mm;
  halfLengthTube       = 45.5/2.*cm;

  thicknessEndCaps     = thicknessOuterTube + 8.*mm;
  halfLengthEndCaps    = 1.*cm/2.;

  G4double startPhi = 0*deg;
  G4double stopPhi  = 360*deg;

  const G4int numZPlanes = 7;

  const G4double zPlane[numZPlanes] = {
                             -halfLengthTube - 2.*halfLengthEndCaps,
                             -halfLengthTube,
                             -halfLengthTube - 0.00001*mm, //hotfix for sharpe edged Polycone
                             0.*cm,
                             halfLengthTube  + 0.00001*mm, //hotfix for sharpe edged Polycone
                             halfLengthTube,
                             halfLengthTube + 2.*halfLengthEndCaps
                             };

  const G4double rInner[numZPlanes] = {
                             rInnerTube,
                             rInnerTube,
                             rInnerTube,
                             rInnerTube,
                             rInnerTube,
                             rInnerTube,
                             rInnerTube,
                             };

  const G4double rOuter[numZPlanes] = {
                             rInnerTube + thicknessOuterTube + thicknessEndCaps,
                             rInnerTube + thicknessOuterTube + thicknessEndCaps,
                             rInnerTube + thicknessOuterTube,
                             rInnerTube + thicknessMiddleTube,
                             rInnerTube + thicknessOuterTube,
                             rInnerTube + thicknessOuterTube + thicknessEndCaps,
                             rInnerTube + thicknessOuterTube + thicknessEndCaps
                             };


  //
  // now define shapes and logical volumes	  
  //

  // Main tube

  solidChamberTube = 
        new G4Polycone("ChamberTube",
                       startPhi,
                       stopPhi,
                       numZPlanes,
                       zPlane,
                       rInner,
                       rOuter);

  logicChamberTube = 
        new G4LogicalVolume(solidChamberTube, 
                            aluminium, 
                            "ChamberTube");



}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void TChamberAuspuff::Placement(G4int copyNo, G4VPhysicalVolume* physiMother, G4bool checkOverlaps)
{
  // Main tube

  physiChamberTube = 
    new G4PVPlacement(0,                  // Rotation
                      G4ThreeVector(),    // Transformation (Rot&Transl)
                      "ChamberTube",      // its logical volume
                      logicChamberTube,   // its name
                      physiMother,   // its physical mother volume
                      false,        // unknown "pMany"; def: false
                      0,      // copy number
                      checkOverlaps);   // checkOverlaps

}
                      
