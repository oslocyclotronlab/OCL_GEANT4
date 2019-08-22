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
//
//

#include "OCLLaBr3Messenger.hh"
#include "OCLParameters.hh"

#include "G4Material.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tokenizer.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIdirectory.hh"
#include "G4UnitsTable.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OCLLaBr3Messenger::OCLLaBr3Messenger()
 : G4UImessenger(),
   fDirectory(0),
   fCoatingAlThicknessFrontCmd(0),
   fCoatingAlThicknessCmd(0),
   fReflectorThicknessCmd(0),
   fShieldingHalfThicknessLidCmd(0),
   fReflectorThickness(1.*mm),
   fCoatingAlThicknessFront(2.*mm),
   fCoatingAlThickness(1.*mm),
   fShieldingHalfThicknessLid(defaultshieldingHalfThicknessLid)

{
  fDirectory = new G4UIdirectory("/OCL/det/");
  fDirectory->SetGuidance("UI commands of OCL (OSCAR)");

  fCoatingAlThicknessFrontCmd = new G4UIcmdWithADoubleAndUnit("/OCL/det/setCoatThickFront",this);
  fCoatingAlThicknessFrontCmd->SetGuidance("Set LaBr3 (aluminum) encapsulation thickness in the front.");
  fCoatingAlThicknessFrontCmd->SetParameterName("t",false);
  fCoatingAlThicknessFrontCmd->SetUnitCategory("Length");
  fCoatingAlThicknessFrontCmd->SetRange("t>=0");

  fCoatingAlThicknessCmd = new G4UIcmdWithADoubleAndUnit("/OCL/det/setCoatThickRad",this);
  fCoatingAlThicknessCmd->SetGuidance("Set LaBr3 /aluminum/ encapsulation radial thickness.");
  fCoatingAlThicknessCmd->SetGuidance("Note: There will be a plexi glass layer to fill up the difference");
  fCoatingAlThicknessCmd->SetParameterName("t",false);
  fCoatingAlThicknessCmd->SetUnitCategory("Length");
  fCoatingAlThicknessCmd->SetRange("t>=0 && t<=36.5"); // 1mm reflector + 2.55mm plastic + 1 mm Al

  fReflectorThicknessCmd = new G4UIcmdWithADoubleAndUnit("/OCL/det/setLaBrRefThick",this);
  fReflectorThicknessCmd->SetGuidance("Set LaBr3 reflector thickness.");
  fReflectorThicknessCmd->SetParameterName("t",false);
  fReflectorThicknessCmd->SetUnitCategory("Length");
  fReflectorThicknessCmd->SetRange("t>=0");

  fShieldingHalfThicknessLidCmd = new G4UIcmdWithADoubleAndUnit("/OCL/det/setLaBrLidHalfThick",this);
  fShieldingHalfThicknessLidCmd->SetGuidance("Set LaBr3 /half/-thickness of the lid of the shielding");
  fShieldingHalfThicknessLidCmd->SetParameterName("t",false);
  fShieldingHalfThicknessLidCmd->SetUnitCategory("Length");
  fShieldingHalfThicknessLidCmd->SetRange("t>=0");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OCLLaBr3Messenger::~OCLLaBr3Messenger()
{
  delete fDirectory;
  delete fCoatingAlThicknessFrontCmd;
  delete fCoatingAlThicknessCmd;
  delete fReflectorThicknessCmd;
  delete fShieldingHalfThicknessLidCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OCLLaBr3Messenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if( command == fCoatingAlThicknessFrontCmd )
   { fCoatingAlThicknessFront = fCoatingAlThicknessFrontCmd->GetNewDoubleValue(newValue); }

  if( command == fCoatingAlThicknessCmd )
   { fCoatingAlThickness = fCoatingAlThicknessCmd->GetNewDoubleValue(newValue); }

  if( command == fReflectorThicknessCmd )
   { fReflectorThickness = fReflectorThicknessCmd->GetNewDoubleValue(newValue); }

  if( command == fShieldingHalfThicknessLidCmd )
   { fShieldingHalfThicknessLid =fShieldingHalfThicknessLidCmd->GetNewDoubleValue(newValue); }
}
