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

#include "OCLDetectorMessenger.hh"
#include "OCLDetectorConstruction.hh"

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

OCLDetectorMessenger::OCLDetectorMessenger(OCLDetectorConstruction* Det)
 : G4UImessenger(),
   fDetector(Det),
   fDetDirectory(0),
   fLaBrDirectory(0),
   fLaBrDistCmd(0),
   fLaBrUseCmd(0),
   fUseCSGOldTCCmd(0),
   fUseCSGOldTargetCmd(0),
   fUseCSGRadSourceCmd(0),
   fUseCSGSiRiCmd(0),
   fUseCSGNiffCmd(0)
{
  fDetDirectory = new G4UIdirectory("/OCL/det/");
  fDetDirectory->SetGuidance("UI commands of OCL detector (geometry)");
  fLaBrDirectory = new G4UIdirectory("/OCL/det/oscar/");
  fLaBrDirectory->SetGuidance("UI commands of OSCAR");

  fLaBrDistCmd = new G4UIcommand("/OCL/det/oscar/setLaBrDist",this);
  fLaBrDistCmd->SetGuidance("Set source/origen to detector housing face");
  fLaBrDistCmd->SetGuidance("[usage] /OCL/det/oscar/setLaBrDist N distance");
  fLaBrDistCmd->SetGuidance("        N:(int) Internal detector (face) number");
  fLaBrDistCmd->SetGuidance("        distance:(double) Distance (in cm)");
  fLaBrDistCmd->SetGuidance("Note:");
  fLaBrDistCmd->SetGuidance("- N=0 and N=32 are the beamline facets");
  fLaBrDistCmd->SetGuidance("- No geometry check is performed, so you're responsible to check ");
  fLaBrDistCmd->SetGuidance("  that there will not be any geometry overlap");

  G4UIparameter* paramDist = new G4UIparameter("N",'i',false);
  paramDist->SetParameterRange("N>=1&&N<=30");
  fLaBrDistCmd->SetParameter(paramDist);
  paramDist = new G4UIparameter("distance",'d',false);
  paramDist->SetParameterRange("distance>=12"); // just approcimate
  fLaBrDistCmd->SetParameter(paramDist);
  fLaBrDistCmd->AvailableForStates(G4State_PreInit);

  fLaBrUseCmd = new G4UIcommand("/OCL/det/oscar/setLaBrUse",this);
  fLaBrUseCmd->SetGuidance("Set usage of detector number N");
  fLaBrUseCmd->SetGuidance("[usage] /OCL/det/oscar/setLaBrUse N bool");
  fLaBrUseCmd->SetGuidance("        N:(int) Internal detector (face) number");
  fLaBrUseCmd->SetGuidance("        bool:(double) True or false");
  fLaBrUseCmd->SetGuidance("Note:");
  fLaBrUseCmd->SetGuidance("- N=0 and N=32 are the beamline facets");

  G4UIparameter* paramUseL = new G4UIparameter("N",'i',false);
  paramUseL->SetParameterRange("N>=1&&N<=30");
  fLaBrUseCmd->SetParameter(paramUseL);
  paramUseL = new G4UIparameter("bool",'b',false);
  fLaBrUseCmd->SetParameter(paramUseL);
  fLaBrUseCmd->AvailableForStates(G4State_PreInit);

  fUseCSGOldTCCmd = new G4UIcmdWithABool("/OCL/det/useCSGOldTC",this);
  fUseCSGOldTCCmd->SetGuidance("Set usage of old (CSG) Target Chamber");
  fUseCSGOldTCCmd->SetGuidance("Note: Due to the useage of the _HP physics list, "
                               "we have to set this command before initialization");
  fUseCSGOldTCCmd->SetParameterName("use",true);
  fUseCSGOldTCCmd->SetDefaultValue(true);

  fUseCSGOldTargetCmd = new G4UIcmdWithABool("/OCL/det/useCSGOldTarget",this);
  fUseCSGOldTargetCmd->SetGuidance("Set usage of CAD Target");
  fUseCSGOldTargetCmd->SetGuidance("Note: Due to the useage of the _HP physics list, "
                               "we have to set this command before initialization");
  fUseCSGOldTargetCmd->SetParameterName("use",true);
  fUseCSGOldTargetCmd->SetDefaultValue(true);

  fUseCSGRadSourceCmd = new G4UIcmdWithABool("/OCL/det/useCSGRadSource",this);
  fUseCSGRadSourceCmd->SetGuidance("Set usage of CSG radioactive calibration source (rectangular)");
  fUseCSGRadSourceCmd->SetGuidance("Note: Due to the useage of the _HP physics list, "
                               "we have to set this command before initialization");
  fUseCSGRadSourceCmd->SetParameterName("use",true);
  fUseCSGRadSourceCmd->SetDefaultValue(true);

  fUseCSGSiRiCmd = new G4UIcmdWithABool("/OCL/det/useCSGSiRi",this);
  fUseCSGSiRiCmd->SetGuidance("Set usage ofold (CSG) SiRi");
  fUseCSGSiRiCmd->SetGuidance("Note: Due to the useage of the _HP physics list, "
                               "we have to set this command before initialization");
  fUseCSGSiRiCmd->SetParameterName("use",true);
  fUseCSGSiRiCmd->SetDefaultValue(true);

  fUseCSGNiffCmd = new G4UIcmdWithABool("OCL/useCSGNiff",this);
  fUseCSGNiffCmd->SetGuidance("Set usage ofCSG NIFF (PPACS)");
  fUseCSGNiffCmd->SetGuidance("Note: Due to the useage of the _HP physics list, "
                               "we have to set this command before initialization");
  fUseCSGNiffCmd->SetParameterName("use",true);
  fUseCSGNiffCmd->SetDefaultValue(true);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OCLDetectorMessenger::~OCLDetectorMessenger()
{
  delete fDetDirectory;
  delete fLaBrDirectory;
  delete fLaBrDistCmd;
  delete fLaBrUseCmd;
  delete fUseCSGOldTCCmd;
  delete fUseCSGOldTargetCmd;
  delete fUseCSGRadSourceCmd;
  delete fUseCSGSiRiCmd;
  delete fUseCSGNiffCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OCLDetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if( command == fLaBrDistCmd ) {
    G4Tokenizer next( newValue );
    // check argument
    G4int detNumber = StoI(next());
    G4double distance = StoD(next())*cm;
    fDetector->SetLaBr3_Distance(detNumber, distance);
  }

  else if( command == fLaBrUseCmd ) {
    G4Tokenizer next( newValue );
    // check argument
    G4int detNumber = StoI(next());
    G4bool usage = StoB(next());
    fDetector->SetUseLaBr3(detNumber, usage);
  }

  else if( command==fUseCSGOldTCCmd )
  { fDetector->SetUseCSGOldTargetChamber(fUseCSGOldTCCmd->GetNewBoolValue(newValue)); }
  else if( command==fUseCSGOldTargetCmd )
  { fDetector->SetUseCSGOldTarget(fUseCSGOldTargetCmd->GetNewBoolValue(newValue)); }
  else if( command==fUseCSGRadSourceCmd )
  { fDetector->SetUseCSGRadSource(fUseCSGRadSourceCmd->GetNewBoolValue(newValue)); }
  else if( command==fUseCSGSiRiCmd )
  { fDetector->SetUseCSGSiRi(fUseCSGSiRiCmd->GetNewBoolValue(newValue)); }
  else if( command==fUseCSGNiffCmd )
  { fDetector->SetUseCSGNiff(fUseCSGNiffCmd->GetNewBoolValue(newValue)); }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// G4String OCLDetectorMessenger::GetCurrentValue(G4UIcommand * command)
// {
//   G4String ans;
//   if( command == fUse ){
//     ans=fRun->GetOutName();
//   }
//   return ans;
// }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
