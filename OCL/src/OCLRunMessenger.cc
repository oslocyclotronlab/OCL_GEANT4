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

#include "OCLRunMessenger.hh"

#include "OCLRunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4Material.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OCLRunMessenger::OCLRunMessenger(OCLRunAction* run)
 : G4UImessenger(),
   fRun(run),
   fDirectory(0),
   fOutNameCmd(0)
{
  fDirectory = new G4UIdirectory("/OCL/");
  fDirectory->SetGuidance("UI commands of OCL (OSCAR)");

  fOutNameCmd = new G4UIcmdWithAString("/OCL/setOutName",this);
  fOutNameCmd->SetGuidance("Select output root (directory and) file name.");
  fOutNameCmd->SetParameterName("choice",false);
  fOutNameCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OCLRunMessenger::~OCLRunMessenger()
{
  delete fDirectory;
  delete fOutNameCmd;
}

void OCLRunMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if( command == fOutNameCmd ) {
    fRun->SetOutName(newValue);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String OCLRunMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String ans;
  if( command == fOutNameCmd ){
    ans=fRun->GetOutName();
  }
  return ans;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
