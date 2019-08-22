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

#ifndef OCLLaBr3Messenger_h
#define OCLLaBr3Messenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADoubleAndUnit;

class OCLLaBr3Messenger: public G4UImessenger
{
  public:
    OCLLaBr3Messenger();
    virtual ~OCLLaBr3Messenger();

    virtual void SetNewValue(G4UIcommand*, G4String);

    G4double GetReflectorThickness() { return fReflectorThickness; };
    G4double GetCoatingAlThicknessFront() { return fCoatingAlThicknessFront; };
    G4double GetCoatingAlThickness() { return fCoatingAlThickness; };
    G4double GetShieldingHalfThicknessLid() { return fShieldingHalfThicknessLid; };

  private:
    G4UIdirectory* fDirectory;

    G4UIcmdWithADoubleAndUnit*   fReflectorThicknessCmd;

    G4UIcmdWithADoubleAndUnit*   fCoatingAlThicknessCmd;
    G4UIcmdWithADoubleAndUnit*   fCoatingAlThicknessFrontCmd;

    G4UIcmdWithADoubleAndUnit*   fShieldingHalfThicknessLidCmd;

  private:
    G4double fReflectorThickness;
    G4double fCoatingAlThicknessFront;
    G4double fCoatingAlThickness;
    G4double fShieldingHalfThicknessLid;

};


#endif
