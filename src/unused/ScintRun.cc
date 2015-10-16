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
// $Id: ScintRun.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file ScintRun.cc
/// \brief Implementation of the ScintRun class

#include "ScintRun.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"

#include "Analysis.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ScintRun::ScintRun()
 : G4Run(), 
   fCollID_scint(-1),
   fPrintModulo(10),
   fGoodEvents(0)
//   fEdep(0.)
{ }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ScintRun::~ScintRun()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ScintRun::RecordEvent(const G4Event* event)
{
  if ( fCollID_scint < 0 ) {
   fCollID_scint 
     = G4SDManager::GetSDMpointer()->GetCollectionID("LaBrScint/edep1");
  
  }

  
  G4int evtNb = event->GetEventID();
  
//  if (evtNb%fPrintModulo == 0) 
  { 
    G4cout << "\n---> end of event: " << evtNb << G4endl;
  }      
  
  //Hits collections
  //  
  G4HCofThisEvent* HCE = event->GetHCofThisEvent();
  if(!HCE) return;
               

   
  G4THitsMap<G4double>* evtMap = 
    static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID_scint));
     
  G4double edep1 = 0.;
            
  std::map<G4int,G4double*>::iterator itr;
  for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++) {
    edep1 = *(itr->second);
  }  


	if(edep1!=0.){
		fGoodEvents++;			// count how many hits particles will pass through both detectors
	}
	


// accumulate statistics
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  // fill histograms
  analysisManager->FillH2(1, evtNb, edep1);
  analysisManager->FillH1(1, edep1);
  
  // fill ntuple
  analysisManager->FillNtupleDColumn(0, edep1);
  analysisManager->FillNtupleIColumn(1, evtNb);
  analysisManager->AddNtupleRow();  
  
  G4Run::RecordEvent(event);      
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ScintRun::Merge(const G4Run* aRun)
{
  const ScintRun* localRun = static_cast<const ScintRun*>(aRun);
  fGoodEvents += localRun->fGoodEvents;


  G4Run::Merge(aRun); 
  
  /* 
  	NOTE: I leave this here only because otherwise I cannot persuade GEANT4 to co-operate. Will try to improve it.
  */
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
