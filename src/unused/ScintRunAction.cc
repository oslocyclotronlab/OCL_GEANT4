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
// $Id: ScintRunAction.cc 71323 2013-06-13 16:54:23Z gcosmo $
//
/// \file ScintRunAction.cc
/// \brief Implementation of the ScintRunAction class

#include "ScintRunAction.hh"
#include "ScintPrimaryGeneratorAction.hh"
#include "ScintRun.hh"

#include "Analysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ScintRunAction::ScintRunAction()
 : G4UserRunAction()
{  
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);     

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in B4Analysis.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories 
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFirstHistoId(1);

  // Book histograms, ntuple
  //
  
  // Creating histograms
  analysisManager->CreateH2("hEdep","Edep", 1000, 0, 1000, 1000, 0., 20.*MeV);	// h2 Id = 1
  analysisManager->SetH2XAxisTitle(1, "Event Number");
  analysisManager->SetH2YAxisTitle(1, "Dep. E [MeV]");
  
  analysisManager->CreateH1("hdEdN","dEdN", 100,0.,20.*MeV);		// h1 Id = 1
  analysisManager->SetH1XAxisTitle(1, "Dep. E [MeV]");
  analysisManager->SetH1YAxisTitle(1, "dE/dN");


  // Creating ntuple
  //
  analysisManager->CreateNtuple("Scint", "Energy");
  analysisManager->CreateNtupleDColumn("edep1");		// column Id = 0
  analysisManager->CreateNtupleIColumn("evtNb");		// column Id = 1
  analysisManager->FinishNtuple();   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ScintRunAction::~ScintRunAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* ScintRunAction::GenerateRun()
{ return new ScintRun; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ScintRunAction::BeginOfRunAction(const G4Run* run)
{ 
  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;
  
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  
  
   // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  //
  G4String fileName = "output";
  analysisManager->OpenFile(fileName);
  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ScintRunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;
  
  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const ScintPrimaryGeneratorAction* generatorAction
    = static_cast<const ScintPrimaryGeneratorAction*>(
        G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String partName;
  if (generatorAction) 
  {
    G4ParticleDefinition* particle 
      = generatorAction->GetParticleGun()->GetParticleDefinition();
    partName = particle->GetParticleName();
  }  
  
  
  // save histograms & ntuple
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  analysisManager->Write();
  analysisManager->CloseFile();
  
  //results
  //
  const ScintRun* scRun = static_cast<const ScintRun*>(run);
  G4int nbGoodEvents = scRun->GetNbGoodEvents();

  {
    G4cout
     << "\n--------------------End of Local Run------------------------"
     << " \n The run was " << nofEvents << " "<< partName;
  }      
  G4cout
     << "\n Hits in both detectors: " << nbGoodEvents
     << "\n------------------------------------------------------------\n"
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
