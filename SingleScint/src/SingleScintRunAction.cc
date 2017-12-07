///////////////////////////////////////////
//
// Oct/2015  Fabio -> RunAction.cc
//
///////////////////////////////////////////

#include "SingleScintRunAction.hh"
#include "SingleScintAnalysis.hh"

#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4Run.hh"
#include "G4ios.hh"

#include <iomanip>

SingleScintRunAction::SingleScintRunAction() : G4UserRunAction()
{
	  // Create analysis manager
	  // The choice of analysis technology is done via selecting of a namespace
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

	  //  G4int CreateH1(const G4String& name, const G4String& title,
	  //                 G4int nbins, G4double xmin, G4double xmax,
	  //                 const G4String& unitName = "none",
	  //                 const G4String& fcnName = "none",
	  //                 const G4String& binSchemeName = "linear");

	  xmin = 0; // in keV
	  xmax = 12e3; // in keV
	  binsize = 2.; // in keV
	  nbins= (int)((xmax-xmin)/binsize);
	  analysisManager->CreateH1("hedep1","Edep in Crystal1", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep2","Edep in Crystal2", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep3","Edep in Crystal3", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep4","Edep in Crystal4", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep5","Edep in Crystal5", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep6","Edep in Crystal6", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep7","Edep in Crystal7", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep8","Edep in Crystal8", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep9","Edep in Crystal9", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep10","Edep in Crystal10", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep11","Edep in Crystal11", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep12","Edep in Crystal12", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep13","Edep in Crystal13", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep14","Edep in Crystal14", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep15","Edep in Crystal15", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep16","Edep in Crystal16", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep17","Edep in Crystal17", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep18","Edep in Crystal18", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep19","Edep in Crystal19", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep20","Edep in Crystal20", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep21","Edep in Crystal21", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep22","Edep in Crystal22", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep23","Edep in Crystal23", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep24","Edep in Crystal24", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep25","Edep in Crystal25", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep26","Edep in Crystal26", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep27","Edep in Crystal27", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep28","Edep in Crystal28", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep29","Edep in Crystal29", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep30","Edep in Crystal30", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep31","Edep in Crystal31", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hedep32","Edep in Crystal32", nbins, xmin*keV, xmax*keV);

	  xmin = 0; //
	  xmax = 12e3; //
	  binsize = 2; //
	  nbins= (int)(xmax-xmin)/binsize;
          analysisManager->CreateH1("Histo2","Absorbed Photons", nbins, xmin, xmax);

      // Here we need some units!
      xmin = 0; // in ns
	  xmax = 500; // in ns
	  binsize = 2; // in ns
	  nbins= (int)(xmax-xmin)/binsize;
          analysisManager->CreateH1("Histo3","Time of Absorption", nbins, xmin, xmax*ns);

 	  xmin = 0; // in keV
	  xmax = 12e3; // in keV
	  binsize = 2.; // in keV
	  nbins= (int)((xmax-xmin)/binsize);
	  analysisManager->CreateH1("Histo4","Folded Edep in Crystal", nbins, xmin*keV, xmax*keV);

	  analysisManager->CreateH1("hefolded1","Edep folded in Crystal1", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded2","Edep folded in Crystal2", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded3","Edep folded in Crystal3", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded4","Edep folded in Crystal4", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded5","Edep folded in Crystal5", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded6","Edep folded in Crystal6", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded7","Edep folded in Crystal7", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded8","Edep folded in Crystal8", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded9","Edep folded in Crystal9", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded10","Edep folded  in Crystal10", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded11","Edep folded  in Crystal11", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded12","Edep folded  in Crystal12", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded13","Edep folded  in Crystal13", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded14","Edep folded  in Crystal14", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded15","Edep folded  in Crystal15", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded16","Edep folded  in Crystal16", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded17","Edep folded  in Crystal17", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded18","Edep folded  in Crystal18", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded19","Edep folded  in Crystal19", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded20","Edep folded  in Crystal20", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded21","Edep folded  in Crystal21", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded22","Edep folded  in Crystal22", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded23","Edep folded  in Crystal23", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded24","Edep folded  in Crystal24", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded25","Edep folded  in Crystal25", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded26","Edep folded  in Crystal26", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded27","Edep folded  in Crystal27", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded28","Edep folded  in Crystal28", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded29","Edep folded  in Crystal29", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded30","Edep folded  in Crystal30", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded31","Edep folded  in Crystal31", nbins, xmin*keV, xmax*keV);
	  analysisManager->CreateH1("hefolded32","Edep folded  in Crystal32", nbins, xmin*keV, xmax*keV);
      

	  // Creating ntuple
	  //
	  analysisManager->CreateNtuple("B4", "Edep and TrackL...");
	  analysisManager->CreateNtupleDColumn("Edep");
	  analysisManager->CreateNtupleDColumn("nAbsPhotons");
	  analysisManager->CreateNtupleDColumn("absTime");
	  analysisManager->CreateNtupleDColumn("EdepFolded");
	  analysisManager->FinishNtuple();

}

SingleScintRunAction::~SingleScintRunAction()
{
	delete G4AnalysisManager::Instance();
}

void SingleScintRunAction::BeginOfRunAction(const G4Run*)
{
	  //inform the runManager to save random number seed
	  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);

	  // Get analysis manager
	  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

	  // Open an output file
	  //
	  G4String fileName = "../data/Edep.root";
	  analysisManager->OpenFile(fileName);
}

void SingleScintRunAction::EndOfRunAction(const G4Run*)
{
	  // print histogram statistics
	  //
	  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	  if ( analysisManager->GetH1(1) ) {
	    G4cout << G4endl << " ----> print histograms statistic ";
	    if(isMaster) {
	      G4cout << "for the entire run " << G4endl << G4endl;
	    }
	    else {
	      G4cout << "for the local thread " << G4endl << G4endl;
	    }

	    G4cout << G4endl << " EAbs : mean = "
	       << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy")
	       << " rms = "
	       << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Energy") << G4endl;

	  }

	  // save histograms & ntuple
	  //
	  analysisManager->Write();
	  analysisManager->CloseFile();

}


