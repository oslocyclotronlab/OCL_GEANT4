/////////////////////////////////////////////////////////////////////////
//
// Oct/2013  E. Nacher -> SteppingAction.cc
//
// This class collects information at three different levels: 
// 1 - Energy deposited in the sensitive volume (LaBr3 crystal)
// 2 - Number of photons absorbed by the photocathode (or just generated!)
// 3 - Time at which each photon is absorbed
//
// This information is passed to the Event Action Class via 
// the eventAction pointer
//
/////////////////////////////////////////////////////////////////////////

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "G4SteppingManager.hh"

#include "fstream"
#include "iomanip"

using namespace std;	 

SteppingAction::SteppingAction(EventAction* EvAct) 
:eventAction(EvAct)
{ }

SteppingAction::~SteppingAction()
{ }

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
	const G4String currentPhysicalName 
    = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
	
	const G4String particleName
	= aStep->GetTrack()->GetDefinition()->GetParticleName();
	
	if (currentPhysicalName == "Crystal"){
		
//        if (particleName == "gamma" )
//        if (particleName == "gamma" || "opticalphoton" )
		if (particleName == "opticalphoton" )
		{
		}
		else 
		{
		G4double EdepStep = aStep->GetTotalEnergyDeposit();
		if (EdepStep > 0.) {
		G4double eventAction->totEnergyDep = eventAction->totEnergyDep + EdepStep;
		//G4cout <<  eventAction->totEnergyDep << endl;
		}
		}
		
//		//count scintillating photons and kill the photons after the first step
//		if (particleName == "opticalphoton"){
//			eventAction->nAbsPhotons++;
//			G4double absTime = aStep -> GetPreStepPoint() -> GetGlobalTime();
//			eventAction-> fillThistogram(absTime);
//			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
//		}
	}
	
	// check if the photon is absorbed in the sensitive volume
	if (currentPhysicalName == "Cathode"){
		const G4String ProcessName = 
		aStep -> GetPostStepPoint() -> GetProcessDefinedStep() -> GetProcessName();
		if (ProcessName == "OpAbsorption"){ 
			
			// get number of absorbed photons
			eventAction->nAbsPhotons++;
			
			// get time
			G4double absTime = aStep -> GetPreStepPoint() -> GetGlobalTime();
			eventAction-> fillThistogram(absTime);
			
//			// get deposited energy
//			G4double EdepStepCathode = aStep->GetTotalEnergyDeposit();
//		    if (EdepStep > 0.) eventAction->totEnergyDepCathod 
//		                                        = eventAction->totEnergyDepCathod + EdepStepCathod;
		} 
	}
}


