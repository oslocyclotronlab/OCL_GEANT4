/////////////////////////////////////////////////////////////////////////
//
//
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

#include "OCLSteppingAction.hh"
#include "OCLEventAction.hh"
#include "G4SteppingManager.hh"
#include "OCLAnalysis.hh"

#include "G4RunManager.hh"

#include "fstream"
#include "iomanip"

using namespace std;	 

OCLSteppingAction::OCLSteppingAction(OCLEventAction* EvAct)
:eventAction(EvAct)
{ }

OCLSteppingAction::~OCLSteppingAction()
{ }

void OCLSteppingAction::UserSteppingAction(const G4Step* aStep)
{

	G4StepPoint* point = aStep->GetPreStepPoint();
	G4TouchableHandle touch = point->GetTouchableHandle();

	const G4String currentPhysicalName 
    = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
	
	const G4String particleName
	= aStep->GetTrack()->GetDefinition()->GetParticleName();

	G4int copyNumber;
	G4int depth = 5; 						// in caes it really is the DetectorVolume
	
	if (currentPhysicalName == "Crystal"){
		G4String GrandMotherPhysicalName = touch->GetVolume(depth)->GetName();
		
		G4double EdepStep = aStep->GetTotalEnergyDeposit();
		copyNumber = touch->GetCopyNumber(depth); // of depth=5 mother: "OCLDetector"
		
		if (EdepStep > 0. &&GrandMotherPhysicalName =="OCLDetector") {

			eventAction->EdepInCrystalTest[copyNumber]+=EdepStep;							}

		//count scintillating photons and kill the photons after the first step
		if (particleName == "opticalphoton"){
			eventAction->nAbsPhotons++;
			eventAction->absTime = aStep -> GetPreStepPoint() -> GetGlobalTime();
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
		}


	}
	
	// check if the photon is absorbed in the sensitive volume
	if (currentPhysicalName == "Cathode"){
		const G4String ProcessName = 
		aStep -> GetPostStepPoint() -> GetProcessDefinedStep() -> GetProcessName();
		if (ProcessName == "OpAbsorption"){ 
			
			eventAction->nAbsPhotons++;
			
			eventAction->absTime = aStep -> GetPreStepPoint() -> GetGlobalTime();
		} 
	}
	  
}

