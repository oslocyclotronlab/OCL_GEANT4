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

#include "SingleScintSteppingAction.hh"
#include "SingleScintEventAction.hh"
#include "G4SteppingManager.hh"
#include "SingleScintAnalysis.hh"

#include "G4RunManager.hh"

#include "fstream"
#include "iomanip"

using namespace std;	 

SingleScintSteppingAction::SingleScintSteppingAction(SingleScintEventAction* EvAct)
:eventAction(EvAct)
{ }

SingleScintSteppingAction::~SingleScintSteppingAction()
{ }

void SingleScintSteppingAction::UserSteppingAction(const G4Step* aStep)
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

			// if (copyNumber==0)eventAction->EdepInCrystal1  = eventAction->EdepInCrystal1 + EdepStep;
			// if (copyNumber==1)eventAction->EdepInCrystal2  = eventAction->EdepInCrystal2 + EdepStep;
			// if (copyNumber==2)eventAction->EdepInCrystal3  = eventAction->EdepInCrystal3 + EdepStep;
			// if (copyNumber==3)eventAction->EdepInCrystal4  = eventAction->EdepInCrystal4 + EdepStep;
			// if (copyNumber==4)eventAction->EdepInCrystal5  = eventAction->EdepInCrystal5 + EdepStep;
			// if (copyNumber==5)eventAction->EdepInCrystal6  = eventAction->EdepInCrystal6 + EdepStep;
			// if (copyNumber==6)eventAction->EdepInCrystal7  = eventAction->EdepInCrystal7 + EdepStep;
			// if (copyNumber==7)eventAction->EdepInCrystal8  = eventAction->EdepInCrystal8 + EdepStep;
			// if (copyNumber==8)eventAction->EdepInCrystal9  = eventAction->EdepInCrystal9 + EdepStep;
			// if (copyNumber==9)eventAction->EdepInCrystal10 = eventAction->EdepInCrystal10 + EdepStep;
			// if (copyNumber==10)eventAction->EdepInCrystal11 = eventAction->EdepInCrystal11 + EdepStep;
			// if (copyNumber==11)eventAction->EdepInCrystal12 = eventAction->EdepInCrystal12 + EdepStep;
			// if (copyNumber==12)eventAction->EdepInCrystal13 = eventAction->EdepInCrystal13 + EdepStep;
			// if (copyNumber==13)eventAction->EdepInCrystal14 = eventAction->EdepInCrystal14 + EdepStep;
			// if (copyNumber==14)eventAction->EdepInCrystal15 = eventAction->EdepInCrystal15 + EdepStep;
			// if (copyNumber==15)eventAction->EdepInCrystal16 = eventAction->EdepInCrystal16 + EdepStep;
			// if (copyNumber==16)eventAction->EdepInCrystal17 = eventAction->EdepInCrystal17 + EdepStep;
			// if (copyNumber==17)eventAction->EdepInCrystal18 = eventAction->EdepInCrystal18 + EdepStep;
			// if (copyNumber==18)eventAction->EdepInCrystal19 = eventAction->EdepInCrystal19 + EdepStep;
			// if (copyNumber==19)eventAction->EdepInCrystal20 = eventAction->EdepInCrystal20 + EdepStep;
			// if (copyNumber==20)eventAction->EdepInCrystal21 = eventAction->EdepInCrystal21 + EdepStep;
			// if (copyNumber==21)eventAction->EdepInCrystal22 = eventAction->EdepInCrystal22 + EdepStep;
			// if (copyNumber==22)eventAction->EdepInCrystal23 = eventAction->EdepInCrystal23 + EdepStep;
			// if (copyNumber==23)eventAction->EdepInCrystal24 = eventAction->EdepInCrystal24 + EdepStep;
			// if (copyNumber==24)eventAction->EdepInCrystal25 = eventAction->EdepInCrystal25 + EdepStep;
			// if (copyNumber==25)eventAction->EdepInCrystal26 = eventAction->EdepInCrystal26 + EdepStep;
			// if (copyNumber==26)eventAction->EdepInCrystal27 = eventAction->EdepInCrystal27 + EdepStep;
			// if (copyNumber==27)eventAction->EdepInCrystal28 = eventAction->EdepInCrystal28 + EdepStep;
			// if (copyNumber==28)eventAction->EdepInCrystal29 = eventAction->EdepInCrystal29 + EdepStep;
			// if (copyNumber==29)eventAction->EdepInCrystal30 = eventAction->EdepInCrystal30 + EdepStep;
			// if (copyNumber==30)eventAction->EdepInCrystal31 = eventAction->EdepInCrystal31 + EdepStep;
			// if (copyNumber==31)eventAction->EdepInCrystal32 = eventAction->EdepInCrystal32 + EdepStep;
			eventAction->EdepInCrystalTest[copyNumber]+=EdepStep;

							}

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

