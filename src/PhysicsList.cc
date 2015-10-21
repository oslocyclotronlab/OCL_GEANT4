//
//
// Physics List is a simplification of the LXePhysicsList 
// ($G4INSTALL/examples/extended/optical/LXe). EM physics 
// just registering G4EmStandardPhysics and no Decay Physics.
//
////////////////////////////////////////////////////////

#include "PhysicsList.hh"

#include "G4EmStandardPhysics.hh"
#include "G4OpticalPhysics.hh"
#include "G4OpticalProcessIndex.hh"

#include "G4ParticleTypes.hh"

#include "G4SystemOfUnits.hh"


PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
	// default cut value  (0.1 mm)
	defaultCutValue = 0.1*mm;
	
	// EM Physics
	RegisterPhysics(new G4EmStandardPhysics());
	
	// Optical Physics
	G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
	RegisterPhysics(opticalPhysics);


	//opticalPhysics->SetScintillationYieldFactor(1);
	opticalPhysics->SetScintillationYieldFactor(0.008);

    opticalPhysics->SetScintillationExcitationRatio(0.);

	opticalPhysics->SetTrackSecondariesFirst(kScintillation,true);
	
}


PhysicsList::~PhysicsList() {}

void PhysicsList::SetCuts(){
	
	SetCutsWithDefault();
}
