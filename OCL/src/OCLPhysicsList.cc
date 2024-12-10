//
//
// Physics List is a simplification of the LXePhysicsList 
// ($G4INSTALL/examples/extended/optical/LXe). EM physics 
// just registering G4EmStandardPhysics and no Decay Physics.
//
////////////////////////////////////////////////////////

#include "OCLPhysicsList.hh"

#include "G4EmStandardPhysics.hh"
#include "G4OpticalPhysics.hh"

#include "G4ParticleTypes.hh"

#include "G4SystemOfUnits.hh"


OCLPhysicsList::OCLPhysicsList() : G4VModularPhysicsList()
{
	// default cut value  (0.1 mm)
	defaultCutValue = 0.1*mm;
	
	// EM Physics
	RegisterPhysics(new G4EmStandardPhysics());
	
	/*
	// this whole block can be commented to speed up the simulation
	// Optical Physics
	G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
	RegisterPhysics(opticalPhysics);


	//opticalPhysics->SetScintillationYieldFactor(1);
	// "articifially" reduce the yield in order to create/display less opt. photons
	opticalPhysics->SetScintillationYieldFactor(0.008);

    opticalPhysics->SetScintillationExcitationRatio(0.);

	opticalPhysics->SetTrackSecondariesFirst(kScintillation,true);
	*/
}


OCLPhysicsList::~OCLPhysicsList() {}

void OCLPhysicsList::SetCuts(){
	
	SetCutsWithDefault();
}
