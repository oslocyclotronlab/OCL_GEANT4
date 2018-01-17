///////////////////////////////////////////////////////////////////////////////////////
//
// Apr/2015  E. Nacher --> SingleScintPrimaryGeneratorAction.cc
//
// Based on the /gps method. This simplifies life, but...
// The SingleScintPrimaryGeneratorAction must be instantiated after initialization of the
// runManager in the main.cc:  
//                          runManager->Initialize();
//                          runManager->SetUserAction(new SingleScintPrimaryGeneratorAction);
//
///////////////////////////////////////////////////////////////////////////////////////

#include "SingleScintPrimaryGeneratorAction.hh"
#include "SingleScintParameters.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "G4ThreeVector.hh"
#include "globals.hh"

#include "G4ios.hh"
#include "fstream"
#include "iomanip"
#include "G4GeneralParticleSource.hh" 

SingleScintPrimaryGeneratorAction::SingleScintPrimaryGeneratorAction()
{
	
	// Default values  
	
	particleGun = new G4GeneralParticleSource();
	particleGun->SetCurrentSourceIntensity (1);
	particleGun->SetParticlePosition(G4ThreeVector());
  	G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
  	particleGun->SetParticleDefinition(particleDefinition);

	// 	// Source position determined from Parameters.hh
	// 	particleGun->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(0., 0., - distSourceHalfCry ));
}

SingleScintPrimaryGeneratorAction::~SingleScintPrimaryGeneratorAction()
{
	delete particleGun;
}

void SingleScintPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	//create vertex
	
	particleGun->GeneratePrimaryVertex(anEvent);
}
