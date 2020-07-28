///////////////////////////////////////////////////////////////////////////////////////
//
// Apr/2015  E. Nacher --> OCLPrimaryGeneratorAction.cc
//
// Based on the /gps method. This simplifies life, but...
// The OCLPrimaryGeneratorAction must be instantiated after initialization of the
// runManager in the main.cc:
//                          runManager->Initialize();
//                          runManager->SetUserAction(new OCLPrimaryGeneratorAction);
//
///////////////////////////////////////////////////////////////////////////////////////

#include "OCLPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "G4ThreeVector.hh"
#include "globals.hh"

#include "G4ios.hh"
#include "fstream"
#include "iomanip"
#include "G4GeneralParticleSource.hh"

OCLPrimaryGeneratorAction::OCLPrimaryGeneratorAction()
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

OCLPrimaryGeneratorAction::~OCLPrimaryGeneratorAction()
{
	delete particleGun;
}

void OCLPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	//create vertex

	particleGun->GeneratePrimaryVertex(anEvent);
}
