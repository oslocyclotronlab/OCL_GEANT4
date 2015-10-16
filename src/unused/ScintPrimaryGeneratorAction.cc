#include "ScintPrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ScintPrimaryGeneratorAction::ScintPrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="gamma");
  fParticleGun->SetParticleDefinition(particle);
//  fParticleGun->SetParticleEnergy(15.*MeV);			// use if you want all particles to have the same energy
  fParticleGun->SetParticleEnergy(170.*keV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ScintPrimaryGeneratorAction::~ScintPrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ScintPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
//	// variable origin
//  G4double size = 10.*mm;
//  G4double x0 = size*(G4UniformRand()-0.5);
//  G4double y0 = size*(G4UniformRand()-0.5);
//  G4double z0 = -5.*cm;

	// fixed origin
  G4double x0 = 0.;
  G4double y0 = 0.;
  G4double z0 = -25.*cm;
  
  
  // fixed momentum vector
  G4double px = 0.;
  G4double py = 0.;
  G4double pz = 1.;

//  // use the settings bellow to get 'isotropically' distributed projectiles emitted from a point-like source
//  
//  G4double phi = 2*pi*G4UniformRand();
//  G4double theta = (pi/2)*G4UniformRand();
//  G4double px = sin(phi)*sin(theta);
//  G4double py = cos(phi)*sin(theta);
//  G4double pz = cos(theta);


//	// VARY ENERGY
//	G4double Emax = 10.*MeV;
//	G4double energy = Emax*G4UniformRand();
//	fParticleGun->SetParticleEnergy(energy);
//	
  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px,py,pz));

  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

