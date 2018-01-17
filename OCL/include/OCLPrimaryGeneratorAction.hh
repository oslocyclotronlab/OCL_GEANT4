#ifndef OCLPrimaryGeneratorAction_h
#define OCLPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4GeneralParticleSource;
class G4Event;

class OCLPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    OCLPrimaryGeneratorAction();
    ~OCLPrimaryGeneratorAction();
	
public:
    void GeneratePrimaries(G4Event* anEvent);
	
private:
    G4GeneralParticleSource* particleGun;
	
};

#endif
