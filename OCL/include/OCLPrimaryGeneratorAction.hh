#ifndef SingleScintPrimaryGeneratorAction_h
#define SingleScintPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4GeneralParticleSource;
class G4Event;

class SingleScintPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    SingleScintPrimaryGeneratorAction();
    ~SingleScintPrimaryGeneratorAction();
	
public:
    void GeneratePrimaries(G4Event* anEvent);
	
private:
    G4GeneralParticleSource* particleGun;
	
};

#endif
