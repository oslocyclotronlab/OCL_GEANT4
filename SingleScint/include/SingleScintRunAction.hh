
#ifndef SingleScintRunAction_h
#define SingleScintRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

#include <fstream>

using namespace std;

class G4Run;

class SingleScintRunAction : public G4UserRunAction
{
  public:
    SingleScintRunAction();
    ~SingleScintRunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);

	G4double xmin;
	G4double xmax;
	G4int binsize;
	G4int nbins;
		
};

#endif

