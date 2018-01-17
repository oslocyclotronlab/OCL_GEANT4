
#ifndef OCLSteppingAction_h
#define OCLSteppingAction_h 1

#include "G4UserSteppingAction.hh"

class OCLEventAction;
// class RunAction;

class OCLSteppingAction : public G4UserSteppingAction
{
  public:
    OCLSteppingAction(OCLEventAction*);
    ~OCLSteppingAction();

    void UserSteppingAction(const G4Step*);
    
  private:
    OCLEventAction* eventAction;

};

#endif
