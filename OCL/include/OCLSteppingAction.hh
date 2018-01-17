
#ifndef SingleScintSteppingAction_h
#define SingleScintSteppingAction_h 1

#include "G4UserSteppingAction.hh"

class SingleScintEventAction;
// class RunAction;

class SingleScintSteppingAction : public G4UserSteppingAction
{
  public:
    SingleScintSteppingAction(SingleScintEventAction*);
    ~SingleScintSteppingAction();

    void UserSteppingAction(const G4Step*);
    
  private:
    SingleScintEventAction* eventAction;

};

#endif
