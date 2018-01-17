 #ifndef OCLPhysicsList_h
 #define OCLPhysicsList_h 1
 
 #include "G4VModularPhysicsList.hh"
 #include "globals.hh"
 
 class OCLPhysicsList: public G4VModularPhysicsList
 {
	   public:
	 
         OCLPhysicsList();
         virtual ~OCLPhysicsList();
	 
	   public:
	 
	     // SetCuts()
	     virtual void SetCuts();
	 
	 };
 
 #endif


