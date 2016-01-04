 #ifndef SingleScintPhysicsList_h
 #define SingleScintPhysicsList_h 1
 
 #include "G4VModularPhysicsList.hh"
 #include "globals.hh"
 
 class SingleScintPhysicsList: public G4VModularPhysicsList
 {
	   public:
	 
         SingleScintPhysicsList();
         virtual ~SingleScintPhysicsList();
	 
	   public:
	 
	     // SetCuts()
	     virtual void SetCuts();
	 
	 };
 
 #endif


