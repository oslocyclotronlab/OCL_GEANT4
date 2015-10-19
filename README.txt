fabiobz 0 2015-10-12 17:08 README.txt


     =========================================================
     Geant4 - Simulation of LaBr3:Ce Scint. Detectors 
     =========================================================

                            Simple Detector
                              -----------

 This simluation implements a single scintillation detector with a PMT. 
 The timing information, energy spectrum and number of absorbed photon per
 event is saved as a histogram. 
 
 //Test of small change
	
 1- GEOMETRY DEFINITION
	
   The geometry is constructed in the DetectorConstruction class.
   The setup consists of a cylinder containing the LaBr3 crystal, and outer ring with
   shielding and a lit in front side of the detector. It is optically coupled to 
   a bialkali photocathode through a Borosilicate PMT Window.
   -- 12/10/15 Currently there is no reflector but the shielding is faked as such.

   The dimensions and materials have been chosen as close as possible to our set-up. Where
   no manufacturer information was available, we used assumptions based e.g. other studies.
   The is an arbitrary mix of materials created with the NIST manager and by hand.
   Additionally, we defined the optical properties and Surfaces and boundary processes for the 
   Scintillation process.
		
 2- PHYSICS LIST
 
  The physic processes are set in the PhysicsList class. Currentøy implemented are
  G4EmStandardPhysics and G4OpticalPhysics (the latter is necessary for scintillation).
  
  One should review/set the ScintillationExcitationRatio, the ratio for the fast&slow excitation
  ratio.
 	
  In addition the build-in interactive command:
               /process/(in)activate processName
  allows to activate/inactivate the processes one by one.
   
 3- AN EVENT : THE PRIMARY GENERATOR
 
  The Primary Generator is defined in the PrimaryGeneratorAction  via 
  the G4GeneralParticleSource. The type of the particle and its energy (and possible 
  biases/shape...) are via macro. Currently run1.mac uses a particle 
  which hits the detector perpendicular to the input face. 
        
 5- DETECTOR RESPONSE

   The PMT response is simulated via UserSteppingAction in the SteppingAction class.
   More precisely, the number and time of absorbed photons in the PMT cathode is recorded;
   the broadening due to the PMT has to be modeled separately and is not implemented here.
   
   The total energy deposited is taken from the crystal volume. As Geant4 does not conserve
   Energy for optical photons, these are ignored in the summation. However, this means one
   still has to carefully check that the crystal reflector is of the right type and all
   scintillation photons are reflected as they should.
   
   This example demonstrates a simple scoring implemented directly
   in the user action classes and B1Run object. 
   Alternative ways of scoring via Geant4 classes can be found in the 
   other examples.
   
   The energy deposited is collected step by step for a selected volume
   in B1SteppingAction and accumulated event by event in B1EventAction.
   
   At end of event, the value acummulated in B1EventAction is added in B1Run
   and summed over the whole run (see B1EventAction::EndOfevent()).
   
   Total dose deposited is computed at B1RunAction::EndOfRunAction(), 
   and printed together with informations about the primary particle.
   
   In multi-threading mode the energy accumulated in B1Run objects per
   workers is merged to the master in B1Run::Merge() and the final
   result is printed on the screen.

   An example of creating and computing new units (e.g., dose) is also shown 
   in the class constructor. 

 The following paragraphs are common to all basic examples

 A- VISUALISATION

   The visualization manager is set via the G4VisExecutive class
   in the main() function in exampleB1.cc.    
   The initialisation of the drawing is done via a set of /vis/ commands
   in the macro vis.mac. This macro is automatically read from
   the main function when the example is used in interactive running mode.

   By default, vis.mac opens an OpenGL viewer (/vis/open OGL).
   The user can change the initial viewer by commenting out this line
   and instead uncommenting one of the other /vis/open statements, such as
   HepRepFile or DAWNFILE (which produce files that can be viewed with the
   HepRApp and DAWN viewers, respectively).  Note that one can always
   open new viewers at any time from the command line.  For example, if
   you already have a view in, say, an OpenGL window with a name
   "viewer-0", then
      /vis/open DAWNFILE
   then to get the same view
      /vis/viewer/copyView viewer-0
   or to get the same view *plus* scene-modifications
      /vis/viewer/set/all viewer-0
   then to see the result
      /vis/viewer/flush

   The DAWNFILE, HepRepFile drivers are always available
   (since they require no external libraries), but the OGL driver requires
   that the Geant4 libraries have been built with the OpenGL option.

   From Release 9.6 the vis.mac macro in example B1 has additional commands
   that demonstrate additional functionality of the vis system, such as
   displaying text, axes, scales, date, logo and shows how to change
   viewpoint and style.  Consider copying these to other examples or
   your application.  To see even more commands use help or
   ls or browse the available UI commands in the Application
   Developers Guide, Section 7.1.

   For more information on visualization, including information on how to
   install and run DAWN, OpenGL and HepRApp, see the visualization tutorials,
   for example,
   http://geant4.slac.stanford.edu/Presentations/vis/G4[VIS]Tutorial/G4[VIS]Tutorial.html
   (where [VIS] can be replaced by DAWN, OpenGL and HepRApp)

   The tracks are automatically drawn at the end of each event, accumulated
   for all events and erased at the beginning of the next run.

 B- USER INTERFACES
 
   The user command interface is set via the G4UIExecutive class
   in the main() function in exampleB1.cc 
   The selection of the user command interface is then done automatically 
   according to the Geant4 configuration or it can be done explicitly via 
   the third argument of the G4UIExecutive constructor (see exampleB4a.cc). 
 
 C- HOW TO RUN

    - Execute exampleB1 in the 'interactive mode' with visualization:
        % ./exampleB1
      and type in the commands from run1.mac line by line:  
        Idle> /control/verbose 2
        Idle> /tracking/verbose 1
        Idle> /run/beamOn 10 
        Idle> ...
        Idle> exit
      or
        Idle> /control/execute run1.mac
        ....
        Idle> exit

    - Execute exampleB1  in the 'batch' mode from macro files 
      (without visualization)
        % ./exampleB1 run2.mac
        % ./exampleB1 exampleB1.in > exampleB1.out
        
  NB. Numbering scheme for histograms:
  layer     : from 1 to NbOfLayers (inclued)
  absorbers : from 1 to NbOfAbsor (inclued)
  planes    : from 1 to NbOfLayers*NbOfAbsor + 1 (inclued)
  
 One can control the binning of the histo with the command:
  /analysis/h1/set   idAbsor  nbin  Emin  Emax  unit 
  where unit is the desired energy unit for that histo (see TestEm3.in).
         
  One can control the name of the histograms file with the command:
  /analysis/setFileName  name  (default testem3)
   
  It is possible to choose the format of the histogram file : root (default),
  xml, csv, by using namespace in HistoManager.hh 
   	
 It is also possible to print selected histograms on an ascii file:
 /analysis/h1/setAscii id
 All selected histos will be written on a file name.ascii  (default testem3)
