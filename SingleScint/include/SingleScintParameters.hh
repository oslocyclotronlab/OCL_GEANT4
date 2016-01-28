/*
 * Parameters.hh
 *
 * Includes all the constants used in the creation of the
 * DetectorGeometry class
 * and in the source position in PrimaryActionGenerator class
 *
 *  Created on: Oct 22, 2015
 *      Author: fabiobz
 */

#ifndef SingleScintPARAMETERS__HH_
#define SingleScintPARAMETERS__HH_

#include "G4SystemOfUnits.hh"



    //
    // World
    //

	const G4double world_sizeXYZ = 80*cm;
    //
    // Detector & Shielding
    //

	const G4double crystalOuterR = 132.*mm/2.;
	const G4double crystalInnerR = 0.0*mm;
	const G4double crystalHalfLength = 89.*0.5*mm;
  	const G4double startPhi = 0.*deg;
  	const G4double deltaPhi = 360.*deg;

  	const G4double reflectorThickness = 0.5*mm; // assumption: 1 mm thick reflector on the front side
	const G4double reflectorHalfLength = crystalHalfLength + 0.5 * reflectorThickness; // assumption: backside doesn't have a reflector
	//const G4double ReflectorInnerR = crystalOuterR;
	const G4double reflectorInnerR = 0.*mm;
	const G4double reflectorOuterR = crystalOuterR + reflectorThickness;

  	const G4double coatingThickness = 0.5*mm; // thickness as in the radius part
  	const G4double coatingThicknessFront = 0.5*mm; // we assume a smaller thickness at the front of the detector
  	const G4double coatingOuterR = 134*mm/2. ;
	// in between reflector and coating, there will be some plastic
	  const G4double coatingPlasticThickness = coatingOuterR - coatingThickness -reflectorThickness - crystalOuterR; // assumption: 2.55 mm plexiglas coating around the reflector before the aluminium
  	const G4double coatingHalfLength = reflectorHalfLength + 0.5 * coatingThicknessFront + 0.5 * coatingPlasticThickness; // backside doesn't have an (Aluminium) coating

  	const G4double shieldingThickness = 5.*mm; 		// thickness of the tube
  	const G4double shieldingHalfThicknessLid = 2.*mm/2.;
  	const G4double shieldingInnerR = 0*mm; 			// as we use it as a mother volume
  	const G4double shieldingOuterR = coatingOuterR + shieldingThickness;

  	//in the front, the shielding tube diameter is reduces. It's later modeled by a conical section
  	const G4double shieldingConeHalfLength = 10.*mm;// in the front, the tube
  	//const G4double shieldingConeInnerRFront = coatingOuterR;
  	const G4double shieldingConeOuterRFront = coatingOuterR + 2.*mm;
  	//const G4double shieldingConeInnerRBack = shieldingConeInnerRFront;
  	const G4double shieldingConeOuterRBack = coatingOuterR + 5.*mm;

  	const G4double shieldingHalfLength = coatingHalfLength - shieldingConeHalfLength; // without conical Section and Lid
  																						 //  we assume no coating at the back side

  	const G4double plexiGlasWindowOuterR = shieldingOuterR; // currently we just assume a flat window on the top.
  	const G4double plexiGlasWindowHalfLength= 0.5 * 1.*mm;


    //
    // PMT
    //

	const G4double PMTWindowHalfLength = 1.0*mm;
	const G4double PMTWindowRadius = 85*0.5*mm;

	const G4double cathodeHalfLength = 0.005*mm;
	const G4double cathodeRadius =85*0.5*mm;


    //
	// Collimator and Source
    //

	// when you change the Collimator length and distance to Source, check that it's still inside the World Volume!

	const G4double distSourceCol = 12.*cm; 		// Distance from source to Collimator (beginning)
	const G4double collimatorHalfLength = 5.*cm; // adapt here for different collimator lengths

	// Distance from collimator End to Crystal Half point (or fraction r in crystal length)
	const G4double distHalfColHalfCry = collimatorHalfLength + 2.*shieldingHalfThicknessLid + coatingThicknessFront
									  + coatingPlasticThickness + reflectorThickness + crystalHalfLength;
	const G4double distSourceHalfCry =  distSourceCol + collimatorHalfLength + distHalfColHalfCry;


	// const G4double ratioInCrystal = 0.5;          // range: [0..1], defines point r from where the gammas can hit the crystal
	// const G4double distColEndPointToRatioCrsytal = distHalfColHalfCry - collimatorHalfLength + ( 2*ratioInCrystal - 1.) * crystalHalfLength;

	//parameters for collimator as a cone as a function of the parameters above
	// 1 is the front (towards source), 2 the backside
	// const G4double colRmin1 = crystalOuterR * (distSourceCol / (distSourceCol + 2.*collimatorHalfLength + distColEndPointToRatioCrsytal) );
	// const G4double colRmax1 = crystalOuterR * (distSourceCol / (distSourceCol + 2.*collimatorHalfLength) );
	// const G4double colRmin2 = crystalOuterR * (distSourceCol + 2*collimatorHalfLength)
	// 						 / ( distSourceCol + 2*collimatorHalfLength + distColEndPointToRatioCrsytal );
	// const G4double colRmax2 = shieldingConeOuterRFront;
  const G4double colRmin1 = 38.*mm/2.;
  const G4double colRmax1 = 72./2.*mm;
  const G4double colRmin2 = 70./2.*mm;
  const G4double colRmax2 = 135./2.*mm;



#endif /* PARAMETERS__HH_ */
