//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B1DetectorConstruction.cc 75117 2013-10-28 09:38:37Z gcosmo $
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
	G4NistManager* nistManager = G4NistManager::Instance();
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Envelope parameters
  //
  //G4double env_sizeXY = 20*cm, env_sizeZ = 30*cm;


//  G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //nistManager->FindOrBuildMaterial("G4_POLYETHYLENE");
  //nistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  //nistManager->FindOrBuildMaterial("G4_Al");
  nistManager->FindOrBuildMaterial("G4_Y");

  //G4Material* HDPE = nistManager->FindOrBuildMaterial("G4_POLYETHYLENE");
  //G4Material * PipeMat = nistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  //G4Material * PipeCovMat = nistManager->FindOrBuildMaterial("G4_Al");
  //G4Material * ColliMat = nistManager->FindOrBuildMaterial("G4_Al");
 // G4Material * AfterColliMat = nist->FindOrBuildMaterial("G4_Galactic");
	G4Material * TarMat = nist->FindOrBuildMaterial("G4_Y");
	G4Material * DetMat = nist->FindOrBuildMaterial("G4_Galactic");
	G4Material * TarAg =  nist->FindOrBuildMaterial("G4_Ag");
	//G4Material * al_mat =  nist->FindOrBuildMaterial("G4_Al");
  //     
	G4double z, a, density ,fractionmass;
	G4String name, symbol;
	G4int ncomponents, natoms;

	//a = 100*g / mole;
	//a = 67 * g / mole;
	//a = 203 * g / mole;
	//G4Element* elMo = new G4Element(name = "Mo", symbol = "Mo", z = 42 , a);
	//G4Element *elZn = new G4Element(name = "Zn", symbol = "Zn", z = 30, a);
	G4Element *elTl = new G4Element(name = "Tl", symbol = "Tl", z = 81, a= 203 * g / mole);

	//density =10.28 *g / cm3;
	//density = 7.14*g / cm3;
	//density = 11.85* g / cm3;
	
	//G4Material*Mo = new G4Material(name = "mo100", density, ncomponents = 1);
	//G4Material * Zn67 = new G4Material(name = "Zn67", density, ncomponents = 1);
	G4Material * Tl203 = new G4Material(name = "Tl203", density = 11.85* g / cm3, ncomponents = 1);
	

	//Mo->AddElement(elMo, natoms = 1);
	//Zn67->AddElement(elZn, natoms = 1);
	Tl203->AddElement(elTl, natoms = 1);

		
	///
		//stainless steel
	G4Element* Cr = nist->FindOrBuildElement("Cr");
	G4Element* Mn = nist->FindOrBuildElement("Mn");
	G4Element* Fe = nist->FindOrBuildElement("Fe");
	G4Element* Ni = nist->FindOrBuildElement("Ni");
	G4Element* Si = nist->FindOrBuildElement("Si");
	G4Element* C = nist->FindOrBuildElement("C");

	G4Material* StainlessSteel = new G4Material("StainlessSteel", density = 8.06*g / cm3, ncomponents = 6);
	StainlessSteel->AddElement(C, fractionmass = 0.001);
	StainlessSteel->AddElement(Si, fractionmass = 0.007);
	StainlessSteel->AddElement(Cr, fractionmass = 0.18);
	StainlessSteel->AddElement(Mn, fractionmass = 0.01);
	StainlessSteel->AddElement(Fe, fractionmass = 0.712);
	StainlessSteel->AddElement(Ni, fractionmass = 0.09);
	
	G4Material* NaI_mat = nist->FindOrBuildMaterial("G4_SODIUM_IODIDE");
	G4Material* ColliMat = nist->FindOrBuildMaterial("G4_Pb");

	// Lead Collimator material
	G4Element* Pb = nist->FindOrBuildElement("Pb");

	G4Material* Lead_Collimater = new G4Material("Lead_Collimater", density = 11.34*g / cm3, ncomponents = 1);
	Lead_Collimater->AddElement(Pb, fractionmass = 1);
	
	G4Element* Al = nist->FindOrBuildElement("Al");
	G4Element* H = nist->FindOrBuildElement("H");

	G4Material* al_mat = new G4Material("al_mat", density = 2.7*g / cm3, ncomponents = 1);
	al_mat->AddElement(Al, fractionmass = 1.0);
	
	
	G4Material* ej309_mat = new G4Material("ej309", density = 0.959*g / cm3, ncomponents = 2);
	ej309_mat->AddElement(H, natoms= 543);
	ej309_mat->AddElement(C, natoms= 435);
	
	G4Material* ej290_mat = new G4Material("ej290", density = 1.02*g / cm3, ncomponents = 2);
	ej290_mat->AddElement(H, natoms= 517);
	ej290_mat->AddElement(C, natoms= 465);	
	
	// World
  //
  //G4double world_sizeXY = 1.2*env_sizeXY;
  //G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
    
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       200*cm, 200*cm, 200*cm);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
               
	
	//
	//
	//
	//
	//
	//
	//											
	//										   |----------------------
	//										   |
	//										   |					  ^
	//       O-------------------------------> |					79.25 mm 
	//										   |					  
	//										   |
	//										   |----------------------
	//
	//                   150 cm								
	//
	//
	//
	//
	//
	//
	//
	//




	//********ej309*********//
	
/*	
  G4double distance_trigger =   30 *mm ;
	G4ThreeVector positionDet_trigger = G4ThreeVector(0, 0, distance_trigger);
	
	G4RotationMatrix* myRotation90 = new G4RotationMatrix();
	myRotation90->rotateX(0.*deg);
	myRotation90->rotateY(90.*deg);
	myRotation90->rotateZ(0.*rad);







	G4Tubs* Tri_DetCasing = 
		new G4Tubs("casing_Tri_Det", 0.0*cm, 52.8/2* mm , 52.8/2 *mm, 0.*deg, 360.*deg);
	G4Tubs * Tri_Det =
		new G4Tubs("Tri_Det", 0.0*cm, 50.8/2* mm , 50.8/2 *mm, 0 * deg, 360 * deg);

	G4LogicalVolume*logic_casing_Tri_Det =
		new G4LogicalVolume(
		Tri_DetCasing,            //its solid
		al_mat,                 //its material
		"casing_Tri_Det");

	G4LogicalVolume*logic_Tri_Det =
		new G4LogicalVolume(Tri_Det,            //its solid
		ej309_mat,             //its material
		"Tri_Det");              //its name

	G4VPhysicalVolume * Tri_DetCasing_phy =
		new G4PVPlacement(
		myRotation90,                       // rotation
		positionDet_trigger,     //at (0,0,0)
		logic_casing_Tri_Det,                //its logical volume
		"casing_Tri_Det",                   //its name
		logicWorld,             //its mother  volume
		false,                   //no boolean operation
		0,                       //copy number
		checkOverlaps);          //overlaps checking

	G4VPhysicalVolume * Tri_Det_phy =
		new G4PVPlacement(0,                       //no rotation
		G4ThreeVector(0, 0, 0),     //at (0,0,0)
		logic_Tri_Det,                //its logical volume
		"Tri_Det",                   //its name
		logic_casing_Tri_Det,             //its mother  volume
		false,                   //no boolean operation
		0,                       //copy number
		checkOverlaps);          //overlaps checking
		/////////////////////////
*/	

/*
	G4Box * Solid_Tri_Det =
		new G4Box("Tri_Det", 100/2*cm, 100/2* mm , 10/2 *mm);

	G4LogicalVolume*logic_Tri_Det =
		new G4LogicalVolume(Solid_Tri_Det,            //its solid
		ej290_mat,             //its material
		"Tri_Det");              //its name

	G4VPhysicalVolume * Tri_Det_phy =
		new G4PVPlacement(0,                       //no rotation
		positionDet_trigger,     //at (0,0,0)
		logic_Tri_Det,                //its logical volume
		"Tri_Det",                   //its name
		logicWorld,             //its mother  volume
		false,                   //no boolean operation
		0,                       //copy number
		checkOverlaps);          //overlaps checking
	*/

	G4double R90 = 90 * deg;
	G4RotationMatrix* myRotation90 = new G4RotationMatrix();
	myRotation90->rotateX(0.*deg);
	myRotation90->rotateY(R90);
	myRotation90->rotateZ(0.*rad);
	
	G4ThreeVector position_tri_det = G4ThreeVector(0,0, 0 );

	G4Tubs* tri_Det_Casing =
		new G4Tubs("casing_tri_det", 0.0*cm, 52.4/2* mm , 52.2/2 *mm, 0.*deg, 360.*deg);
	G4Tubs * tri_Det =
		new G4Tubs("tri_det", 0, 50.4/2* mm , 50.4/2 *mm, 0 * deg, 360 * deg);

	G4LogicalVolume*logic_tri_Det_Casing =
		new G4LogicalVolume(tri_Det_Casing,            //its solid
		al_mat,                 //its material
		"casing_tri_det");

	G4LogicalVolume*logic_tri =
		new G4LogicalVolume(tri_Det,            //its solid
		ej309_mat,             //its material
		"tri_det");              //its name
	
	G4VPhysicalVolume * tri_Det_Casing_phy =
		new G4PVPlacement(myRotation90,                       //no rotation
		position_tri_det,     //at (0,0,0)
		logic_tri_Det_Casing,                //its logical volume
		"casing_tri_det",                   //its name
		logicWorld,             //its mother  volume
		false,                   //no boolean operation
		0,                       //copy number
		checkOverlaps);          //overlaps checking
	G4VPhysicalVolume * tri_det_phy =
		new G4PVPlacement(0,                       //no rotation
		G4ThreeVector(0, 0, 0),     //at (0,0,0)
		logic_tri,                //its logical volume
		"tri_det",                   //its name
		logic_tri_Det_Casing,             //its mother  volume
		false,                   //no boolean operation
		0,                       //copy number
		checkOverlaps);          //overlaps checking

	G4double distance_main = 1470;
	G4double R30 = 45 * deg;  // angel of recording detector

	//G4ThreeVector positionDet30deg = G4ThreeVector(distance_main * std::sin(R30)*mm, 0, distance_main * std::cos(R30)*mm);
	G4ThreeVector positionDet30deg = G4ThreeVector(-distance_main * std::sin(R30)*mm, 0,-distance_main * std::cos(R30)*mm);

	G4RotationMatrix* myRotation30 = new G4RotationMatrix();
	myRotation30->rotateX(0.*deg);
	myRotation30->rotateY(-R30);
	myRotation30->rotateZ(0.*rad);

	G4Tubs* deg30DetCasing =
		new G4Tubs("casing@30deg", 0.0*cm, 78.2/2* mm , 78.2/2 *mm, 0.*deg, 360.*deg);
	G4Tubs * deg30Det =
		new G4Tubs("30deg", 0, 76.2/2* mm , 76.2/2 *mm, 0 * deg, 360 * deg);

	G4LogicalVolume*logic30DetCasing =
		new G4LogicalVolume(deg30DetCasing,            //its solid
		al_mat,                 //its material
		"casing@30deg");

	G4LogicalVolume*logic30deg =
		new G4LogicalVolume(deg30Det,            //its solid
		ej309_mat,             //its material
		"30deg");              //its name
	
	G4VPhysicalVolume * deg30DetCasing_phy =
		new G4PVPlacement(myRotation30,                       //no rotation
		positionDet30deg,     //at (0,0,0)
		logic30DetCasing,                //its logical volume
		"casing@30deg",                   //its name
		logicWorld,             //its mother  volume
		false,                   //no boolean operation
		0,                       //copy number
		checkOverlaps);          //overlaps checking
	G4VPhysicalVolume * deg30Det_phy =
		new G4PVPlacement(0,                       //no rotation
		G4ThreeVector(0, 0, 0),     //at (0,0,0)
		logic30deg,                //its logical volume
		"30deg",                   //its name
		logic30DetCasing,             //its mother  volume
		false,                   //no boolean operation
		0,                       //copy number
		checkOverlaps);          //overlaps checking
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
