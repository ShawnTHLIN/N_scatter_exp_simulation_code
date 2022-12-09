/*
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
// $Id: B1PrimaryGeneratorAction.cc 69565 2013-05-08 12:35:31Z gcosmo $
//
/// \file B1PrimaryGeneratorAction.cc
/// \brief Implementation of the B1PrimaryGeneratorAction class

	fGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Plane");
	fGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Circle");
	fGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(0 * mm, 0 * mm, -15 * cm));
  fGPS->GetCurrentSource()->GetPosDist()->SetRadius(1.0*cm);

	//fParticleGun->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(0., 0.,0));
	//fParticleGun->GetCurrentSource()->GetPosDist()->SetBeamSigmaInX(0.667*cm);
	//fParticleGun->GetCurrentSource()->GetPosDist()->SetBeamSigmaInY(0.667*cm);

	fGPS->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1));
	//if (index < total - 1){
	//  fParticleGun->AddaSource(1);
	// }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
  delete fGPS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
*/
#include "B1PrimaryGeneratorAction.hh"
#include <iostream>
using namespace std;

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"	
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"


G4double rnd_number,possible_x;
G4int length_of_xaxis;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
G4String particleName;

B1PrimaryGeneratorAction::B1PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0),
  fGPS(0),fgamma(0)
  //fEnvelopeBox(0)
{
	//G4int n_particle = 9;///////1->3
  fGPS = new G4GeneralParticleSource();
  fgamma = particleTable->FindParticle("neutron");
  //fGPS->SetNumberOfParticles();
  //fGPS->SetParticlePosition(G4ThreeVector(0,0, -2.54 * cm));

	  fGPS->GetCurrentSource()->SetParticleDefinition(fgamma);
	  //fGPS->SetCurrentSourceIntensity(weighting[index]);
	  //fGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Gauss");
	  //....oooOO0OOooo........////....oooOO0OOooo........////....oooOO0OOooo........////....oooOO0OOooo........//
	  //fGPS->GetCurrentSource()->GetEneDist()->SetMonoEnergy(10.0 * MeV);///////
	  //....oooOO0OOooo........////....oooOO0OOooo........////....oooOO0OOooo........////....oooOO0OOooo........//
	  //fGPS->GetCurrentSource()->GetEneDist()->SetBeamSigmaInE(stdd[groupindex[index]]*MeV);
	  fGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(0 * mm, 0 * mm, 3 * cm));//3cm
	  fGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Beam");
	  fGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Circle");
	  fGPS->GetCurrentSource()->GetPosDist()->SetRadius(0*cm);

	  //fGPS->GetCurrentSource()->GetEneDist()->SetBeamSigmaInE(stdd[groupindex[index]]*MeV);

	  //fGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Beam");
    
	  //fGPS->GetCurrentSource()->GetPosDist()->SetBeamSigmaInX(xstdd[groupindex[index]]*mm);
	  //fGPS->GetCurrentSource()->GetPosDist()->SetBeamSigmaInY(ystdd[groupindex[index]]*mm);

	  //fGPS->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(0,0,-1));  
   	//fGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Plane");
	//fGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Circle");
	//fGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(0 * mm, 0 * mm, -15 * cm));
    //fGPS->GetCurrentSource()->GetPosDist()->SetRadius(1.0*cm);
    

 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
  //delete fParticleGun;
  delete fGPS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of ecah event
  //

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get Envelope volume
  // from G4LogicalVolumeStore.
    // 0.4367 Y //4657


	//if (index < total - 1){
	//fGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Plane");
	//fGPS->GetCurrentSource()->GetPosDist()->ConfineSourceToVolume("True");
	//fGPS->GetCurrentSource()->GetPosDist()->SetBeamSigmaInR(0.2 * cm);

	//fGPS->GetCurrentSource()->GetPosDist()->SetBeamSigmaInX(0.2 * cm);
	//fGPS->GetCurrentSource()->GetPosDist()->SetBeamSigmaInY(0.2 * cm);

	//no correct but approach
	//G4double theta = 69*G4UniformRand()*deg; //tan^-1(2.54/5.54) = 24.6 degree
	//G4double phi   = 360 *G4UniformRand()*deg;
	//G4cout<< "theta:  "<< theta << "phi: "<< phi<<G4endl;
	//fgetc(stdin);

	//G4ThreeVector dir(std::sin(theta)*std::cos(phi),std::sin(theta)*std::sin(phi),-1*abs(std::cos(theta)));
	//G4ThreeVector dir(0,0,-1);

	//****** This part is for generate AmBe neutron spectrum *********////
	//G4cout <<" x_axis size at primary generator " <<x_axis.size()  << G4endl;
	rnd_number = G4UniformRand(); // random seed for produce spectrum 
	length_of_xaxis=x_axis.size(); // for checking the number
	//G4cout <<" rnd_number " << rnd_number<< G4endl;

	for (int i4=0; i4<length_of_xaxis; i4++){
		if (rnd_number>ISO_cdf[i4] && rnd_number<ISO_cdf[i4+1]){
			//G4cout <<" i4 " << i4 << G4endl;
    		possible_x = x_axis[i4]+ ((rnd_number-ISO_cdf[i4])/(ISO_cdf[i4+1]-ISO_cdf[i4]))*(x_axis[i4+1]-x_axis[i4]);// interpolation
			break;
		}
	}
	fGPS->GetCurrentSource()->GetEneDist()->SetMonoEnergy(possible_x);
	//G4cout <<"AmBe emitter neutron Energy " << possible_x<< G4endl;
	//

	G4double theta = pi*G4UniformRand(); // maximum theta = 38.1/1500 = 0.0254
	G4double psi   = twopi*G4UniformRand();  //psi uniform in [0, 2*pi] 
	//G4ThreeVector dir(std::cos(psi)*std::tan(theta),std::sin(psi)*std::tan(theta),1);
	//G4ThreeVector dir(abs(sinAlpha*std::cos(psi)), sinAlpha*std::sin(psi), cosAlpha);
	
	G4ThreeVector dir(std::sin(theta)*std::cos(psi),std::sin(theta)*std::sin(psi),-1*abs(std::cos(theta)));// force to z minus
	
	fGPS->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(dir);

	fGPS->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

