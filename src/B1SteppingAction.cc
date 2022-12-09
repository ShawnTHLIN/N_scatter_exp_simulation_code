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
// $Id: B1SteppingAction.cc 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file B1SteppingAction.cc
/// \brief Implementation of the B1SteppingAction class
#include <iostream>
#include <fstream>

#include "B1SteppingAction.hh"
#include "B1EventAction.hh"
#include "B1DetectorConstruction.hh"
#include "B1RunAction.hh"

#include "G4Threading.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4MPImanager.hh"
#include "G4SystemOfUnits.hh"
std::vector<double> Trigger_position_array;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::B1SteppingAction(B1EventAction* eventAction, B1RunAction *runAction)
: G4UserSteppingAction(),
fEventAction(eventAction),
frunAction(runAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::~B1SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1SteppingAction::UserSteppingAction(const G4Step* step)
{


#if 1
	if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume() != 0){
		if(step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ==("tri_det")){
			if ((step->GetTrack()->GetDefinition()->GetParticleName() == "proton")
			&&	(step->GetTrack()->GetCurrentStepNumber()==1 )){
				


				G4double En = step-> GetPreStepPoint()->GetKineticEnergy();
				//G4cout<< "proton energy" << En << G4endl;
				fEventAction-> cpaddNum(En);

				G4double Light_output = 0.62*En-1.3*(1-G4Exp(-1*0.39*(pow(En,0.97))))  ;
				//G4cout<< "Light output" << Light_output << G4endl;
	
				if (Light_output >= 0){
					fEventAction-> bpaddNum(Light_output);
					fEventAction-> dpaddNum(-1);
					fEventAction-> epaddNum(-1);
				}
			}
		}
	}
/*
	if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume() != 0){
		if(step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ==("tri_det")){
			fEventAction-> dpaddNum(-1);
		}
	}
*/



	if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume() != 0){
		if((step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ==("World"))
			&&(step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ==("casing@30deg"))){
			if (step->GetTrack()->GetDefinition()->GetParticleName() == "neutron"){

				//G4cout<<step->GetTrack()->GetCurrentStepNumber()<< G4endl;
				//G4double En = step-> GetPreStepPoint()->GetKineticEnergy();
				//fEventAction-> dpaddNum(En);
				//frunAction-> output(En,0,0);
				//G4cout<< En << G4endl;
				//if (En > 2.9 && En<3){
				//	fgetc(stdin);
				//}
				
				//if (step->GetTrack()->GetTrackID()==1 ){
					G4double En_1 = step-> GetPreStepPoint()->GetKineticEnergy();
					//frunAction-> output(0,En_1,0); 
					fEventAction-> dpaddNum(En_1);   //return the energy on the neutron detector surface
				//}

				//frunAction->output(0,0,En_dep);

				//fgetc(stdin);

			}
		}
	}
	if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume() != 0){
		if(step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ==("30deg")){
			if ((step->GetTrack()->GetDefinition()->GetParticleName() == "proton")
			&&	(step->GetTrack()->GetCurrentStepNumber()==1 )){

				G4double En_n = step-> GetPreStepPoint()->GetKineticEnergy();
				G4double Light_output_3inches = 0.817*En_n-2.63*(1-G4Exp(-1*0.297*(pow(En_n,1))))  ;
				//G4cout<< "Light output" << Light_output_3inches << G4endl;
	
				if (Light_output_3inches >= 0){
					fEventAction-> epaddNum(Light_output_3inches);
					fEventAction-> bpaddNum(-1);
					fEventAction-> cpaddNum(-1);
				}
			}
		}
	}
}

#else
//Check the neutron spectrum
if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume() != 0){
		if(step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ==("tri_det")){
			if ((step->GetTrack()->GetDefinition()->GetParticleName() == "neutron")&&
				(step->GetTrack()->GetCurrentStepNumber()==1 ))
				{
					G4double En_neutron = step-> GetPreStepPoint()->GetKineticEnergy();
					//G4cout<< "intial neutron energy:   "<<  En_neutron <<G4endl;
					frunAction -> output(En_neutron,0,0,0);
				}
		}
		}
}

#endif



	/*
	if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume() != 0){
		if(step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ==("Tri_Det")){
			if ((step->GetTrack()->GetDefinition()->GetParticleName() == "proton")
				&&	(step->GetTrack()->GetCurrentStepNumber()==1 )){
				
				
				G4double En_dep_Trigger = step-> GetPreStepPoint()->GetKineticEnergy();
				fEventAction-> bpaddNum(En_dep_Trigger);

				G4double LO_Tri_Det_index = (En_dep_Trigger)*1000; // intput data is from 0.01MeV with a step 0.001 MeV: convert to index
				G4double LO_Tri_Det  = LO_Birks_array[LO_Tri_Det_index]; 

				fEventAction-> cpaddNum(LO_Tri_Det);
				//G4cout<< "Position: " << Trigger_position[0]<<' '<<Trigger_position[1]<<' '<<Trigger_position[2] << G4endl;      
				//frunAction -> bp << fcount <<' '<<En_Track<< ' '<< Trigger_position[0]<<' '<<Trigger_position[1]<<' '<<Trigger_position[2] << G4endl;
				fEventAction-> dpaddNum(-1); // tag that the trigger detector recieved signals 
				
				
				G4ThreeVector Trigger_position = step-> GetPreStepPoint()->GetPosition();
				Trigger_position_array = {Trigger_position[0],Trigger_position[1],Trigger_position[2]};
				fEventAction-> vp_event(Trigger_position_array);

			}
		}
	}
	*/
// The signal detector recieved
/*
	if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume() != 0){
		if(step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ==("0deg")){
			if ((step->GetTrack()->GetDefinition()->GetParticleName() == "proton")
				&&	(step->GetTrack()->GetCurrentStepNumber()==1 )){

				//G4cout<<step->GetTrack()->GetCurrentStepNumber()<< G4endl;
				G4double En_dep = step-> GetPreStepPoint()->GetKineticEnergy();


				G4int LO_Birks_index = (En_dep)*1000; // intput data is from 0.01MeV with a step 0.001 MeV: convert to index
				G4double LO_Birks_single  = LO_Birks_array[LO_Birks_index];
				fEventAction-> dpaddNum(LO_Birks_single);
				
				G4ThreeVector Trigger_position = step-> GetPreStepPoint()->GetPosition();
				Trigger_position_array = {Trigger_position[0],Trigger_position[1],Trigger_position[2]};
				fEventAction-> vp_event_2(Trigger_position_array);

				//frunAction->output(0,0,En_dep);

				//fgetc(stdin);

			}
		}
	}
}
*/
	/*
if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume() != 0){
	if (step->GetTrack()->GetDefinition()->GetParticleName() == "neutron"){
		if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() == "Det0deg"){
			frunAction->output(0,step->GetPostStepPoint()->GetKineticEnergy(),0);
		}
		if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() == "Det90deg"){
			frunAction->output(0,0,step->GetPostStepPoint()->GetKineticEnergy());
		}
	}
}
*/
// tally photon insert detector
/*
if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume() != 0){
	if (step->GetTrack()->GetDefinition()->GetParticleName() == "opticalphoton"){
		if ((step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() == "Detector")
			&&(step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() == "VolEJ290")){
				
			fEventAction->dpaddNum();
		}
	}
}

	if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume() != 0){
		if(step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ==("NaI@90deg")){
			if ((step->GetTrack()->GetDefinition()->GetParticleName() == "e-")or(step->GetTrack()->GetDefinition()->GetParticleName() == "e+")){
				G4double En = step->GetTotalEnergyDeposit() ;
				//frunAction->output(En,0,0);
				fEventAction-> bpaddNum(En);
				//G4cout<< "e-&e+" << En << G4endl;
			}
		}
	}

if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume() != 0){
	if (step->GetTrack()->GetDefinition()->GetParticleName()=="gamma"){
		if((step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() !=("NaI@90deg"))&&
			(step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ==("NaI@90deg"))){
				frunAction->output(0,step->GetPostStepPoint()->GetKineticEnergy(),0);
			}
	}
}

	
	if(step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ==("NaI@0deg")){
		if (step->GetTrack()->GetDefinition()->GetParticleName() == "e-"){
			G4double En = step->GetTotalEnergyDeposit() ;
			
			//frunAction->output(En,0,0);
			fEventAction-> bpaddNum(En);
			//G4cout << En << G4endl;
			//fgetc(stdin);
		}
	}	
	
	if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume() != 0){
		if(step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ==("NaI@90deg")){
			if ((step->GetTrack()->GetDefinition()->GetParticleName() == "e-")or(step->GetTrack()->GetDefinition()->GetParticleName() == "e+")){
				G4double En = step->GetTotalEnergyDeposit() ;
				//G4double En = step->GetPostStepPoint()->GetKineticEnergy() ;
				//frunAction->output(En,0,0);
				fEventAction-> cpaddNum(En);
				//G4cout<< "e-&e+" << En << G4endl;
			}
		}
	}
	*/

	/*
if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume() != 0){
	if ((step->GetTrack()->GetDefinition()->GetParticleName() == "e-")){
		if(step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ==("NaI@0deg")){
		G4double En = step->GetTotalEnergyDeposit() ;
		//frunAction->output(En,0,0);
		fEventAction-> cpaddNum(En);
		//G4cout<<"e-in NaI" << En << G4endl;
		}

	}
}
*/
	/*
	if(step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ==("NaI@0deg")){
		if (step->GetTrack()->GetDefinition()->GetParticleName() == "gamma"){
			G4double Eg = step->GetTotalEnergyDeposit() ;
			
			//frunAction->output(0,Eg,0);
			//fEventAction-> dpaddNum(Eg);
			//G4cout << En << G4endl;
			//fgetc(stdin);
		}
	}	
*/
/*
if (step->GetTrack()->GetDefinition()->GetParticleName() == "e-"){
	G4double En = step->GetTotalEnergyDeposit() ;
	//frunAction->output(En,0,0);
	fEventAction-> cpaddNum(En);
	//G4cout << En << G4endl;
}
*/
/*
if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume() != 0){
	if ((step->GetTrack()->GetDefinition()->GetParticleName()=="gamma")&&(step->GetTrack()->GetTrackID()==1)){
		if((step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() !=("NaI@90deg"))&&
			(step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ==("NaI@90deg"))){
				fEventAction -> dpaddNum(1);
			}
	}
}
*/
				
				
				/*
if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume() != 0){
	if ((step->GetTrack()->GetDefinition()->GetParticleName()=="gamma")&&(step->GetTrack()->GetTrackID()==1)){
		if ((step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="compt")||
			(step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="phot")){
			G4double Ei = step->GetPreStepPoint()->GetKineticEnergy();
			G4double Ef = step->GetPostStepPoint()->GetKineticEnergy();
			fEventAction -> dpaddNum(Ei-Ef);
			//G4cout << Ei-Ef << G4endl;
			//fgetc(stdin);
		}
	}
}


	if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume() != 0){
		if((step->GetTrack()->GetDefinition()->GetParticleName()=="gamma")&&(step->GetTrack()->GetTrackID()!=1)){
			if(step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ==("NaI@0deg")){
				if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="phot"){
					G4double Ei2 = step->GetPreStepPoint()->GetKineticEnergy();
					G4double Ef2 = step->GetPostStepPoint()->GetKineticEnergy();
					fEventAction-> dpaddNum(Ei2-Ef2);
					//G4cout << Ei2-Ef2 << G4endl;
					//fgetc(stdin);
					}
			else {
				fEventAction -> dpaddNum(0);
				}
			}
		}
	}
*/

	/*
	if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="phot"){
		G4double Ei = step->GetPreStepPoint()->GetKineticEnergy();
		G4double Ef = step->GetPostStepPoint()->GetKineticEnergy();
		fEventAction-> bpaddNum(Ei-Ef);
	}
	*/



	/*
	if (step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() == "Volume1"){
		if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() == "Volume2"){
			
			G4ThreeVector InitialDir = step->GetPreStepPoint()->GetMomentumDirection();
			G4double Inix, Iniy, Iniz, IniTheta;
			Inix = InitialDir.getX();
			Iniy = InitialDir.getY();
			Iniz = InitialDir.getZ();
			IniTheta = InitialDir.getTheta()/CLHEP::degree;

			G4ThreeVector FinalDir = step->GetPostStepPoint()->GetMomentumDirection();
			G4double Finalx, Finaly, Finalz, FinalTheta;
			Finalx = FinalDir.getX();
			Finaly = FinalDir.getY();
			Finalz = FinalDir.getZ();
			FinalTheta = FinalDir.getTheta()/CLHEP::degree;

			if (Finalz > 0){
				frunAction->fp << Inix << " " << Iniy << " " << Iniz << " " << IniTheta << " " << sin(IniTheta) << G4endl;
				frunAction->fp << Finalx << " " << Finaly << " " << Finalz << " " << FinalTheta << " " << sin(FinalTheta) << G4endl<<G4endl;
				frunAction->Refraction++;
			}
			if (Finalz < 0){
				frunAction->Reflection++;
			}
		}
	}
	
	
	if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume() != 0){
		if (step->GetTrack()->GetDefinition()->GetParticleName() == "gamma"){
			if (step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() == "Volume1"){
				if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() == "World"){
					//system("pause");
					step->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
					frunAction->NUMofSkip++;
				}
			}
		}
		if (step->GetTrack()->GetDefinition()->GetParticleName() == "e-"){
			if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() == "Volume1"){
				step->GetTrack()->SetTrackStatus(fStopButAlive) ;
				system("pause");
			}
		}
		if (step->GetTrack()->GetDefinition()->GetParticleName() == "opticalphoton"){
			frunAction->fp << step->GetTrack()->GetKineticEnergy() / CLHEP::eV << G4endl;

			frunAction->NUM++;
			//G4cout << frunAction->NUM << G4endl;

			//G4cout << "STEPPING" << G4endl;
			//system("pause");
		}
	}
	*/
/*
	if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume() != 0){
		//G4cout << frunAction -> QMeVproton <<"1"<< G4endl;

		
		
			if (step->GetTrack()->GetDefinition()->GetParticleName() == "opticalphoton"){
				if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() == "VolEJ290"){
					if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() == "initStep"){
						G4cout << G4endl << G4endl;
						//G4cout << fEventAction->NumofPhotons<< G4endl;
						//G4cout << frunAction -> QMeVproton<<"  3" << G4endl;
						fEventAction->addNum();
					}
				}
			}
		
	}
*/
/*
	if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume() != 0){
		//G4cout << frunAction -> QMeVproton <<"1"<< G4endl;
		if (step->GetTrack()->GetDefinition()->GetParticleName() == "proton"){
			G4double ProtonEnergy;
				ProtonEnergy = step->GetTrack()->GetKineticEnergy()/CLHEP::MeV;
				
				step->GetTrack()->SetTrackStatus(fStopAndKill);
				
				G4cout <<"STEPPING "<< ProtonEnergy << G4endl;
				fEventAction-> fEnprotons=fEventAction-> fEnprotons+ProtonEnergy;
				G4cout << "ACCUMLATING "<<fEventAction-> fEnprotons<<G4endl;

				

			//if (ProtonEnergy>0*MeV){
			//	frunAction-> output(ProtonEnergy);
			//}
		}
	}	

*/
/*
	if (step->GetTrack()->GetDefinition()->GetParticleName() == "neutron"){
		G4String ProcessName;
		G4String Trans="Transportation";
		G4String Elas="hadElastic";
		
		ProcessName =  step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();

		G4cout << G4endl << step->GetPostStepPoint()->GetKineticEnergy() << "MeV"<<G4endl;
		G4cout << ProcessName << G4endl;
		//if ((ProcessName!= Trans)&&(ProcessName!=Elas)){
		//	G4cout <<123 <<G4endl;
		//}
	}
*/


	



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

