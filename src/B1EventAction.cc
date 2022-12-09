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
// $Id: B1EventAction.cc 75117 2013-10-28 09:38:37Z gcosmo $
//
/// \file B1EventAction.cc
/// \brief Implementation of the B1EventAction class

#include "B1EventAction.hh"
//#include "B1Run.hh"
#include "B1RunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

int Trigger_tag, recieve_tag_LO,recieve_tag_proton,Trigger_tag_for_LO;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::B1EventAction(B1RunAction* runAction)
: G4UserEventAction(), 
frunAction(runAction),
frawPhotons(0),
fdetPhotons(0),
fElectron(0)
//fcount(0)
  //fEdep(0.)



{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::~B1EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::BeginOfEventAction(const G4Event*)
{    
frawPhotons = 0;
fdetPhotons = 0;
fElectron = 0;
lightoutput =0;
lightoutput_temp=0;
proton_energy=0;
proton_energy_temp =0;
neutron_energy=0;
LO_in_n_detector=0;
//fcount = 0;
Trigger_tag =0;
recieve_tag_LO = 0;
recieve_tag_proton = 0;
Trigger_tag_for_LO = 0;
vect2 ={0,0,0};
vect3 ={0,0,0};
//frunAction->QMeVproton=false;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::EndOfEventAction(const G4Event*)
{   		//G4cout <<fNumofPhotons<< G4endl;
		//frunAction-> output(fNumofPhotons);
		//frunAction->outnumber(fNumofPhotons);

/*
if (InifEnprotons/CLHEP::MeV - 3 > 0.5){
	G4cout << "ENDOFEVENT"<<G4endl;
	G4cout << InifEnprotons << G4endl;
	G4cout << "ERR"<<G4endl;
	fgetc(stdin);
}
*/
/*
		if (fNumofPhotons != 0){
			frunAction-> output(fNumofPhotons);
			//G4cout << fNumofPhotons << G4endl;
		}
*/
	//G4cout << fNumofprotons << "protons are produced "<< G4endl;
	/*
	if (fcount == 1){
		frunAction->output(fElectron,0,frawPhotons);
	}
	if (fcount > 1){
		frunAction->output(0,fElectron,frawPhotons);
	}
*/

if ((Trigger_tag==1)&&(recieve_tag_LO==1)&&(recieve_tag_proton==1)&&(Trigger_tag_for_LO==1)){
	//G4cout<< "======"<< lightoutput<<"   "<< proton_energy<<"     "<< neutron_energy<<"   "<<G4endl;
	//fgetc(stdin);
	frunAction->output(lightoutput,proton_energy,neutron_energy,LO_in_n_detector);
	/*
	if ((lightoutput != 0) && (proton_energy != 0) && (neutron_energy != 0)&&(LO_in_n_detector != 0) ){
 		frunAction->fp << lightoutput<<" "<<proton_energy<<" "<<neutron_energy<<" "<<LO_in_n_detector<<" "<<" \n";
	}
	*/
}

	//frunAction->voutput(vect2);
	//frunAction->voutput_2(vect3);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::bpaddNum(G4double En)
{
	if (En>0){
		lightoutput = lightoutput+En;
	}
	if (En<0){
		recieve_tag_LO =1;
	}

/*
	if (En>0){
		recieve_tag_LO =1;
		lightoutput_temp = lightoutput_temp + En; // temporary save the accumulated lightoutput
		G4cout<< "====================recieve_tag_LO: "<< lightoutput_temp <<"================================  "<<G4endl;
		//fgetc(stdin);
	}
	if ((En<0)&&(recieve_tag_LO ==1)){
		
		lightoutput = 0+lightoutput_temp;   // if we get the signal that the neutron detector recieve signal, grather the saved lihgoutput.
		G4cout<< "==================== ===============================================  "<<G4endl;
		G4cout<< "==================== ===============================================  "<<G4endl;
		G4cout<< "==================== ===============================================  "<<G4endl;
		G4cout<< "===================== Eventaction lightoutput recorded:  "<< lightoutput<<"  =================  "<<G4endl;
		G4cout<< "==================== ===============================================  "<<G4endl;
		G4cout<< "==================== ===============================================  "<<G4endl;
		G4cout<< "==================== ===============================================  "<<G4endl;
	}

*/
}

void B1EventAction::cpaddNum(G4double En)
{
	if (En>0){
		proton_energy = proton_energy+En;
	}
	if (En<0){
		recieve_tag_proton =1;
	}

/*
	if (En>0){
		recieve_tag_proton =1;
		proton_energy_temp = proton_energy_temp + En; // temporary save the accumulated lightoutput
		G4cout<< "====================recieve_tag_proton================================  "<<G4endl;
		//fgetc(stdin);

	}
	if ((En<0)&&(recieve_tag_proton ==1)){
		
		proton_energy = 0+proton_energy_temp;   // if we get the signal that the neutron detector recieve signal, grather the saved lihgoutput.
		G4cout<< "=========================================================================================================================  "<<G4endl;
		G4cout<< "=========================================================================================================================  "<<G4endl;
		G4cout<< "==============================================================================================  "<<G4endl;
		G4cout<< "===================== Eventaction proton_energy recorded:  "<< proton_energy<<"  =================  "<<G4endl;
		G4cout<< "==============================================================================================  "<<G4endl;
		G4cout<< "==============================================================================================  "<<G4endl;
		G4cout<< "==============================================================================================  "<<G4endl;
	}
*/
}
void B1EventAction::dpaddNum(G4double En)
{

	if (En < 0){
		Trigger_tag = 1;
		//G4cout << " Trigger detector was triggered " << G4endl;
	}

	if (En > 0) {
		//fgetc(stdin);
		//G4cout << " Lightout in detector: "<< En << G4endl;

		//fdetPhotons =fdetPhotons +En;
		//G4cout << " Energy in detector "<< En  << G4endl;
		//fgetc(stdin);
		//if (En >2.5 &&En<3){
			//G4cout << " Energy on detector surface   "<< En  << G4endl;
			//fgetc(stdin);

		//}
		neutron_energy = 0+En ;

		
	}
	
//fdetPhotons =fdetPhotons +En;

}

void B1EventAction::epaddNum(G4double En)
{
	if (En>0){
		LO_in_n_detector = LO_in_n_detector+En;
	}
	if (En<0){
		Trigger_tag_for_LO =1;
	}
	
}



void B1EventAction::vp_event(std::vector<double> event_vector){
	vect2 = event_vector;

}

void B1EventAction::vp_event_2(std::vector<double> event_vector){
	vect3 = event_vector;

}

void B1EventAction::output(G4int En){


}




