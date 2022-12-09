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
// $Id: B1EventAction.hh 75216 2013-10-29 16:08:11Z gcosmo $
//
/// \file B1EventAction.hh
/// \brief Definition of the B1EventAction class

#ifndef B1EventAction_h
#define B1EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"


class B1RunAction;
/// Event action class
///


class B1EventAction : public G4UserEventAction
{
  public:
	B1EventAction(B1RunAction* runAction);
    virtual ~B1EventAction();
    
    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

    //void AddEdep(G4double edep) { fEdep += edep; }


	//G4int frawPhotons,fdetPhotons;
	G4double fdetPhotons;
	G4double frawPhotons;
	G4double fElectron;
	G4int fcount;
	G4double lightoutput;
	G4double proton_energy;
	G4double neutron_energy;
	G4double proton_energy_temp;
	G4double lightoutput_temp;
	G4double LO_in_n_detector;
	//G4int Trigger_tag;

	void bpaddNum(G4double);
	void cpaddNum(G4double);
	void dpaddNum(G4double);
	void epaddNum(G4double);
	
	void output(G4int En);

	std::vector<double> event_vector;
	std::vector<double>  vect2;
	std::vector<double>  vect3;
	void vp_event(std::vector<double>);
	void vp_event_2(std::vector<double>);
  private:
    //G4double  fEdep;
      B1RunAction* frunAction;
	


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
