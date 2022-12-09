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
// $Id: B1RunAction.cc 87359 2014-12-01 16:04:27Z gcosmo $
//
/// \file B1RunAction.cc
/// \brief Implementation of the B1RunAction class
#include <iostream>
#include <fstream>
#include <string>
#include "B1RunAction.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "B1DetectorConstruction.hh"
//#include "B1Run.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include "G4MPImanager.hh"
#include <G4Run.hh>

#include <sstream> //for std::stringstream 
#include <ctime>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<double> ISO_cdf,x_axis;
int count =0;
int test_counts=0;
G4int runID;
std::vector< std::vector<double> > vec1;
std::vector<double> tempVec;
G4int length_Of_ISO_cdf;
//std::vector<std::vector<double>> vect;
B1RunAction::B1RunAction()
: G4UserRunAction()
{ 
  // add new units for dose
  // 
	G4RunManager::GetRunManager()->SetPrintProgress(1000);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::~B1RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::BeginOfRunAction(const G4Run* )
{ 
  G4int rank_file = G4MPImanager::GetManager()-> GetRank();
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  //fp.open("Particle.txt", std::ios_base::app); 
  G4String TargetEn ;
  //TargetEn = "20220921_EJ309_N_scatter_3inches";//////
  TargetEn = "Output_exp_simulation_AmBe_source_45deg_2E9/202210_EJ309_N_scatter_exp_45deg_2E9_";
  //TargetEn = "test_energy";
  G4String Applied_En ;
  //....oooOO0OOooo........//
  //Applied_En = "10.0MeV";
  Applied_En = "AmBe";
  
  G4String numberOfrank = std::to_string(rank_file);
  //....oooOO0OOooo........//
  //TargetEn = "EJ309_TOF_trS_response_3.0MeV";
  G4String 	//bpFile = TargetEn + "b_LO_in_trigger_" + Applied_En + ".txt",
	          //cpFile = TargetEn + "c_proton_in_trigger_" + Applied_En + ".txt",
	          //dpFile = TargetEn + "d_E_neutron_detector_" + Applied_En + ".txt",
            //epFile = TargetEn + "e_LO_in_n_detector_" + Applied_En + ".txt",
            fpFile = TargetEn + "All_in_one_" + Applied_En +"_" + numberOfrank+ ".txt";
            //gpFile = "Output_exp_simulation/" + TargetEn + "All_in_one_postwrite_" + Applied_En +".txt";
  //bp.open(bpFile , std::ios_base::app);
  //cp.open(cpFile , std::ios_base::app);
  //dp.open(dpFile , std::ios_base::app);
  //ep.open(epFile , std::ios_base::app);
  fp.open(fpFile , std::ios_base::app);
  //gp.open(gpFile , std::ios_base::app);
  G4int read_count=0;
  ISO_cdf={};

    std::ifstream ifs1("ISO_spectrum_cdf.txt", std::ios::in);
    if (!ifs1.is_open()) {
          G4cout << "Failed to open file.\n"<<G4endl;
      } 
	  if (ifs1.is_open()) {
          float temp;
          while (ifs1 >> temp) {
            ISO_cdf.push_back(temp);
            G4cout <<" temp " << temp<< G4endl;
            read_count++;
         }
          ifs1.close();
          G4cout <<" read_count  at runaction" << read_count << G4endl;
          //fgetc(stdin);
      }
    G4cout <<" ISO_cdf  at runaction" << ISO_cdf.size() << G4endl;
    
    x_axis={};
	  length_Of_ISO_cdf = ISO_cdf.size();
	  for (int i=0;i<length_Of_ISO_cdf;i++){
    	x_axis.push_back(0.01*i);
		}
		G4cout <<" x_axis size at runaction" <<x_axis.size()  << G4endl;
		//fgetc(stdin);

    //ifs1.close();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::EndOfRunAction(const G4Run* )
{
	//fp << "There are " << Refraction << "s of photons are refracted" << G4endl;
	//bp << "There are " << Reflection << "s of photons are reflected" << G4endl;
	//bp << NUMofSkip << " gammas are skipped" << G4endl;

//fp.close();
/*
if (test_counts==0){
  for (int ii=0;ii < vec1.size();ii++){ 
    for (int iii=0;iii < vec1[0].size();iii++){
      gp << vec1[ii][iii]<<" ";
      //G4cout << vec1[ii][iii] << "  " ;
    }
    gp<< "\n";
    //G4cout  << G4endl ;
  }
}

gp << " ============== "<< "\n";
gp << " test_counts "<< test_counts<<"\n";
G4cout <<" test_counts "<< test_counts << G4endl ;

test_counts++;

gp.close();
*/
//bp.close();
//cp.close();
//dp.close();
//ep.close();
fp.close();

}

void B1RunAction::output(G4double bptotal,G4double cptotal,G4double dptotal,G4double eptotal){
//G4int thr = G4Threading::G4GetThreadId();
//G4int cor = G4Threading::G4GetNumberOfCores ();
//G4int runID = G4MPImanager::GetManager()-> GetSize();
G4int rank = G4MPImanager::GetManager()-> GetRank();
count++;
//if (QMeVproton==true){
//fp << total << "\n";
//Result = initial - Final ;
//bp << initial << " " << Final << " " << Result << "\n";
/*
if (bptotal != 0){
	bp  << " "<< count<<" "<<bptotal << "\n" ;
  //bp  <<bptotal << "\n" ;
  //bp << 1 << "\n" ;
}

if (cptotal != 0){
  cp << " "<< count<<" "<< cptotal << "\n" ;
  //cp << 2 << "\n" ;
}

if (dptotal != 0){
	dp << " "<< count<<" "<< dptotal << "\n" ;
  //dp << 3 << "\n" ;
}

if (eptotal != 0){
	ep  << " "<< count<<" "<< eptotal << "\n" ;
  //ep << 4 << "\n" ;
  /*
  if (dptotal == 0){
	  G4cout << dptotal << G4endl ;
    fgetc(stdin);
  } 
  */
//}

if ((bptotal != 0) && (cptotal != 0) && (dptotal != 0)&&(eptotal != 0) ){
  fp <<rank <<" "<<count<<" "<< bptotal<<" "<<cptotal<<" "<<dptotal<<" "<<eptotal<<" "<<" \n";
  //fp<<" "<< rank << " "<< count<<" "<< 1 <<" "<<2<<" "<<3<<" "<<4<<"\n";
/*
  tempVec={};
  tempVec.push_back(rank);
  tempVec.push_back(count);
  tempVec.push_back(bptotal);
  tempVec.push_back(cptotal);
  tempVec.push_back(dptotal);
  tempVec.push_back(eptotal);

  G4cout <<"tempVec:  "<< tempVec[0]<<"  " <<tempVec[1]<<"  "<<tempVec[2]<<"  "<<tempVec[3]<<"  "<<tempVec[4]<<"  "<<tempVec[5]<<"  "<< G4endl ;
  
  vec1.insert(vec1.end(), tempVec);
  */
 /*
  G4cout <<"vec1 size:  "<< vec1.size()<< G4endl ;
  G4cout <<"vec1[0] size:  "<< vec1[0].size()<< G4endl ;
  G4cout <<"vec1   "<<" ";
  for (int i1=0;i1<6;i1++){
     G4cout << vec1[vec1.size()-1][i1] << " " ;
  }
  */
}





//bp << "this" << "\n";
//G4cout <<"I am on "<< this <<G4endl;
//}
}
void B1RunAction::voutput(std::vector<double> vector_p){

//ep << vector_p[0]<<" "<<vector_p[1]<<" "<<vector_p[2] << "\n" ;

}

void B1RunAction::voutput_2(std::vector<double> vector_p_2){

//fp << vector_p_2[0]<<" "<<vector_p_2[1]<<" "<<vector_p_2[2] << "\n" ;

}


//**************************//

/*
void B1RunAction::Poutput(G4String Particle){
  fp << Particle << "\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
*/