#include "PhysicsList.hh"

#include "G4IonPhysics.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4ProcessManager.hh"

#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4LossTableManager.hh"
#include "G4EmSaturation.hh"

#include "G4EmStandardPhysics_option4.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ios.hh"
#include "G4StepLimiter.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList()
: G4VUserPhysicsList(),
fVerboseLebel(0)  //Tan: choose the level of details about each interaciton displayed on the screen
{
	// default values.
	//G4Cerenkov::SetMaxNumPhotonsPerStep(20);
	//G4Cerenkov::SetMaxBetaChangePerStep(10.0);
	//G4Cerenkov::SetTrackSecondariesFirst(true);
	//G4Scintillation::SetScintillationYieldFactor(1.);
	//G4Scintillation::SetTrackSecondariesFirst(true);
	//G4Scintillation::SetScintillationByParticleType(true);


	decPhysicsList = new G4DecayPhysics();
	radDecayPhysicsList = new G4RadioactiveDecayPhysics();
	emPhysicsList = new G4EmStandardPhysics_option4(1);

	G4HadronElasticPhysicsHPList = new G4HadronElasticPhysicsHP(1);
	HadronPhysicsQGSP_BIC_HPList = new G4HadronPhysicsQGSP_BIC_HP(1);
	G4StoppingPhysicsList = new G4StoppingPhysics();
	G4IonPhysicsList = new G4IonPhysics();
	G4EmExtraPhysicsList = new G4EmExtraPhysics();

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList() {
	delete decPhysicsList;
	delete radDecayPhysicsList;
	delete emPhysicsList;
	delete G4HadronElasticPhysicsHPList;
	delete HadronPhysicsQGSP_BIC_HPList;
	delete G4StoppingPhysicsList;
	delete G4IonPhysicsList;
	delete G4EmExtraPhysicsList;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
	// In this method, static member functions should be called
	// for all particles which you want to use.
	// This ensures that objects of these particle types will be
	// created in the program.
	decPhysicsList->ConstructParticle();
	//G4BosonConstructor bConstructor;
	// bConstructor.ConstructParticle();

	//G4LeptonConstructor lConstructor;
	// lConstructor.ConstructParticle();

	//G4MesonConstructor mConstructor;
	//mConstructor.ConstructParticle();

	//G4BaryonConstructor rConstructor;
	//rConstructor.ConstructParticle();

	//G4IonConstructor iConstructor;
	//iConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
	AddTransportation();
	//ConstructDecay();
	//ConstructEM();
	emPhysicsList->ConstructProcess();
	decPhysicsList->ConstructProcess();
	radDecayPhysicsList->ConstructProcess();
	G4HadronElasticPhysicsHPList->ConstructProcess();
	HadronPhysicsQGSP_BIC_HPList->ConstructProcess();
	G4StoppingPhysicsList->ConstructProcess();
	G4IonPhysicsList->ConstructProcess();
	G4EmExtraPhysicsList->ConstructProcess();
	ConstructOp();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Decay.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4Threading.hh"

void PhysicsList::ConstructOp()
{
	//G4Cerenkov* cerenkovProcess = new G4Cerenkov("Cerenkov");
	//cerenkovProcess->SetMaxNumPhotonsPerStep(20);
	//cerenkovProcess->SetMaxBetaChangePerStep(10.0);
	//cerenkovProcess->SetTrackSecondariesFirst(true);
	G4Scintillation* scintillationProcess = new G4Scintillation("Scintillation");
	scintillationProcess->SetScintillationYieldFactor(1.);
	scintillationProcess->SetTrackSecondariesFirst(false);
	scintillationProcess->SetScintillationByParticleType(true);

	G4OpAbsorption* absorptionProcess = new G4OpAbsorption();
	G4OpRayleigh* rayleighScatteringProcess = new G4OpRayleigh();
	G4OpMieHG* mieHGScatteringProcess = new G4OpMieHG();
	G4OpBoundaryProcess* boundaryProcess = new G4OpBoundaryProcess();

	//cerenkovProcess->SetVerboseLevel(fVerboseLebel);
	scintillationProcess->SetVerboseLevel(fVerboseLebel);
	absorptionProcess->SetVerboseLevel(fVerboseLebel);
	rayleighScatteringProcess->SetVerboseLevel(fVerboseLebel);
	mieHGScatteringProcess->SetVerboseLevel(fVerboseLebel);
	boundaryProcess->SetVerboseLevel(0); //Tan: choose the level of details about each interaciton displayed on the screen

	// Use Birks Correction in the Scintillation process
	if (!G4Threading::IsWorkerThread())
	{
		G4EmSaturation* emSaturation =
			G4LossTableManager::Instance()->EmSaturation();
		//G4Scintillation::AddSaturation(emSaturation);
		scintillationProcess->AddSaturation(emSaturation);
	}
	// G4cout << "************************************************************************************" << G4endl;
	
	auto theParticleIterator=GetParticleIterator();

  	theParticleIterator->reset();
 	
	 while( (*theParticleIterator)() ){
		G4ParticleDefinition* particle = theParticleIterator->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();
		G4String particleName = particle->GetParticleName();
		//if (cerenkovProcess->IsApplicable(*particle)) {
		//	pmanager->AddProcess(cerenkovProcess);
		//	pmanager->SetProcessOrdering(cerenkovProcess, idxPostStep);
		//}
		if (scintillationProcess->IsApplicable(*particle)) {
			pmanager->AddProcess(scintillationProcess);
			pmanager->SetProcessOrderingToLast(scintillationProcess, idxAtRest);
			pmanager->SetProcessOrderingToLast(scintillationProcess, idxPostStep);
		}
		if (particleName == "opticalphoton") {
			G4cout << " AddDiscreteProcess to OpticalPhoton " << G4endl; 
//Tan: This portion runs at the beginning of the simulation
			pmanager->AddDiscreteProcess(absorptionProcess);
			pmanager->AddDiscreteProcess(rayleighScatteringProcess);
			pmanager->AddDiscreteProcess(mieHGScatteringProcess);
			pmanager->AddDiscreteProcess(boundaryProcess);
		}
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetVerbose(G4int verbose)
{
	fVerboseLebel = verbose;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetNbOfPhotonsCerenkov(G4int MaxNumber)
{
	//cerenkovProcess->SetMaxNumPhotonsPerStep(MaxNumber);
	//fMaxNumPhotonStep = MaxNumber;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts()
{
	//  " G4VUserPhysicsList::SetCutsWithDefault" method sets
	//   the default cut value for all particle types
	//
	SetCutsWithDefault();

	if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
