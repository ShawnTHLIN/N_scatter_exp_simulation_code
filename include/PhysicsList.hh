#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "G4EmConfigurator.hh"
#include "globals.hh"

class G4VPhysicsConstructor;


class PhysicsList : public G4VUserPhysicsList
{
public:
	PhysicsList();
	virtual ~PhysicsList();

public:
	virtual void ConstructParticle();
	virtual void ConstructProcess();

	virtual void SetCuts();

	//these methods Construct physics processes and register them
	//void ConstructDecay();
	//void ConstructEM();
	void ConstructOp();

	//for the Messenger 
	void SetVerbose(G4int);
	void SetNbOfPhotonsCerenkov(G4int);

private:
	G4int    fVerboseLebel;
	G4VPhysicsConstructor* emPhysicsList;
	G4VPhysicsConstructor* decPhysicsList;
	G4VPhysicsConstructor* radDecayPhysicsList;

	G4VPhysicsConstructor* G4HadronElasticPhysicsHPList;
	G4VPhysicsConstructor* HadronPhysicsQGSP_BIC_HPList;
	G4VPhysicsConstructor* G4StoppingPhysicsList;
	G4VPhysicsConstructor* G4IonPhysicsList;
	G4VPhysicsConstructor*  G4EmExtraPhysicsList;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

