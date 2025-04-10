#ifndef _BANDGAP_H_
#define _BANDGAP_H_

#include "lattikem.h"
#include "sourcecontrol.h"



#define optics_llimit 1.0
#define optics_ulimit 2.5
typedef struct OptoelectronicState{
	
    double *energystates;
    double *weights;    
    double *chargeDensity;    
    double *emission;
	double *absorption;
    double *DOS;
    double *current;
	double *photonenergies;
	//These are pointers which are used for compact Data copying and saving
	double **scalars;
	
	int *coord;
	//Additional network which has been collapsed.
	crystalnetwork *bandgapNetwork;
	
	int numScalars;
	double *fermiEnergy;
	double *photocarrierEnergy;
	double interatomicEnergy;
	double voxelSize;
#ifdef __CUDACC__	
	double *d_energyLevels, *d_weights, *d_energyWeights, *d_weightSum, *d_fermiEnergy, *d_energy;
#endif
} OptoelectronicState;

//enum bandgapStates {LIGHT_OFF, LIGHT_ON};
//enum bandgapGPU {GPU_OFF, GPU_ON};
enum thermalDistributions {BOLTZMANN_DISTRIBUTION,FERMI_DISTRIBUTION};
enum pulsedModes { CONSTANT_WAVE, PULSED, CROSSOVER};

void absorptionshoulder(double broadeningenergy, double Eb, double *shoulder, double *shoulderenergies, double *shoulder_llimit, double *shoulder_ulimit);
void OS_save(Configuration *config);
void OS_allocate(Configuration *config);
void OS_free(Configuration *config);
void OS_energy(crystal *crys,int *coord,Configuration *config);
double bg_energy(Configuration *config);

void bg_setup(double newbandgapFunction(int numBandgapAlteringElements, int *numEachElement),char *newcoordElement,char *newbandgapAlteringElements,int newnumBandgapAlteringElements,int *coordinationShells,int newnumStates,
char *newrepulsiveElements,int newnumRepulsiveElements, double *newrepulsiveEnergies,int newrepulsiveShellSize, int newmaxRepulsiveCoordination,
void traj_generator(Trajectory *traj));

void bg_registerSettings();
void bg_trajectory(Trajectory *traj);
void OS_combineWeighted(Configuration *configs, int numCombine, Configuration *outconfig, double *weights);
void OS_printCoord(OptoelectronicState *OS);
double interatomic_energy(Configuration *config);
void bg_networkSwap(Configuration *config, int atom1, int atom2);
void bg_coordSwap(Configuration *config, int atom1, int atom2, int forward_reverse);
double getFermiConvergence();
#endif
