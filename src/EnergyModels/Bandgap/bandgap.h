#ifndef _BANDGAP_H_
#define _BANDGAP_H_

#include "lattikem.h"

//numOpticsEnergies can have a major impact on performance for small systems with breakpoints=all or breakpoints=bias
//#define numOpticsEnergies 1000

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
	int numScalars;
	double fermiEnergy;
	double photocarrierEnergy;
	double interatomicEnergy;
	double voxelSize;
} OptoelectronicState;

enum bandgapStates {LIGHT_OFF, LIGHT_ON};
enum thermalDistributions {BOLTZMANN_DISTRIBUTION,FERMI_DISTRIBUTION};

void absorptionshoulder(double broadeningenergy, double Eb, double *shoulder, double *shoulderenergies, double *shoulder_llimit, double *shoulder_ulimit);
void OS_save(Configuration *config);
void OS_allocate(OptoelectronicState *OS);
void OS_free(OptoelectronicState *OS);
double OS_energy(crystal *crys,int *coord,Configuration *config);
double bg_energy(Configuration *config);
void bg_setup(double *newbandgapFunction(int numBandgapAlteringElements, int *numEachElement),char *newcoordElement,char *newbandgapAlteringElements,int newnumBandgapAlteringElements,int *coordinationShells,int newnumStates,
char *newrepulsiveElements,int newnumRepulsiveElements, double *newrepulsiveEnergies,int newrepulsiveShellSize, int newmaxRepulsiveCoordination,
void *traj_generator(Trajectory *traj));
void bg_registerSettings();
void bg_trajectory(Trajectory *traj);
void OS_combineWeighted(Configuration *configs, int numCombine, Configuration *outconfig, double *weights);
double interatomic_energy(Configuration *config);
#endif
