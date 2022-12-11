#ifndef _ELECTROSTATICS_H_
#define _ELECTROSTATICS_H_

#include "sourcecontrol.h"

#define coulombconstant 14.40 //eV * angstrom / (electron charge)^2
#define permittivity 0.005524 //electron charge / angstrom / volt
//coulomb constant = 1/(4*pi*permittivity)

enum electrostaticExternalStates {NO_BIAS, ION_ONLY, ION_FORWARD_FIELD, ION_REVERSE_FIELD};
enum voltageMode {VOLTAGE_SQUARE_WAVE, VOLTAGE_SINE_WAVE};

typedef struct VoltageProfile{
	
    double *voltage;
    double *field;
    double *charge; 
    double *ionicCharge; 
    double *electronCharge; 
    double *holeCharge; 
	double polarization;
	double ionicPolarization;
	double debyePolarization;
	double *ionsPerLayer;
} VoltageProfile;

void es_registerSettings();
//void es_loadSettings(FILE *settingsfile);
//void es_saveSettings(FILE *settingsfile);
void VP_allocate(VoltageProfile *VP);
void VP_free(VoltageProfile *VP);
double es_Energy(Configuration *config);
void VP_save(Configuration *config);
void VP_combine(Configuration *configs, int numcombine,Configuration *outconfig);
void es_setup(double *newzValues, int newnumzValues, char *newionNames,int newnumIons,double *newionCharge, void (*traj_generator)(Trajectory *));
void es_Trajectory(Trajectory *traj);
double chargeBalance(crystal *crys, double ionicCharge, double interfaceCharge);

#endif
