#ifndef _POLARIZATION_H_
#define _POLARIZATION_H_

#include "pmc.h"

#define coulombconstant 14.40 //eV * angstrom / (electron charge)^2
#define permittivity 0.005524 //electron charge / angstrom / volt
//coulomb constant = 1/(4*pi*permittivity)
#define ITERATION_LIMIT 100

enum electrostaticExternalStates {NO_BIAS, ION_ONLY, ION_FORWARD_FIELD, ION_REVERSE_FIELD};
enum voltageMode {VOLTAGE_SQUARE_WAVE, VOLTAGE_SINE_WAVE};

typedef struct FieldProfile{
	
    //In contrast to electrostatics.c and electrostatics.h, this is a vector field for every ion or dipole in the system.
    double *field;
    double *voltage;
    //expectation value of the dipole moment of each polar species. 
    double *moment;
    
    double *layerCharge;
    double *layerMoment;
    double *layerField;
    double *layerVoltage;
     
	double polarization;
	double *ionsPerLayer;
	double *dipolesPerLayer;
} FieldProfile;

void pz_loadSettings(FILE *settingsfile);
void pz_saveSettings(FILE *settingsfile);
void FP_allocate(FieldProfile *VP);
void FP_free(FieldProfile *VP);
double pz_Energy(Configuration *config);
void FP_save(Configuration *config);
void FP_combine(Configuration *configs, int numcombine,Configuration *outconfig);
void pz_setup(double *newzvalues, int newnumzvalues, char *newionNames,int newnumIons,double *newionCharge, char *newdipoleNames,  int newnumDipoles, double *newdipolePolarizability, void (*traj_generator)(Trajectory *));
void pz_Trajectory(Trajectory *traj);
#endif
