#ifndef _STACD_H_
#define _STACD_H_

#include "sourcecontrol.h"

#define coulombconstant 14.40 //eV * angstrom / (electron charge)^2
#define permittivity 0.005524 //electron charge / angstrom / volt
#define majorityDensityCutoff 1e-12 //if density exceeds 1-majorityDensity, then it will be given as bottomed out
#define minorityDensityCutoff 1e-160  //if density is less than minorityDensityCutoff, it will bottom out.

//2 makes since for perovskites since the unit cell is this may layers thick, but other materials many require different strides
#define saveStride 1 //How many layers to be combined during saving
#define initializationTemperature 0.1
#define fermiLevelConvergence 0.001 //eV 
#define hbarSquaredOver2m 7.62 //eV/angstrom

enum electrostaticExternalStates {NO_BIAS, ION_ONLY, ION_FORWARD_FIELD, ION_REVERSE_FIELD};
enum voltageMode {VOLTAGE_SQUARE_WAVE, VOLTAGE_SINE_WAVE};
//In carrierProfile many data are doubled. When they are doubled, there is first and electron part and then a whole part forming a contiguous array/matrix for each. 
//This makes them easier to use with array/matrix based functions.
typedef struct CarrierProfile{
	
	double polarization;
	double chargeInjection;//negative for extraction
	double ionFieldEnergy;
	double ionIonEnergy;
	double electronEnergy;
	double holeEnergy;
	double totalEnergy;
	double quasiFermiLevelSplitting;
	double potentialDifference;
	double carrierGeneration;
	
	//These are calculated for every electronic state
	double *stateCharge;//This is one for a fully occupied defect, or scales uniformly for a band state. Goes electrons, holes for each electronic state
	double *quasiFermiLevels;//The quasiFermiLevel of the for each electron, hole state
	
	//All of these are calculated only for the bands/zValues
    double *potential;
    double *field;//Only the Azimuthal Field
    double *zValues;
    double *layerCharge; 
    double *ionicCharge; 
    double *electronCharge; 
    double *holeCharge; 
	double *ionsPerLayer;
	double *generation;
	
	//These are calculated for each defect
	double *effectiveCharge;//The charge of the defect considering it's nominal charge and trap state occupancy
	
	//These are reductions of the tensors to vectors in the carrier change equations
	double *stateRecombination;
	double *stateCurrent;
	double *stateCurrent2;
	double *stateChange;
	
	//These are reductions to per-zValue for saving and plotting
	double *layerCurrent;
	double *layerCurrent2;
	double *layerBandCharge;
	double *layerRecombination;
	double *layerChange;
	
	//These are pointers which are used for compact Data copying and saving
	double **scalars;
	int numScalars;
	double **layerVectors;
	int numLayerVectors;
	double dampingParameter; //Plays the same role as timeStep, but for equilibrium equations.
} CarrierProfile;

//Every part of the deviceStack is constant during the simulation. It just needs to be set at the beginning.
typedef struct DeviceStack{
double *zValues;
int numzValues;
int numStates;
int numDefectLevels;
int numDefectStates;

//It is assumed that only 1 type of defect are in the bands
double *defectDensities;
double *bandDefectLevels;
double *bandDefectCharge;

char *ionNames;
double *ionCharge;	
int numIonTypes;
int numIons;
double numAbsorberLayers;
int *ionAtoms;//The index is the ion as in the numIons above, the value is the atom in the crystal.
int *ionEles;//The index is the ion as in the numIons above, the value is the element in the ionNames, ionCharge.

int *materialTransitions;	//Lists the z index where the device transitions to each material.
int *absorbers; //This marks the material as an absorber or not an absorber.
char **materialNames; 
double *bandOffsets; //electron, hole for each material
int *numIonDefectLevels; //Number of defect levels for each ion type. 
double *defectLevels; //Energy level of each defect state relative to its corresponding band.
double *defectLevelCharge; //Associated charge upon filling defect level. 
double *dielectricConstants; //Dielectric constant for each material.
double *formulaUnits; //The number of formulation units of each material per z value. The most-important material will likely be integral with periodic boundaries, whilst the other materials will be fractional and only harbor charge carriers
double *thickness;
int numMaterials;

double volumePerUnitCell;

crystal *initializationCrystal;
} DeviceStack;

void CP_allocate(CarrierProfile *CP);
void CP_free(CarrierProfile *CP);
void CP_combineWeighted(Configuration *configs, int numCombine,Configuration *outconfig, double *weights);
void CP_save(Configuration *config);
void stacd_Trajectory(Trajectory *traj);
void stacd_registerSettings();
void stacd_setup(DeviceStack *newArchitecture, void (*traj_generator)(Trajectory *));
void printDeviceStack(DeviceStack *stack);
double stacd_Energy(Configuration *config);
double stacd_charge(Configuration *config);
double stacd_potential(Configuration *config, double chargeSum, double potentialDifference);
void stacd_chargeExtraction(Configuration *config);
double bandOccupancyToFermiLevel(double occupancy,double stateVolume, int invert);
double bandFermiLevelToOccupancy(double fermiLevel,double stateVolume, int iCharge);
double defectOccupancyToFermiLevel(double occupancy, int invert);
double defectFermiLevelToOccupancy(double fermiLevel, int iCharge);
double bandDefectOccupancyToFermiLevel(double fermiLevelGuess, double occupancy, double defectDensity, double volume, double temperature, double *defectLevels, int numDefectLevels, int invert);
double bandDefectFermiLevelToOccupancy(double fermiLevel, double defectDensity, double volume, double temperature, double *defectLevels, int numDefectLevels, int iCharge);
double overlap(double distance);
int zIndex(DeviceStack *architecture, double z);
int zMaterial(DeviceStack *architecture, int iz);
void fillIonAtoms(DeviceStack *architecture);
double stateStateDistance(crystal *crys, DeviceStack *architecture, int iState1, int iState2);
void compressZProfile(DeviceStack *architecture);
void CP_initialize(CarrierProfile *CP);
void stacd_stateToIndex(crystal *crys,int *iState, int *iBand, int *iCharge, int *iz, int *iMaterial, int *iIon, int *iEle, int *iAtom, int *iDL, int *ionDL, int *iDS);
void stacd_iDSToState(crystal *crys,int *iState, int *iBand, int *iCharge, int *iz, int *iMaterial, int *iIon, int *iEle, int *iAtom, int *iDL, int *ionDL, int *iDS);
void stacd_iIonToState(crystal *crys,int *iState, int *iBand, int *iCharge, int *iz, int *iMaterial, int *iIon, int *iEle, int *iAtom, int *iDL, int *ionDL, int *iDS);
void stacd_iBandToState(crystal *crys,int *iState, int *iBand, int *iCharge, int *iz, int *iMaterial, int *iIon, int *iEle, int *iAtom, int *iDL, int *ionDL, int *iDS);
void zIndexWeights(double *weight1, double *weight2, int *iNeighbor, double zValue, int iz);
#endif
