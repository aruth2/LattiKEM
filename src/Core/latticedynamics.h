#ifndef _LATTICEDYNAMICS_H_
#define _LATTICEDYNAMICS_H_

#include "crystalnetwork.h"
#include "settings.h"

enum breakpoints {ALL_BREAKPOINTS,LOGARITHMIC_BREAKPOINTS,BIAS_BREAKPOINTS};
enum trajectoriestate {TRAJECTORY_STATE_FINISHED, TRAJECTORY_STATE_RUNNING};
enum hopModes {HOP_MODE_METROPOLIS, HOP_MODE_KINETIC};

//A LatticeDynamics object describes allowed changes which can be made to a crystal
typedef struct LatticeDynamics{
    //Describes the allowed moves
    char *hopPairs;
    double *hopPairEnergies;
    int numHopPairs;

} LatticeDynamics;


typedef struct MarkovChainMonteCarlo{
	
	int *moves;
    double *moveBarriers;
    int numMoves; 
    double *moveEnthalpies;
    
	int chosenMove;
	double initialenergy;
	
	//This is only used in kmc
    double *moveRates;
    double *moveProbabilityRanges;
    double timestep;
    double rate;
	
} MarkovChainMonteCarlo;

typedef struct Configuration{
	crystal *crys;
	int savedata;
    char dataFileName[1000];
    char crysFileName[1000];
    
    double energy;
    void *data;
    int enthalpyState;//0 means no energy biasing is used. 
    double *externalConditions;//A double precision number which is passed to determine how the energy is calculated. //E.g. this could be an applied voltage
	//Series data is a double precision number for each configuration. 
	//These are attached to the trajectory and saved by the trajectory itself
    double *seriesData;
} Configuration;

//A Trajectory object contains a list of all moves performed during simulation.
//It also describes breakpoints for when data is saved and the steps when external conditions are changed
//It is started as a thread on the performTrajectory function 
typedef struct Trajectory{
	
    int step;    
    //int numSteps;
	int numExternalConditions;

	int iTraj;

	double *energies;
	//This is used to encode information which gets passed to the energyfunction to describe what energetic terms are used
	//E.g. turn on/turn off electric field. If the enthalpyState is zero then no enthalpy will be calculated.
	
	//Breakpoints are when data is saved
	//They may become obsolete if a postprocessing utility is made which can quickly produce the data at any step/configuration
	int *breakpoints;
	int numbreakpoints;
	
	double *timeseries;
	double *timeSteps;
	double *externalConditions;
	//Series data is a double precision number for each configuration. 
	//These are attached to the trajectory and saved by the trajectory itself
	double *seriesData;
	int numSeriesData;
	char **seriesDataNames;
	
	//double temperature; 
    //int hopMode;

    char *dir;
    int *selectedMoves;
    
    MarkovChainMonteCarlo *mcmc;
    LatticeDynamics *LD;
    NearestNeighborDescriptor nnd[maxnnds];
    int numnnds;
    
    crystal *crys;
    //This is only used in metropolis
    int numrejectedhops;
	int state;
	int willBeLoaded;
} Trajectory;

void LD_setup(double (*newEnergyFunction)(Configuration *), void (*newSavingFunction)(Configuration *), void (*newbreakpoint_allocator)(void *),
 void (*newbreakpoint_freer)(void *), void (*newAveragingFunction)(Configuration *,int, Configuration *), void (*newWeightedAveragingFunction)(Configuration *,int, Configuration *, double *), 
 void (*newtraj_generator)(Trajectory *), int newDataSize);

void LD_listMoves(crystal *crys, LatticeDynamics *LD, MarkovChainMonteCarlo *mcmc);
void mcmc_printMoves(crystal *crys, MarkovChainMonteCarlo *mcmc);
void mcmc_allocate(MarkovChainMonteCarlo *mcmc, int maxMoves);

void config_allocate(Configuration *config, Trajectory *traj);
void config_free(Configuration *config);
void config_save(Configuration *config, Trajectory *traj, int step, int maxSteps);
double config_energy(Configuration *config);
void config_average(Configuration * configs,int numConfigs, Configuration *outConfig);
int config_toBeSaved(Configuration *config, Trajectory *traj, int step, int maxSteps);
void config_savecrys(Configuration *config, Trajectory *traj, int step, int maxSteps);

//void traj_loadSettings(FILE *settingsfile);
//void traj_saveSettings(FILE *settingsfile);
void traj_registerSettings();
void traj_logarithmicBreakpoints(Trajectory *traj);
void traj_allBreakpoints(Trajectory *traj);
void traj_biasBreakpoints(Trajectory *traj);
void traj_crysAtStep(crystal *crys, int *usedMoves, int initialStep, int targetStep, int hasNetwork);
void traj_average(Trajectory *trajectories, int numRuns, Trajectory *outTraj);
void traj_weightedAverage(Trajectory *trajectories, int numRuns, Trajectory *outTraj, int longestRun);
void traj_loadSeries(Trajectory *traj);
int traj_load(Trajectory *traj);
void traj_saveSeries(Trajectory *traj);
void traj_free(Trajectory *traj);
void traj_saveEnergy(Trajectory *traj);
void traj_loadEnergy(Trajectory *traj);
void traj_saveMoves(Trajectory *traj);
void traj_loadMoves(Trajectory *traj);
void traj_saveSeriesData(Trajectory *traj);
void traj_loadSeriesData(Trajectory *traj);
//void traj_copySettings(Trajectory *dest);
void traj_initialize(Trajectory *traj);
Trajectory * traj_getSettings();
void config_weightedAverage(Configuration * configs,int numConfigs, Configuration *outConfig, double *weights);
void traj_identifySelf(Trajectory *traj, int iTraj, char *dir, int willBeLoaded);
int getNumSteps();
int getHopMode();
double getTemperature();
#endif
