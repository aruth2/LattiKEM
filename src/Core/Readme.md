## Core Libraries of LattiKEM
This document is a guide for those who wish to create their own Energy Model in LattiKEM. This is a difficult undertaking that requires understanding what the core functionality of LattiKEM will provide and what the programmer must provide in order to simulate the energy model. Learning how to use an existing energy model is a much easier task, and each of them provides their own self-contained guides.

# Primary Libraries
 - [Crystal](#crystal)
 - [Crystal Network](#crystal-network)
 - [Lattice Dynamics](#lattice-dynamics)
 - [Parallel Kinetic Monte Carlo](#parallel-kinetic-monte-carlo)
 - [LattiKEM](#lattikem)
# Support Libraries
 - [Support](#support)
 - [Matrix](#matrix)
 - [Settings](#settings)
 - [Source Control](#source-control)
 - [Perovskite](#perovskite)



# Crystal

The crystal module defines a cell of an infinite lattice. The crystal module defines the "crystal" structure:

```
typedef struct crystal
{
	char *elements; //List of the elements used in the crystal
	int *totalEachSpecies; //The number of atoms of each element
	int numElements;
	double *positions; //List of x,y, and z values of every atom in the crystal
	char *species; //The element of each atom
	int totalAtoms;
	double latticeVectors[9];
	crystalnetwork *network;
	
} crystal;
```

A crystal contains a list of atoms and their positions. Periodic boundary conditions are assumed for the crystal via the lattice vectors. When creating a "crystal", the crystal->totalAtoms attribute should be set first, then call crystal_allocate(crys). This will allocate the arrays so the data can be populated.

The crystal module also contains numerous methods for creating and manipulating crystals. A full listing of the functions of the crystal module is below:

```
void crys_allocate(crystal *crys);
void crys_free(crystal *crys);
crystal * crys_multiply(struct crystal *base, int m1, int m2, int m3);
crystal *crys_duplicate(crystal *crys);
crystal * crys_combine(struct crystal *crys1, struct crystal *crys2, double *offset, int addvectors);
crystal * crys_linearInterpolation(crystal *crys1, crystal *crys2, int numimages);
void crys_removeAtom(struct crystal *crys, int removalsite);
int crys_atomsOfElement(crystal *crys, char *element);
int crys_elementOffset(crystal *crys, char *element);
void crys_elementBoundsArray(crystal *crys, int *offsetArray);
int crys_elementCount(crystal *crys, char *element);
char * crys_elementString( int num, ... );
int crys_elementInString(char *elementlist, int numElementsinlist,char *element);
double crys_atomDistance(crystal *crys, int index1, int index2);
int crys_atomDirection(crystal *crys, int i, int j);
void crys_printAtomPosition(crystal *crys, int atom);
void crys_makexyz( crystal *crys,char *name);
void crys_makeVASP(struct crystal *crys, char *name, int cart);
void jmol(char *filename);
void vmdRender(char *dir, int id);
struct crystal * crys_readVASP(char *filename);
void crys_interatomicDistances(crystal *crys, double *intercrys_atomDistances);
double minsetdistances(double *intercrys_atomDistances, int numtoremove,int *atomnumbers,int totalAtoms);
crystal * crys_maxDefectSpacing(crystal *inputcrys, char **atoms,int numremove);
void cn_free(crystalnetwork *network);
crystalnetwork *crys_duplicateNetwork(crystal *crys);
void makeimage(crystal *crys,char *dir, int id);
void crys_replaceRandomAtoms(crystal *crys, int numReplace, char *replaceElements, int numReplaceElements, char *newElement);
void crys_printElements(crystal *crys);
char * crys_formatElementString(char *elementString, int numElements);
void crys_zValues(crystal *crys, double *list, int *numitems);
void crys_setupMolecules(char *newMolecules, int newNumMolecules, crystal *newMolecularSubstituitions);
void crys_printAllAtoms(crystal *crys);
crystal *crys_readxyz(char *name);
double crys_atomVector(crystal *crys, int iatom1, int iatom2, double *r);
crystal *crys_simpleCubic(char *element, double latticeConstant);
void crys_addAtom(struct crystal *crys, char *element, double x, double y, double z);
void crys_allocatedSize(crystal *crys);
int crys_elementIndex(crystal *crys, char *element);
```


# Crystal Network

The crystal network module adds graph-based network functionality to a crystal. In the crystal network, adjacent atoms are linked to one another. Each atom has it's own adjacency list and two atoms are linked if they appear in each other's adjacency lists. Storing adjacency information in this way is a means to contain the computational scaling of chemistry codes. The number of pairwise interactions scales quadratically with the number of atoms in the system, but the number of neighboring pairs scales linearly. 

The crystal network is defined in crystal.h. 

```
typedef struct crystalnetwork {
int *adjacencyList; //Indicies are iAtom * maxconnections + jNeighbor
int *numAdjacent;
} crystalnetwork;
```

Adjacency information could represent a chemical bond, but it could also represent any other information needed to evaluate the interaction between atoms or the possible chemical evolution of the crystal.

Filling of adjacency information is aided by a structure which describes what elements can be linked to which other elements, and the maximum distance between them. The "NearestNeighborDescriptor" is defined in crystalnetwork.h:

```
typedef struct NearestNeighborDescriptor {
double nndistance;
char *elementList;
int numEle;
} NearestNeighborDescriptor;
```

Possibly the most important function in the crystal network module is cn_swap which performs the delicate task of swapping two atoms on the lattice and reconstituing their network adjacency information.

```
void cn_swap(crystal *crys, int atom1, int atom2, int noNetwork);
```

a full list of functions in the crystal network module is provided below:

```
void cn_allocateCrystalNetwork(crystal *crys);
void cn_fillNetwork(crystal *crys, double nndistance, char *elementList, int numEle);
void cn_printAdjacencyList(crystal *crys);
void cn_nearestNeighbors(crystal *crys, int *neighbors, int *numneighbors, int index);
void cn_nthNearestNeighbors(crystal *crys, nneighbors, int *numnneighbors, int index, int n);
void cn_integratednthNearestNeighbors(crystal *crys, int *integratednneighbors, int *numintegratednneighbors, int index, int nmax);
void cn_shellComposition(crystal *crys, int index, double *shells, char *shellelements, int numshellelements, double *distances, int *numshells);
void cn_coordination(crystal *crys, char *element, char *coordElements, int numcoordElements, int *coord, int coordinationNumber, int maxCoordination);
void numberOfCoordNeighbors(crystal *crys, char *coordElement, char *countedElements, int numCountedElements);
void cn_swap(crystal *crys, int atom1, int atom2, int noNetwork);
void cn_fillFromnnd(crystal *crys, NearestNeighborDescriptor *nnd);
void cn_addLink(crystal *crys, int i, int j);
int cn_areConnected(crystal *crys, int i, int j);
void cn_bucketFillNetwork(crystal *crys, double nndistance, char *elementList, int numele, int bucketDirection);
void cn_allocatedSize(crystal *crys);
```

# Lattice Dynamics

Dynamic processes in the crystal are managed by the lattice dynamics module. The lattice dynamics module provides means for associating changes in energy with specific moves that could be made in the lattice. At present, the only allowed moves are a swap betweeen two adjacent atoms on the lattice. 

Pairwise elemental energy barriers are defined in the lattice dynamics structure:

```
typedef struct LatticeDynamics{
    //Describes the allowed moves
    char *hopPairs;
    double *hopPairEnergies;
    int numHopPairs;

} LatticeDynamics;
```

The lattice dynamics module provides means for enumerating and comparing different possible moves. The structure responsible for this is called "MarkovChainMonteCarlo". MarkovChainMonteCarlo is defined in latticedynamics.h:

```
typedef struct MarkovChainMonteCarlo{
	
    int *moves; //List of atom swap pairs
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
```

The MCMC object stores the list of all possible moves from the current crystal. It has places for storing the energy differences from each move as well as arrays for the data generated in the kinetic monte carlo algorithm: move rates and probability ranges.

An MCMC object can be initialized from a crystal using a description of possible moves (a LatticeDynamics object) using the LD_list moves function:

```
void LD_listMoves(crystal *crys, LatticeDynamics *LD, MarkovChainMonteCarlo *mcmc)
```

In higher level modules, the crystal is passed around as a data packet containing the crystal, external conditions applied to the crystal, and derived properties of the crystal. The structure of this data packet is called "Configuration" and is defined in latticedynamics.h

```
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
```

Finally, the lifetime of the crystal throughout the simulation is described by the "Trajectory" structure. The trajectory structure contains the present information about the crystal's current state as well as past information about what moves were performed to evolve the crystal. It also contains kinetic information including the time and total energy at each step of the evolution.

```
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
```

# Parallel Kinetic Monte Carlo

Having defined a KMC step with the MCMC structure, we now perform the KMC step. The KMC step is performed using the first level of parallism in LattiKEM: parallelization over the possible moves that could be performed in the KMC step. The step itself is performed in the traj_parallelKineticHopping function.

```
void traj_parallelKineticHopping(Trajectory *traj);
```

while the computationally-costliest part, evaluating the energy of a single Configuration, is offloaded to worker threads in a thread pool. pkmc.c contains functions for intializing the worker threadpool, adding jobs to the queue, and the operation of worker threads. 

A single job for a worker thread is defined in the kineticJob structure:

```
//A kinetic job is a single task which any worker thread can take on and complete. Following completion of the task,
//The energy variable is filled in, and if requested the data is also filled in. The worker thread may be tasked
//with saving the data and also possibly the crystal following its completion of the job.
typedef struct kineticjob{
	int iTraj;
	int step;
	double *energy;
	Configuration *dataSwap;
	int additionalMoveAtom1, additionalMoveAtom2;
} kineticjob;
```

The full list of functions in pkmc.h is below:

```
void traj_parallelKineticHopping(Trajectory *traj);
void threadwait();
void energyjob(double *energy, int iTraj, int step,int additionalMoveAtom1, int additionalMoveAtom2, Configuration *dataswap);
void kineticthread(workerData *wd);
void startWorkerThreadPool(int numthreads, Trajectory *newtraj, int newNumTrajectories);
void traj_parallelInitialization(Trajectory *trajectories, int numRuns, int maxThreads, char *dir, int jobType);
```

## LattiKEM
The LattiKEM module performs a Markov-Chain Monte Carlo simulation of ion migration on a lattice using either the kinetic Monte Carlo or Metropolis Monte Carlo algorithms. 
The LattiKEM module provides a range of quality of life features to these simulations including: repeating simulations **numRuns** times and combining the results, fault tolerance via saving and loading of simulation state, and exposing 2 levels of thread-based parallelism via calls to pkmc.

# Fault Tolerance, saving and loading trajectories 
Each running trajectory will be checkpointed every **savePeriod** steps. When the saving is performed all threads must be paused, so saving too frequently will negatively impact runtime. On the other hand, frequent saving can preserve data in the event the program is terminated for any reason. LattiKEM is generally highly fault tolerant because the amount of data that needs to be preserved is minimal. This data only consists of the initial crystal structure, the ion swaps performed at each step, and time series data including the energy, time, and time step.
The settings and simulation-wide data is stored in the highest level directory given by the **dir** variable. Each trajectory is stored in its own subdirectory using a three digit zero-indexed name. The trajectory will have an “initial.xyz” file, a “moves” file, an “energy” file, and files for each type of **seriesData** defined for the energy model. For example, the bandgap module describes photocarrier energy and interatomic energy as series data. 

During initialization, the directory for each trajectory is examined and previous crystal structures and moves along with their associated time series data are loaded. Each trajectory will then proceed from the last loaded move unless the **startAtStep** setting is provided in which case all trajectories will begin from the specified step. Loading of previous trajectories can be skipped using “job = restart”. Note that changing parameters of the energy model prior to loading previously generated trajectories can create an incompatibility between the trajectories and the expected simulation. Besides the trajectories that were found and loaded, additional trajectories will be created until the total number equals **numRuns**. In this way, additional runs can be added to a simulation and LattiKEM will only process the newly requested trajectories.

Each trajectory will be ran until numSteps have completed. The external conditions applied to the simulation during each step are defined within the corresponding energy models. Note that at present, the Metropolis Monte Carlo option has been disabled so only “hopMode = kinetic” should be used.


# Post Processing
Both the simulation and the trajectories may have additional data associated with specific steps or step ranges and these data are stored in the “/traj” directories. The “/traj” directories are for larger files that are created during post processing. If “saveCrys = 1”, a .xyz file will be generated for every step of the zeroth trajectory in the “000/traj” directory. This can use an enormous amount of disk space. The bandgap module stores spectra and spatially-distributed charge carrier density in the “/traj” directories.  
Multiple trajectories are averaged during post processing to produce the response of the simulation. This averaging is only performed if “saveData=1”. The energy model is called to calculate “Configurations” at specific steps of specific trajectories. The energy model must provide a function to the latticeDynamics module which performs weighted averaging of “Configurations”. The trajectories can be averaged using either step-based alignment using “averagingMethod=step” or time-based alignment using “averagingMethod=time”. The number of steps which are combined to produce response data in the /traj directories is given by **averagingSteps** 

In step-based alignment, the averaged trajectory at each step is simply the weighted average of each constituent trajectory at the same corresponding step. Step-based alignment is simple to use, however it is not fully consistent with the macroscopic view providing by the kinetic Monte Carlo method. This is because KMC associates a time with each step, and the amount of time between steps is variable. Therefore, LattiKEM also provides the time-based averaging method. In the time-based averaging method, first the times for the trajectories are shifted so that the step given by “equivocationStep” is time zero. This provides a means to align the trajectories with a change in the external conditions of the simulation such as the flipping of a switch or turning on a light. 
Next, the time series from the trajectory with the greatest timespan is used as the basis to construct an effective average trajectory. The timeseries of each trajectory is compared to the basis time series so that the two steps surrounding the targeted time are identified and each Configuration is calculated and weights are generated for linear interpolation. This process is repeated for every trajectory and the results are averaged giving equal total weight to each trajectory. Finally, the steps of the effective average trajectory are averaged by **averagingSteps**. 

# Listing of all functions in the LattiKEM module:
```
void saveTrajectoryCrystals(Trajectory *traj);
void performTrajectory(Trajectory *traj);
void simulateTrajectories();
void stepAverage(Trajectory *trajectories);
void timeAverage(Trajectory *trajectories);
void postProcess(Trajectory *trajectories);
void lattikem_registerSettings();
void lattikem_startupMessage();
char *getDir();
```

# Support

# Matrix

## Settings
The settings module implements a single settings file for loading of all global settings in modules compiled together with LattiKEM. Each setting should be defined globally within the .c file of corresponding module and registered with the settings module during runtime. 

Before registering the first setting, the allocateSettings() function should be called to create internal space for the settings module.

Registering a setting consists of calling the register function of the appropriate data type, passing it a descriptive name which will be used to identify the setting within the settings file, passing a pointer to the settings variable, and providing a default value for the setting variable. 

Registering an enum type additionally requires providing the number of different options available and descriptive names for each option. During loading the enum type will be used to fill an integer value into the associated settings variable enumerated from the descriptive names. 

Once all settings have been registered, the loadSettings(FILE *filename) function should be called passing in the settings file. Then all settings variables from all modules will be loaded simultaneously. If a setting is listed multiple times in the settings file, only the first instance will be used.

The settings can be saved using the saveSettings(FILE *filename) function. Saving the settings in the same directory as the output of a simulation is good practice for ensuring repeatability of a simulation. Only settings which differ from their default values will be saved unless the setting saveDefaults = 1 is included in the input settings file. 

The settings module allows for comments in the settings file. All text after any of the characters "!@#$%^&*`~';:" will be treated as a comment. These characters cannot be used in the descriptor for any settings variable. 

# Example implementation of the settings module from lattikem.c
In the preamble of the lattikem.c file, the savePeriod setting is defined globally:
    int savePeriod;
In the lattiekem_registerSettings() function, the savePeriod settings variable is registered using a descriptor of “savePeriod” and a default value of 1000:
    registerInt(&savePeriod,"savePeriod",1000);

# Listing of all functions in the settings module:
```
void allocateSettings();
void registerInt(int *pointer, char *descriptor, int defaultValue);
void registerDouble(double *pointer, char *descriptor, double defaultValue);
void registerString(char *pointer, char *descriptor, char *defaultValue);
void registerEnum( int num, int *pointer, char *descriptor, int defaultValue,  ... );
void loadSettings(FILE *settingsFile);
void saveSettings(FILE *settingsFile);
```


# Source Control

# Perovskite
