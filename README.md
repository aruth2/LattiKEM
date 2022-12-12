# LattiKEM
Simulate chemical reactions on a lattice using Kinetic Monte Carlo Methods

##Getting started
Go to the source directory and build the executables
```
make mixedhalide
```
The mixed halide executable should have been created in the src/Executables directory. 

Run the program by passing it a settings file. Here we will use the Examples/mixedhalide/medium_fewsteps_4runs file

```
src/Executables/mixedhalide Examples/mixedhalide/medium_fewsteps_4runs
```

The program should run through stages: 
1. Loading settings
2. Initializing administrator and worker threads
3. Performing each trajectory
4. Post-processing

After running the program, the output directory should contain details of each trajectory including its initial structure, final structure, and every move that was made. The structures are in the commonly used '.xyz' format which can be read by Jmol and VMD. The energies of the trajectories over time is also listed. The average energies of the combined trajectory is also included in the root directory. Finally, the 'traj' directory includes spectral response functions over the course of the simulation including absorption spectra, photoluminescence spectra, and phase fractions. 

The settings were organized to emphasize the most-performance critical settings first. A list of the settings, the module they control, and their purpose are listed below:

- dir - from lattikem.c - value: ./test_run/ - target directory for simulations. If if does not exist, it will be created, but the parent - directory must exist.
- numRuns - from lattikem.c - value: 4 - the number of trajectories to perform. Each trajectory uses it's own administrator thread and each administrator thread submits jobs for worker threads to complete.
- numSteps - from latticedynamics.c - value: 50 - how many KMC steps to perform for each trajectory
- hardcodedMaxThreads - from lattikem.c - value: 8 - sets the number of worker threads instead of using sysconf to obtain number of processors.
- vacancyRatio - from mixedhalide.c - value: 0.01 - what fraction of ions are replaced by vacancies. Determines the number of available hops at each step of the trajectory.
- sizex/sizey/sizez - from mixedhalide.c - value: 6 - determines the size of the supercell
coordinationNumber - from bandgap.c - value: 3 - determines how deep to perform a Breadth-First Search in the cluster approximation.
- biasTurnOn - from bandgap.c - value: 0 - what step the photocarriers are turned on, equivalent to beginning illumination in the experiment
- biasTurnOff - from bandgap.c - value: 45 - what step the photocarriers are turned off, equivalent to ending illumination in the experiment
- numExcitations - from bandgap.c - value: 5 - how many photocarriers are placed in the supercell
- thermalDistribution - from bandgap.c - value: fermi - whether to use the Maxwell-Boltzmann or Fermi-Dirac distribution for the photocarriers
- numOpticsEnergies - from bandgap.c - value: 200 - how many points to use when exporting optical data. Has a strong effect on post-processing time.
- saveDefaults - from settings.c - value: 1 - whether to list settings that are on their default value in the outputted settings file.

Here the modules and their purpose are listed. They are organized from the deepest core modules to the most high level modules. All core modules are compiled into the Executables, but only modules within the same EnergyModel are compiled together.

###Core:
1. supp.c - based support functions, list manipulations, string manipulations, and some basic math
2. matrix.c - a few basic matrix operations.
3. settings.c - loads the settings for all modules from a single settings file, also assigns default values.
4. crystal.c - defines a crystal as a set of atoms within a cell, organizes atoms by element, and includes methods for manipulating crystal objects.
5. crystalnetwork.c - defines the adjacency list which connects atoms in the crystal. This is the lattice which subsequent modules use for manipulation.
6. perovskite.c - defines a perovskite as a crystal and provides methods for making perovskite crystals.
7. latticedynamics.c - defines the changes allowed to the crystalnetwork, for instance a halide vacancy and a halide atom can swap places. Each reaction is given an energy barrier to affect its reaction rate. Defines the trajectory.
8. pkmc.c - manages a threadpool model of POSIX threads. Performs a parallel kinetic hop to advance the trajectory forward 1 step. This involves enumerating possible moves, creating jobs to evaluate each move, and determining the chosen move.
9. latticekem.c - manages trajectories, creates administrator threads to propagate each trajectory, includes methods for post-processing.

###EnergyModels:
bandgap.c - assumes atoms in the lattice move to reduce the energy of charge carriers. Calculates the energy of the charge carriers using a cluster approximation.
mixedhalide.c - implements bandgap.c, calls for a lattikem simulation of photosegregation in mixed-halide perovskites using the bandgap model
