## LATTIce chemical Kinetic Evolution Modeling (LattiKEM) v2.0

Welcome to the second open-source release of LattiKEM! 

This release realized two substantial performance boosts. (1) We implemented the partial sums method developed by Belois and (2) we produced a GPU version which performs the self-consistent state occupancy using CUDA.

We have made it a priority to release LattiKEM to the public. As this is the second release, there most certainly will be bugs, errors, and omissions and the documentation is lacking. Please be patient as we work towards filling in the details. 

# Publications that used LattiKEM

 - 2018 https://doi.org/10.1021/acsenergylett.8b01369 Vacancy-Mediated Anion Photosegregation Kinetics in Mixed Halide Hybrid Perovskites: Coupled Kinetic Monte Carlo and Optical Measurements
 - 2021 https://doi.org/10.1021/acsenergylett.1c00790 Distinguishing models for mixed halide lead perovskite photosegregation via terminal halide stoichiometry
 - 2022 https://doi.org/10.1021/acsnano.2c10781 Excitation Intensity- and Size-Dependent Halide Photosegregation in CsPb(I0.5Br0.5)3 Perovskite Nanocrystals
 - 2023 https://pubs.acs.org/doi/10.1021/acs.jpcc.3c04708 Thermodynamic Band Gap Model for Photoinduced Phase Segregation in Mixed-Halide Perovskites
 - 2023 https://pubs.acs.org/doi/10.1021/acsnano.3c07165 Modeling the Photoelectrochemical Evolution of Lead-based, Mixed-halide Perovskites due to Photosegregation

An example simulation using LattiKEM from the last publication above is illustrated below. 

https://github.com/aruth2/LattiKEM/assets/8935880/bb96256a-af1b-4944-9e39-662d6c600e1c

Many of the core ideas of LattiKEM originated in an earlier publication on disintegration of a graphene oxide lattice, however almost all of the code was rewritten: https://doi.org/10.1038/ncomms14521

 - Notably, a significant part of the algorithms used in LattiKEM were reproduced in https://doi.org/10.1021/acsaem.1c00707 
 - New accelerations were made to some of the core algorithms in https://research.tue.nl/en/studentTheses/kmc-modelling-of-light-induced-halide-segregation-in-perovskite-s 

# Functionality
LattiKEM is able to simulate migration and reactions of chemical species on a lattice. The only requirement from the user is to provide a function (*Energy Model*) which can take any configuration of the lattice as an argument and will return the total energy. LattiKEM will then simulate the chemical evolution of a given initial configuration using the Kinetic Monte Carlo (KMC) method. 

While many researchers often write their own Monte Carlo simulations, these typically lack quality of life features which can greatly accelerate the scientific development process. LattiKEM was designed so that an *Energy Model* can be constructed with minimal code and will still enjoy all of the inherent features of LattiKEM. These features include: acceleration of the KMC algorithm through thread-based parallelism, fault tolerance and checkpointing, performing multiple simulations and combining results to quantify stochastic effects, as well as saving, loading, and managing data. 


# Contributing to LattiKEM

If you have any requests for features or any questions including anything which should be covered in the documentation, but is not there, please add an "Issue" on Github, describe the problem as best you can, and we will work to address it. You should get a response within 24 hours for any new issues.

If you work on any project using LattiKEM, we ask that you please fork a branch and provide your working code there. The code you write can be incorporated into LattiKEM and help future researchers. Prototyping should be performed within the "Energy Models", and resuable pieces are subsequently migrated into the "Core".

## Getting started
In order to run simulations of photosegregation in mixed-halide perovskites: 

Starting in version 2.0, the program exerts intelligent control over Non-Uniform Memory Access (NUMA). For this function, users will need to install the numa libraries. This can be done on Debian-based distros using: 

```
sudo apt install libnuma-dev
```

Go to the source directory and build the executables
```
make mixedhalide
```
The mixed halide executable should have been created in the src/Executables directory. 

To use the GPU version, instead run:
```
make mixedhalide_GPU
```

Run the program by passing it a settings file. Here we will use the Examples/mixedhalide/4x4x4 file

```
src/Executables/mixedhalide Examples/mixedhalide/4x4x4
```

The Examples have been recently revamped to be more instructive and showcase the usage and behavior of many key settings.

Further information on the "bandgap" energy model and the "mixedhalide" executable can be found here: https://github.com/aruth2/LattiKEM/tree/main/src/EnergyModels/Bandgap  

The program should run through stages: 
1. Loading settings
2. Initializing administrator and worker threads
3. Performing each trajectory
4. Post-processing

After running the program, the output directory should contain details of each trajectory including its initial structure, final structure, and every move that was made. The structures are in the commonly used '.xyz' format which can be read by Jmol and VMD. The energies of the trajectories over time is also listed. The average energies of the combined trajectory is also included in the root directory. Finally, the 'traj' directory includes spectral response functions over the course of the simulation including absorption spectra, photoluminescence spectra, and density of states. 

The settings were organized to emphasize the most-performance critical settings first. A list of the settings, the module they control, and their purpose are listed below:

- **dir** - from lattikem.c - value: ./4x4x4/ - target directory for simulations. If if does not exist, it will be created, but the parent directory must exist.
- **sizex/sizey/sizez** - from mixedhalide.c - value: 4 - determines the size of the supercell
- **numSteps** - from latticedynamics.c - value: 300 - how many KMC steps to perform
- **biasTurnOn** - from bandgap.c - value: 150 - what step the photocarriers are turned on, equivalent to beginning illumination in the experiment
- **biasSwitch** - from bandgap.c - value: 200 - what step the photocarriers are turned off, equivalent to ending illumination in the experiment
- **numExcitations** - from bandgap.c - value: 10 - how many photocarriers are placed in the supercell
- **vacancyRatio** - from mixedhalide.c - value: 0.01 - what fraction of ions are replaced by vacancies. Determines the number of available hops at each step of the trajectory.
- **job** - from lattikem.c - value: main - what tasks to perform. This determines the 300 ionic steps, but does not run postprocessing. 

Here the modules and their purpose are listed. They are organized from the deepest core modules to the most high level modules. All core modules are compiled into the Executables, but only modules within the same EnergyModel are compiled together.

### Core: https://github.com/aruth2/LattiKEM/tree/main/src/Core
1. **supp.c** - based support functions, list manipulations, string manipulations, and some basic math
2. **matrix.c** - a few basic matrix operations.
3. **settings.c** - loads the settings for all modules from a single settings file, also assigns default values.
4. **crystal.c** - defines a crystal as a set of atoms within a cell, organizes atoms by element, and includes methods for manipulating crystal objects.
5. **crystalnetwork.c** - defines the adjacency list which connects atoms in the crystal. This is the lattice which subsequent modules use for manipulation.
6. **perovskite.c** - defines a perovskite as a crystal and provides methods for making perovskite crystals.
7. **latticedynamics.c** - defines the changes allowed to the crystalnetwork, for instance a halide vacancy and a halide atom can swap places. Each reaction is given an energy barrier to affect its reaction rate. Defines the trajectory.
8. **pkmc.c** - manages a threadpool model of POSIX threads. Performs a parallel kinetic hop to advance the trajectory forward 1 step. This involves enumerating possible moves, creating jobs to evaluate each move, and determining the chosen move.
9. **lattikem.c** - manages trajectories, creates administrator threads to propagate each trajectory, includes methods for post-processing.

### EnergyModels:
- **bandgap.c** - assumes atoms in the lattice move to reduce the energy of charge carriers. Calculates the energy of the charge carriers using a cluster approximation.
- **mixedhalide.c** - implements bandgap.c, calls for a lattikem simulation of photosegregation in mixed-halide perovskites using the bandgap model
