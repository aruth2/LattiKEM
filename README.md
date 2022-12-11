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

The settings were organized to emphasize the most-performance critical settings first. A list of the settings, the module they control, and their purpose are listed below:

dir - 
