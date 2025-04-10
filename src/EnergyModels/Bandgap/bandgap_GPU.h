#ifndef _BANDGAP_GPU_H_
#define _BANDGAP_GPU_H_

//CUDA .cu files are compiled as C++ code. This wrapper allows inclusion of c header in a C++ code
//extern "C" {
#include "bandgap.h"
//}
#include <cuda_runtime.h>
#include <cuda.h>
#include <nvtx3/nvToolsExt.h>


__device__ void fermiIteration(double *energyLevels, double fermiGuess, double temperature, int numEnergyLevels, double *weights, double *energyWeights, double *weightSum);
__global__ void fermiKernel(double *energyLevels, double *fermiEnergy, double fermiGuess, double currentExcitations, double fermiConvergence, double temperature, int numEnergyLevels,  double *weights, double *energyWeights, double *d_energy, double *weightSum);
__host__ void bgg_runKernel(OptoelectronicState *OS, int copyWeights, int num_states, double currentExcitations);
__host__ void bgg_allocate(Configuration *config, int num_states);
__host__ void bgg_free(OptoelectronicState *OS);
__host__ void bgg_registerSettings();
__host__ void bindGPU(int threadnumber);

#endif
