#ifndef _GPUCLUSTER_H_
#define _GPUCLUSTER_H_

#include "lattikem.h"
#include <stdint.h>
#include <cuda.h>

void ACC_nearestNeighbors(crystal *crys, int *neighbors, int *numneighbors, int index);
void ACC_integratednthNearestNeighbors(crystal *crys, uint32_t *restrict integratednneighbors, int *restrict numintegratednneighbors, int *restrict shellatoms, int *restrict lastshellatoms, int index, int n);
void ACC_clusterApprox(int numGangs, int vectorLength, crystal *restrict crys, char *restrict coordElement, char *restrict countedElements, int numCountedElements, int *restrict coord, int coordinationNumber, int maxCoordination);
void  SetBit( uint32_t A[],  int k );
int TestBit( uint32_t A[],  int k );
void printBits(void const * const ptr,size_t const size);
void CPU_integratednthNearestNeighbors(crystal *crys, uint32_t *restrict integratednneighbors, int *restrict numintegratednneighbors, int *restrict shellatoms, int *restrict lastshellatoms, int index, int n);
void CPU_clusterApprox(crystal *restrict crys, char *restrict coordElement, char *restrict countedElements, int numCountedElements, int *restrict coord, int coordinationNumber, int maxCoordination);
void CUDA_clusterApprox(int numGangs, int vectorLength, crystal *restrict crys, char *restrict coordElement, char *restrict countedElements, int numCountedElements, int *restrict coord, int coordinationNumber, int maxCoordination);

#endif
