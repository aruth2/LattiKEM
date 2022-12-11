#ifndef _CRYSTALNETWORK_H_
#define _CRYSTALNETWORK_H_

#include <pthread.h>
#include "crystal.h"

#define maxnnds 20
//The nearestNeighborDescriptor is used to reduce computational cost in certain instances. 
//cn_swap and cn_fillNetwork are commutative and both are costly. Therefore swaps can be performed
//on a crystal without an established network and the network is only filled at the end.
typedef struct NearestNeighborDescriptor {
double nndistance;
char *elementList;
int numEle;
} NearestNeighborDescriptor;

enum cn_hasNetwork {CN_NO_NETWORK, CN_HAS_NETWORK};

void cn_allocateCrystalNetwork(crystal *crys);
void cn_fillNetwork(crystal *crys, double nndistance, char *elementList, int numEle);
void cn_printAdjacencyList(crystal *crys);
void cn_nearestNeighbors(crystal *crys, int *neighbors, int *numneighbors, int index);
void cn_nthNearestNeighbors(crystal *crys, int *cn_nearestNeighborseighbors, int *numcn_nearestNeighborseighbors, int index, int n);
void cn_integratednthNearestNeighbors(crystal *crys, int *integratedcn_nearestNeighborseighbors, int *numintegratedcn_nearestNeighborseighbors, int index, int nmax);
void cn_shellComposition(crystal *crys, int index, double *shells, char *shellelements, int numshellelements, double *distances, int *numshells);
void cn_coordination(crystal *crys, char *element, char *coordElements, int numcoordElements, int *coord, int coordinationNumber, int maxCoordination);
void numberOfCoordNeighbors(crystal *crys, char *coordElement, char *countedElements, int numCountedElements);
void cn_swap(crystal *crys, int atom1, int atom2, int noNetwork);
void cn_fillFromnnd(crystal *crys, NearestNeighborDescriptor *nnd);
void cn_addLink(crystal *crys, int i, int j);
int cn_areConnected(crystal *crys, int i, int j);
void cn_bucketFillNetwork(crystal *crys, double nndistance, char *elementList, int numele, int bucketDirection);
#endif
