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

# TODO Organize and add comments to these functions

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



# Parallel Kinetic Monte Carlo

# LattiKEM

# Support

# Matrix

# Settings

# Source Control

# Perovskite
