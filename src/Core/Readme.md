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

# Lattice Dynamics

# Parallel Kinetic Monte Carlo

# LattiKEM

# Support

# Matrix

# Settings

# Source Control

# Perovskite
