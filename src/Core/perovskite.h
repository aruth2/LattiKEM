#ifndef _PEROVSKITE_H_
#define _PEROVSKITE_H_

#include "crystal.h"

//Give molecules or defects names from the periodic table to enable visualization.
#define FA "Fe"
#define VX "Xe"
#define VA "V"
#define MA "Mo"

crystal * perovskite_positiveLayer(char *element1,char *element2, double a);
crystal * perovskite_negativeLayer(char *element1,char *element2, char *element3, double a);
crystal * perovskite_newCrys(char *element1, char *element2, char *element3, char *element4,char *element5, double a);
int perovskite_plane(crystal *crys, int atom,double latticeConstant);
int perovskite_direction(crystal *crys, int i, int j, double latticeConstant);
void perovskite_initializeMolecules();
#endif
