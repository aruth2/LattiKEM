/*
 * perovskite.c - provides functions for creating and manipulating perovskite crystals
 * */

#include "perovskite.h"

crystal * perovskite_positiveLayer(char *element1,char *element2, double a)
{
	/*************************************************************************************
	* Creates an AX layer of a cubic perovskite structure
	*************************************************************************************/
	//printf("making a positive layer with elements %s %s\n",element1,element2);
	crystal *base = (crystal *)malloc(sizeof(crystal));
	base->numElements = 2;
	base->totalAtoms = 2;
	
	crys_allocate(base);
	*(base->totalEachSpecies) = 1;
	*(base->totalEachSpecies+1) = 1;

	double *lv = base->latticeVectors;
	*(lv+0) = a, *(lv+1) = 0,*(lv+2) = 0;
	*(lv+3) = 0,*(lv+4) = a,*(lv+5) = 0;
	*(lv+6) = 0, *(lv+7) = 0, *(lv+8)=0;


	sprintf(base->elements,"%s",element1);
	sprintf(base->elements+namelength,"%s",element2);
	sprintf(base->species,"%s",element1);
	sprintf(base->species+namelength,"%s",element2);
		
	
	double *bp = base->positions;
	*(bp+0) = *(bp+1) = *(bp+2) = 0;
	*(bp+3) = a/2, *(bp+4) = a/2, *(bp+5) = 0;
	
	
	return base;
}

crystal * perovskite_negativeLayer(char *element1,char *element2, char *element3, double a)
{
	/*************************************************************************************
	 * Creates a BX2 layer of a cubic perovskite structure. Seperate elements can be
     * chosen for the two different X-sites 
	*************************************************************************************/
	//printf("making a negative layer with %s, %s, and %s\n",element1,element2,element3);
	crystal *base = (crystal *)malloc(sizeof(crystal));
	base->numElements = 2;
	base->totalAtoms = 3;
	
	crys_allocate(base);
	*(base->totalEachSpecies) = 1;
	*(base->totalEachSpecies+1) = 2;

	double *lv = base->latticeVectors;
	*(lv+0) = a, *(lv+1) = 0,*(lv+2) = 0;
	*(lv+3) = 0,*(lv+4) = a,*(lv+5) = 0;
	*(lv+6) = 0, *(lv+7) = 0, *(lv+8)=0;


	sprintf(base->elements,"%s",element1);
	sprintf(base->elements+namelength,"%s",element2);
	sprintf(base->species,"%s",element1);
	sprintf(base->species+namelength,"%s",element2);
	sprintf(base->species+namelength*2,"%s",element3);
		
	
	double *bp = base->positions;
	*(bp+0) = a/2,*(bp+1) = a/2,*(bp+2) = 0;
	*(bp+3) = a/2, *(bp+4) = 0, *(bp+5) = 0;
	*(bp+6) = 0, *(bp+7) = a/2, *(bp+8) = 0;

	
	return base;
}


crystal * perovskite_newCrys(char *element1, char *element2, char *element3, char *element4,char *element5, double a)
{
	// Makes ABX3 perovskite
	//printf("Creating crystal %s %s %s %s %s\n",element1,element2,element3,element4,element5);
	crystal *crys1 = perovskite_positiveLayer(element1,element3,a);
	crystal *crys2 = perovskite_negativeLayer(element2,element4,element5,a);
	
	double offset[3];
	offset[0] = offset[1] = 0;
	offset[2] = a/2;
	
	crystal *crys = crys_combine(crys1,crys2,offset,0);
	*(crys->latticeVectors+8) = a;
	return crys;
}

//This function is unused.
int perovskite_plane(crystal *crys, int atom,double latticeConstant)
{
	//This will use the z value of the atom to determine if the atom is in an AC plane or a BC2 plane.
	
	double z = *(crys->positions+3*atom+2);
	double epsilon = 0.1;
	while (z > latticeConstant-epsilon)
	z-=latticeConstant;
	
	if(fabs(z) < epsilon)
	return 0;
	else
	return 1;
}

int perovskite_direction(crystal *crys, int i, int j, double latticeConstant)
{
	//0: B->B, 1: B->A, 2: A->B
	int planei = perovskite_plane(crys,i,latticeConstant);
	int planej = perovskite_plane(crys,j,latticeConstant);
	
	//printf("Atom %d is in plane %d and atom %d is in plane %d\n",i,planei,j,planej);


	
	if(planei == 0 && planej == 0)//This should only ever occur for H and I or H and Br
	{
	//printf("Somehow they are both in plane A\n");
	//crys_printAtomPosition(crys,i);
	//crys_printAtomPosition(crys,j);
	//printf("Distance between is %g\n",crys_atomDistance(crys,i,j));
	}
	if( planei == 1 && planej == 1)
	return 0;
	if( planei == 1 && planej == 0)
	return 1;
	if( planei == 0 && planej == 1)
	return 2;
	
	return -1;
	
}

void perovskite_initializeMolecules()
{
	crystal *FAmol = crys_readVASP("molecules/FA.vasp");
	crys_printAllAtoms(FAmol);
	crystal *MAmol = crys_readVASP("molecules/MA.vasp");
	crys_printAllAtoms(MAmol);

	crystal *molecules = (crystal *)malloc(2*sizeof(crystal));
	*(molecules)=*FAmol;
	*(molecules+1)=*MAmol;
	char *moleculeString = crys_elementString(2,FA,MA);
	crys_setupMolecules(moleculeString,2,molecules);
}
