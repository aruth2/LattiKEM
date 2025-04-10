
/* crystal.c - A header file that defines a crystal and contains functions for manipulating it
 * 
 * */
#include "crystal.h"
char *molecules;
int numMolecules=0;
crystal *molecularSubstitutions;

void crys_allocate(crystal *crys)
{
	crys->species = (char *)calloc(crys->totalAtoms * namelength, sizeof(char));
	crys->positions = (double *)calloc(3*crys->totalAtoms,sizeof(double));
	crys->totalEachSpecies = (int *)calloc(maxelements,sizeof(int));
	crys->elements= (char *)calloc(maxelements*namelength,sizeof(char));
	crys->network = (crystalnetwork *)calloc(1,sizeof(crystalnetwork));
}

void crys_allocatedSize(crystal *crys)
{
	int size;
	size = crys->totalAtoms * namelength* sizeof(char) + 3*crys->totalAtoms*sizeof(double) + maxelements*sizeof(int)
	+ maxelements*namelength*sizeof(char) + 1*sizeof(crystalnetwork);
	printf("Expected size of crystal is %d\n",size);
}

void crys_free(crystal *crys)
{
	free(crys->species);
	free(crys->positions);
	free(crys->totalEachSpecies);
	free(crys->elements);
	
	if(crys->network != NULL)
	cn_free(crys->network);
	free(crys);
}

crystal * crys_multiply(struct crystal *base, int m1, int m2, int m3)
{
	/*********************************************************************************/
	/*
	 * */
	 if(base == NULL)
	 {
		 printf("Crystal not initialized, cannot multiply\n");
		 return NULL;
	 }
	 
	crystal *output = (crystal *)malloc(sizeof(crystal));
	output->totalAtoms = m1*m2*m3*base->totalAtoms;
	output->numElements = base->numElements;
	crys_allocate(output);

	int element;
	for(element = 0;element<output->numElements;element++)
	{
		strcpy(output->elements+namelength*element,base->elements+namelength*element);
	}
		
	int species;
	for(species = 0;species<base->numElements;species++)
	*(output->totalEachSpecies+species) = *(base->totalEachSpecies+species)*m1*m2*m3;
	
	double *lv = output->latticeVectors;
	double *lvi = base->latticeVectors;
	*(lv+0) = m1 * *(lvi+0);
	*(lv+1) = m1 * *(lvi+1);
	*(lv+2) = m1 * *(lvi+2);
	*(lv+3) = m2 * *(lvi+3);
	*(lv+4) = m2 * *(lvi+4);
	*(lv+5) = m2 * *(lvi+5);
	*(lv+6) = m3 * *(lvi+6);
	*(lv+7) = m3 * *(lvi+7);
	*(lv+8) = m3 * *(lvi+8);
	
	
	int i,j,k,l,h;
	int atomcount = 0;
	for(h=0;h<base->totalAtoms;h++)//This looping system assumes that in the basis the atoms are arranged in order I.E. Mo S S not S Mo S and this order matches the order in base->elements
	for(i=0;i<m1;i++)
	for(j=0;j<m2;j++)
	for(k=0;k<m3;k++)
		{
		for(l=0;l<3;l++)
		*(output->positions+3*atomcount + l) = *(lvi+0+l) * i + *(lvi+3+l) * j + *(lvi+6+l) *k + *(base->positions+3*h+l); 
		sprintf(output->species+namelength*atomcount,"%s",base->species + namelength*h);
		atomcount++;
		}
		
	return output;
}

crystal *crys_duplicate(crystal *crys)
{
	//This duplicates a crystal and its crystalnetwork
	crystal *outcrys = (crystal *)malloc(sizeof(crystal));
	outcrys->totalAtoms = crys->totalAtoms;
	outcrys->numElements = crys->numElements;
	crys_allocate(outcrys);

	int element;
	for(element = 0;element<crys->numElements;element++)
	{
		strcpy(outcrys->elements+namelength*element,crys->elements+namelength*element);
	}
		
	int species;
	for(species = 0;species<crys->numElements;species++)
	*(outcrys->totalEachSpecies+species) = *(crys->totalEachSpecies+species);
	
	double *lv = outcrys->latticeVectors;
	double *lvi = crys->latticeVectors;
	*(lv+0) = *(lvi+0);
	*(lv+1) = *(lvi+1);
	*(lv+2) = *(lvi+2);
	*(lv+3) = *(lvi+3);
	*(lv+4) = *(lvi+4);
	*(lv+5) = *(lvi+5);
	*(lv+6) = *(lvi+6);
	*(lv+7) = *(lvi+7);
	*(lv+8) = *(lvi+8);
	
	
	int l,h;
	int atomcount = 0;
	for(h=0;h<crys->totalAtoms;h++)
	{
	for(l=0;l<3;l++)
	*(outcrys->positions+3*atomcount + l) = *(crys->positions+3*h+l); 
	strcpy(outcrys->species+namelength*atomcount,crys->species + namelength*h);
	atomcount++;
	}

	outcrys->network = crys_duplicateNetwork(crys);
	return outcrys;
}

crystal * crys_combine(struct crystal *crys1, struct crystal *crys2, double *offset, int addvectors)
{
	/*******************************************************************/
	/* combines two crystals together with offset x,y,z                */
	/*                                                                 */
	/*******************************************************************/
	if(crys1 == NULL || crys2 == NULL)
	{
		printf("Crystal not initialized, cacn_nearestNeighborsot crys_combined\n");
		return NULL;
	}
	crystal *newcrystal = (crystal *)malloc(sizeof(crystal));


	newcrystal->totalAtoms = crys1->totalAtoms+crys2->totalAtoms;
	newcrystal->numElements = crys1->numElements +crys2->numElements;
	crys_allocate(newcrystal);
	
	int iele1,iele2, ieleout,iatomout,iatom1, iatom2,idim;	
	for(iele1=0;iele1<crys1->numElements;iele1++)
	{
		sprintf(newcrystal->elements+iele1*namelength,"%s",crys1->elements+iele1*namelength);
		*(newcrystal->totalEachSpecies+iele1) = *(crys1->totalEachSpecies+iele1);
	}
	
	newcrystal->numElements = crys1->numElements;
	for(iele2=0;iele2<crys2->numElements;iele2++)
	{
		for(iele1=0;iele1<crys1->numElements && strcmp(crys1->elements+iele1*namelength,crys2->elements+iele2*namelength);iele1++);
		if(!(iele1<crys1->numElements))//Element does not match any element in crystal 1
		{
		sprintf(newcrystal->elements+namelength*newcrystal->numElements,"%s",crys2->elements+iele2*namelength);
		*(newcrystal->totalEachSpecies+newcrystal->numElements) = *(crys2->totalEachSpecies+iele2);
		(newcrystal->numElements)++;
		}
		else//Element is is both crystals
		{
			*(newcrystal->totalEachSpecies+iele1) += *(crys2->totalEachSpecies+iele2);
		}
	}
	
	for(ieleout=0,iatomout=0;ieleout<newcrystal->numElements;ieleout++)
	{
	for(iatom1=0;iatom1<crys1->totalAtoms;iatom1++)
	{
		if(!strcmp(crys1->species+iatom1*namelength,newcrystal->elements+ieleout*namelength))
		{
		sprintf(newcrystal->species+iatomout*namelength,"%s",crys1->species+iatom1*namelength);
			for(idim=0;idim<3;idim++)
			{
			*(newcrystal->positions+3*iatomout+idim) = *(crys1->positions+3*iatom1+idim);
			}
		iatomout++;
		}
	}
	
	for(iatom2=0;iatom2<crys2->totalAtoms;iatom2++)
	{
		if(!strcmp(crys2->species+iatom2*namelength,newcrystal->elements+ieleout*namelength))
		{
		sprintf(newcrystal->species+iatomout*namelength,"%s",crys2->species+iatom2*namelength);
			for(idim=0;idim<3;idim++)
			{
			*(newcrystal->positions+3*iatomout+idim) = *(crys2->positions+3*iatom2+idim) + *(offset+idim);
			}
		iatomout++;
		}
	}
	}
	for(idim=0;idim<9;idim++)
	{
		*(newcrystal->latticeVectors+idim) = *(crys1->latticeVectors+idim);
	}//Should be improved to constructively add vectors
	if(addvectors) //If this is not used the new structure and the original with have the same dimensions //only works for orthorhombic crystals
	{
		double epsilon = 0.1;
		int direction = 0;
	if(fabs(*(offset)) >epsilon)
		direction = 0;
		else if(fabs(*(offset+1)) >epsilon)
		direction = 3;
		else if(fabs(*(offset+2)) >epsilon)
		direction = 6;
	printf("Adding crystals along direction %d\n",direction);	
	for(idim=0;idim<3;idim++)
	if(addvectors ==1)
	*(newcrystal->latticeVectors+direction+idim) += fabs(*(offset+idim));
	else
	*(newcrystal->latticeVectors+direction+idim) += *(crys2->latticeVectors+direction+idim);
	}
	//printf("num atoms %d num elements %d\n",newcrystal->totalAtoms,newcrystal->numElements);
	return newcrystal;
}

//This function is messy and could use some cleanup, but it has not been used in awhile so I will postpone.
crystal * crys_linearInterpolation(crystal *crys1, crystal *crys2, int numimages)
{
	//If given two crystals with identical atoms but different positions and/or lattice vectors, this creates a series of crystals
	//which interpolates between the two original crystals.
	
	if(numimages == 0)
	return NULL;
	crystal *newcrystals = (crystal *)malloc(numimages*sizeof(crystal));
	int index,index2, crystalindex;
	double factor1, factor2; // What portion of the image is coming from crystal1/crystal2

	for(crystalindex=0;crystalindex<numimages;crystalindex++)
	{

	(newcrystals+crystalindex)->totalAtoms = crys1->totalAtoms;
	(newcrystals+crystalindex)->numElements = crys1->numElements;
	crys_allocate((newcrystals+crystalindex));
		
		factor1 = (double)(numimages-crystalindex)/(numimages+1);
		factor2 = (double)(crystalindex+1)/(numimages+1);
	for(index=0;index<crys1->numElements;index++)
	{
		sprintf((newcrystals+crystalindex)->elements+index*namelength,"%s",crys1->elements+index*namelength);
		*((newcrystals+crystalindex)->totalEachSpecies+index) = *(crys1->totalEachSpecies+index);
	}
	
	
	for(index=0;index<9;index++)
	{
		*((newcrystals+crystalindex)->latticeVectors+index) = factor1 * *(crys1->latticeVectors+index) + factor2 * *(crys2->latticeVectors+index);
	}//Should be improved to constructively add vectors
	
	double cutofflength = 3;
	double possiblevalue;
	double latticevector;
	for(index=0;index<crys1->totalAtoms;index++)
	{
		sprintf((newcrystals+crystalindex)->species+index*namelength,"%s",crys1->species+index*namelength);
				for(index2=0;index2<3;index2++)
				{
					if(abs(*(crys1->positions+3*index+index2)-*(crys2->positions+3*index+index2))<cutofflength)
				*((newcrystals+crystalindex)->positions+3*index+index2) = factor1 * *(crys1->positions+3*index+index2) + factor2 * *(crys2->positions+3*index+index2);
					else
					{
					latticevector = *((newcrystals+crystalindex)->latticeVectors+index2 + 3*index2); // Assumes lattice vectors are [{x,0,0}.{0,y,0},{0,0,z}]
					possiblevalue = factor1 * (*(crys1->positions+3*index+index2) + latticevector)  + factor2 * (*(crys2->positions+3*index+index2));
					if(possiblevalue > latticevector)
					possiblevalue = factor1 * (*(crys1->positions+3*index+index2))  + factor2 * (*(crys2->positions+3*index+index2) + latticevector);
					if(possiblevalue > latticevector)
					possiblevalue = factor1 * (*(crys1->positions+3*index+index2) - latticevector)  + factor2 * (*(crys2->positions+3*index+index2));
					if(possiblevalue<0)
					possiblevalue = factor1 * (*(crys1->positions+3*index+index2))  + factor2 * (*(crys2->positions+3*index+index2) - latticevector);
					if(possiblevalue<0)
					printf("no value found could be \n");
					*((newcrystals+crystalindex)->positions+3*index+index2) = possiblevalue;
					}
				}
	}
	}
	return newcrystals;
}

void crys_removeAtom(struct crystal *crys, int removalsite)
{	
	//printf("Removing atom from site %d\n",removalsite);
	if(crys == NULL)
	 {
		 printf("Crystal not initialized, cannot remove atom\n");
		 return;
	 }

	//char buffer[1000];
	//strcpy(buffer,crys->species+namelength*removalsite);
	int iele,iatom;
	for(iele=0;iele<crys->numElements;iele++)
	{
		if(!strcmp(crys->species+namelength*removalsite,crys->elements + namelength*iele))
		{
			//printf("The removed atom is %s\n",crys->species+namelength*removalsite);
			(*(crys->totalEachSpecies+iele))--;
			break;
		}
	}
	
	for(iatom = removalsite;iatom<crys->totalAtoms-1;iatom++)
	{
		*(crys->positions+3*iatom) = *(crys->positions + 3*(iatom+1));
		*(crys->positions+3*iatom+1) = *(crys->positions + 3*(iatom+1)+1);
		*(crys->positions+3*iatom+2) = *(crys->positions + 3*(iatom+1)+2);
		
		strcpy(crys->species + iatom*namelength,crys->species + (iatom+1)*namelength);
	}
	
	(crys->totalAtoms)--;
	//printf("There are now %d atoms\n",crys->totalAtoms);
}

void crys_removeAtoms(struct crystal *crys, int *atoms,int numRemove)
{	
	//When passed a ascending list of atoms to remove, this removes each one.
	//This is more efficient than calling crys_removeAtom n times because it only requires a single pass to compress the empty space.
	
	//printf("Removing atom from site %d\n",removalsite);
	if(crys == NULL)
	 {
		 printf("Crystal not initialized, cannot remove atom\n");
		 return;
	 }

	
	int iEle,iRemove,iAtom;
	for(iRemove=0;iRemove<numRemove;iRemove++)
		for(iEle=0;iEle<crys->numElements;iEle++)
		{
			if(!strcmp(crys->species+namelength* *(atoms+iRemove),crys->elements + namelength*iEle))
			{
				//printf("The removed atom is %s\n",crys->species+namelength*removalsite);
				(*(crys->totalEachSpecies+iEle))--;
				break;
			}
		}
	
	iRemove=1;//Counts how many atoms where removed with an index below the current index. This is how many positions to shift the data.
	int found;
	for(iAtom = *atoms;iAtom<crys->totalAtoms-numRemove;iAtom++)
	{
		found = 0;
		if(iRemove < numRemove && iAtom == (*(atoms+iRemove)-iRemove))
		{
			iRemove++;
			found = 1;
		}
		*(crys->positions+3*iAtom) = *(crys->positions + 3*(iAtom+iRemove));
		*(crys->positions+3*iAtom+1) = *(crys->positions + 3*(iAtom+iRemove)+1);
		*(crys->positions+3*iAtom+2) = *(crys->positions + 3*(iAtom+iRemove)+2);
		
		strcpy(crys->species + iAtom*namelength,crys->species + (iAtom+iRemove)*namelength);
		if(found)//If multiple consecutive atoms are to be removed, we need to recheck iAtom after each removal.
			iAtom--;
	}
	
	(crys->totalAtoms)-=numRemove;
	//printf("There are now %d atoms\n",crys->totalAtoms);
}

void crys_addAtom(struct crystal *crys, char *element, double x, double y, double z)
{	
	int iatom,iatom2;
	int iele;
	iele = crys_elementInString(crys->elements,crys->numElements,element);
	if(iele==-1)//Element is not in list Add element to end
	{
		iele = crys->numElements;
		iatom = crys->totalAtoms;
		(crys->numElements)++;
		*(crys->totalEachSpecies+iele) = 1;
		strcpy(crys->elements+iele*namelength,element);
	}
	else
	{
		iatom = crys_elementOffset(crys,element)+crys_elementCount(crys,element);
		//move all other atoms after this up by one
		for(iatom2 = crys->totalAtoms;iatom2>iatom;iatom2--)
		{
		*(crys->positions+3*(iatom2+1)) = *(crys->positions + 3*iatom2);
		*(crys->positions+3*(iatom2+1)+1) = *(crys->positions + 3*iatom2+1);
		*(crys->positions+3*(iatom2+1)+2) = *(crys->positions + 3*iatom2+2);
		
		strcpy(crys->species + (iatom2+1)*namelength,crys->species + iatom2*namelength);
		}
		
		(*(crys->totalEachSpecies+iele))++;
	}
	//strcpy((crys->elements+iele*namelength),element);
	strcpy((crys->species+iatom*namelength),element);
	*(crys->positions+3*iatom) = x;
	*(crys->positions+3*iatom+1) = y;
	*(crys->positions+3*iatom+2) = z;
	crys->totalAtoms++;
}

void crys_addAtoms(struct crystal *crys, char *element, int numAdd, double *x, double *y, double *z)
{	
	//Adds multiple atoms of a single element type at once.
	//This is faster than n calls to crys_addAtom because rearrangement to make room only requires a single pass.
	int iatom,iatom2;
	int iele;
	iele = crys_elementInString(crys->elements,crys->numElements,element);
	if(iele==-1)//Element is not in list Add element to end
	{
		iele = crys->numElements;
		iatom = crys->totalAtoms;
		(crys->numElements)++;
		*(crys->totalEachSpecies+iele) = numAdd;
		strcpy(crys->elements+iele*namelength,element);
	}
	else
	{
		iatom = crys_elementOffset(crys,element)+crys_elementCount(crys,element);
		//move all other atoms after this up by one
		for(iatom2 = crys->totalAtoms;iatom2>iatom;iatom2--)
		{
		*(crys->positions+3*(iatom2+numAdd)) = *(crys->positions + 3*iatom2);
		*(crys->positions+3*(iatom2+numAdd)+1) = *(crys->positions + 3*iatom2+1);
		*(crys->positions+3*(iatom2+numAdd)+2) = *(crys->positions + 3*iatom2+2);
		
		strcpy(crys->species + (iatom2+numAdd)*namelength,crys->species + iatom2*namelength);
		}
		
		(*(crys->totalEachSpecies+iele))+= numAdd;
	}
	//strcpy((crys->elements+iele*namelength),element);
	int iAdd;
	for(iAdd=0;iAdd<numAdd;iAdd++)
	{
	
	strcpy((crys->species+iatom*namelength),element);
	*(crys->positions+3*iatom) = *(x+iAdd);
	*(crys->positions+3*iatom+1) = *(y+iAdd);
	*(crys->positions+3*iatom+2) = *(z+iAdd);
	crys->totalAtoms++;
	iatom++;
	}
}

//This appears to be a duplicate of crys_elementCount. I prefer the method used in elementCount to this method
/*int crys_atomsOfElement(crystal *crys, char *element)
{
	//Returns the number of atoms of a specific element
	int iele;
	int atomsthiselement=-1;
	
	for(iele=0;iele<crys->numElements;iele++)
	{
	//printf("Comparing elements %s and %s\n",element,crys->elements+iele*namelength);
	if(!strcmp(element,crys->elements+iele*namelength))
	atomsthiselement=*(crys->totalEachSpecies+iele);
	}
	
	return atomsthiselement;
}*/

int crys_elementIndex(crystal *crys, char *element)
{
	//Returns the start index of a specific element
	int elementtype;
	
	for(elementtype=0;elementtype<crys->numElements;elementtype++)
	{
	//printf("Comparing elements %s and %s\n",element,crys->elements+elementtype*namelength);
	if(!strcmp(element,crys->elements+elementtype*namelength))
	//break;
	return elementtype;
	}
	
	printf("element %s not found in list\n",element);
	return -1;
}

int crys_elementOffset(crystal *crys, char *element)
{
	//Returns the start index of a specific element
	int elementtype;
	int crys_elementOffset=0;
	
	for(elementtype=0;elementtype<crys->numElements;elementtype++)
	{
	//printf("Comparing elements %s and %s\n",element,crys->elements+elementtype*namelength);
	if(!strcmp(element,crys->elements+elementtype*namelength))
	break;
	crys_elementOffset+=*(crys->totalEachSpecies+elementtype);
	}
	
	return crys_elementOffset;
}

void crys_elementBoundsArray(crystal *crys, int *offsetArray)
{
	int iEle;
	*offsetArray=0;
	for(iEle=0;iEle<crys->numElements;iEle++)
	{
		*(offsetArray+iEle+1) = *(offsetArray+iEle) + *(crys->totalEachSpecies+iEle);
	}
}

int crys_elementCount(crystal *crys, char *element)
{
	//Returns the number of atoms of a specific element
	int iele = crys_elementInString(crys->elements,crys->numElements,element);
	
	if(iele == -1)
	{
	//This may be a useful debugging message sometimes. However the return 0 is sometimes used on purpose
	//And this can make tremendously large logs.
	//printf("Element %s not found in list\n",element);
	return 0;
	}
	
	return *(crys->totalEachSpecies+iele);
}

//This function likely causes a memory leak if used repeatedly as the string is copied over to a new buffer and modified there. 
char * crys_appendElementString(char *elementString, int numElements,char *newEle)
{
	char *output = (char *)calloc(1000,sizeof(char));
	memcpy(output,elementString,1000*sizeof(char));
	sprintf(output + namelength*numElements,"%s",newEle);
	return output;
}

char * crys_elementString( int num, ... )
{
    va_list arguments;                     
    char *output = (char *)calloc(1000,sizeof(char));
	sprintf(output,"");
	
    va_start ( arguments, num );           
    int i;
    for ( i = 0; i < num; i++ )        
    {
        sprintf(output+namelength*i,"%s", va_arg ( arguments, char * )); 
    }
    va_end ( arguments );                  // Cleans up the list

    return output;                      
}


int crys_elementInString(char *elementList, int numElementsInList,char *element)
{
	//Returns the element's position in list or -1 if not found.
	int i;
	for(i=0;i<numElementsInList;i++)
	if(!strcmp(elementList+i*namelength,element))
	return i;
	
	return -1;
}

double crys_atomDistance(crystal *crys, int iatom1, int iatom2)
{
	//Finds shortest distance between two atoms in periodic boundary conditions of an orthorhomic crystal
	
	double x,y,z, lx,ly,lz,distance;
	lx = *(crys->latticeVectors);
	ly = *(crys->latticeVectors+4);
	lz = *(crys->latticeVectors+8);
	
	
	x = fabs(*(crys->positions+3*iatom1) - *(crys->positions+3*iatom2));
	y = fabs(*(crys->positions+3*iatom1+1) - *(crys->positions+3*iatom2+1));
	z = fabs(*(crys->positions+3*iatom1+2) - *(crys->positions+3*iatom2+2));
	
	if(x>lx/2)
	x-=lx;
	if(y>ly/2)
	y-=ly;
	if(z>lz/2)
	z-=lz;
	distance = sqrt(x*x+y*y+z*z);
	
	return distance;
	
}


double crys_atomVector(crystal *crys, int iatom1, int iatom2, double *r)
{
	//Finds shortest vector between two atoms with periodic boundary conditions in an orthorhomic crystal
	
	double lx,ly,lz,distance;
	lx = *(crys->latticeVectors);
	ly = *(crys->latticeVectors+4);
	lz = *(crys->latticeVectors+8);
	
	
	r[0] = fabs(*(crys->positions+3*iatom1) - *(crys->positions+3*iatom2));
	r[1] = fabs(*(crys->positions+3*iatom1+1) - *(crys->positions+3*iatom2+1));
	r[2] = fabs(*(crys->positions+3*iatom1+2) - *(crys->positions+3*iatom2+2));
	
	if(r[0]>lx/2)
	r[0]-=lx;
	if(r[1]>ly/2)
	r[1]-=ly;
	if(r[2]>lz/2)
	r[2]-=lz;
	distance = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
	
	return distance;
	
}

int crys_atomDirection(crystal *crys, int iatom1, int iatom2)
{
	double displacement[3];
	displacement[0] = *(crys->positions+3*iatom1)-*(crys->positions+3*iatom2);
	displacement[1] = *(crys->positions+3*iatom1+1)-*(crys->positions+3*iatom2+1);
	displacement[2] = *(crys->positions+3*iatom1+2)-*(crys->positions+3*iatom2+2);
	double epsilon = 1e-1;
	if(fabs(displacement[0])<epsilon)
	return 2;
	if(fabs(displacement[1])<epsilon)
	return 1;
	if(fabs(displacement[2])<epsilon)
	return 0;
	printf("The atom direction could not be found\n");
	printf("Displacement %g %g %g\n",displacement[0],displacement[1],displacement[2]);
	return 0;
}

void crys_printAtomPosition(crystal *crys, int iatom)
{
	printf("Atom %d Element %s (%g, %g, %g)\n",iatom,crys->species+iatom*namelength,*(crys->positions+3*iatom),*(crys->positions+3*iatom+1),*(crys->positions+3*iatom+2));
}

void crys_printAllAtoms(crystal *crys)
{
	int iatom;
	for(iatom=0;iatom<crys->totalAtoms;iatom++)
	crys_printAtomPosition(crys,iatom);
}

void crys_zValues(crystal *crys, double *list, int *numitems)
{
	*numitems = 0;
	int iatom,iz;
	double z;
	double epsilon = 1e-3;
	for(iatom=0;iatom<crys->totalAtoms;iatom++)
	{
		z = *(crys->positions+3*iatom+2);
		iz = binarysearchDouble(list,z,epsilon,*numitems);
		if(iz>=0)
		{
		listinsertdouble(list,*numitems,iz,z-epsilon/2);
		(*numitems)++;
		}
	}
	//Add 2 extra points at either end so the entire charged area is encapsulated.
	///padListLeft(list,numitems);
	//padListLeft(list,numitems);
	//padListRight(list,numitems);
	//padListRight(list,numitems);
	
	
}

char * crys_formatElementString(char *elementString, int numElements)
{
	char *buffer = (char *)calloc(1000,sizeof(char));
	strcpy(buffer,"");
	char buffer2[1000];
	//sprintf(buffer,"There are %d elements:",numElements);
	int iele;
	for(iele=0;iele<numElements;iele++)
	{
	sprintf(buffer2,"%s%s\t",buffer,elementString+iele*namelength);
	strcpy(buffer,buffer2);
	}
	return buffer;
}

//The crystal is now saved with 2 decimal points to save disk space
void crys_makexyz( crystal *crys, char *name)
{
	//FILE *outfile = fopen(strcat2(dir,name),"w");
	FILE *outfile = fopen(name,"w");
    if(outfile == NULL)
    printf("Unable to open file %s for writing\n",name);

	int imol,molcount,substitutionCount=0;
	for(imol=0;imol<numMolecules;imol++)
	{
		molcount = crys_elementCount(crys,molecules+imol*namelength);
		substitutionCount += molcount * ((molecularSubstitutions+imol)->totalAtoms-1);
		}
	
	fprintf(outfile,"%d\n",crys->totalAtoms+substitutionCount);//print number of atoms
	
	
	int ilv,iatom,iatommol;
	for(ilv = 0;ilv<9;ilv++)
	fprintf(outfile,"%.16f\t",*(crys->latticeVectors+ilv));//print lattice vectors
	fprintf(outfile,"\n");//add new line after lattice parameters
	for(iatom = 0;iatom<crys->totalAtoms;iatom++)
	//Substitute in the molecule during printing
	if((imol = crys_elementInString(molecules,numMolecules,(crys->species+namelength*iatom))) != -1)
	{
		//printf("%s is in the list of molecules\n",(crys->species+namelength*iatom));
		for(iatommol=0;iatommol<(molecularSubstitutions+imol)->totalAtoms;iatommol++)
		{
		//printf("Placing %s in system\n",((molecularSubstitutions+imol)->species+namelength*iatommol));
		fprintf(outfile,"%s\t%.2f\t%.2f\t%.2f\t\n",((molecularSubstitutions+imol)->species+namelength*iatommol),
		*(crys->positions+3*iatom)+*((molecularSubstitutions+imol)->positions+3*iatommol),
		*(crys->positions+3*iatom+1)+*((molecularSubstitutions+imol)->positions+3*iatommol+1),
		*(crys->positions+3*iatom+2)+*((molecularSubstitutions+imol)->positions+3*iatommol+2));//print positions
		}
	}
	else
	fprintf(outfile,"%s\t%.2f\t%.2f\t%.2f\t\n",(crys->species+namelength*iatom),*(crys->positions+3*iatom),*(crys->positions+3*iatom+1),*(crys->positions+3*iatom+2));//print positions
	fclose(outfile);
}


//If the file was saved with molecular substitution there is currently no way to recover the original.
//Lets add the option of saving/reading a crystal network at some point
crystal *crys_readxyz(char *name)
{
	FILE *infile = fopen(name,"r");
    if(infile == NULL)
    {
    printf("Unable to open file %s for reading\n",name);
	return NULL;
	}	

	crystal *crys = (crystal *)malloc(sizeof(crystal));
	char elements[maxelements*namelength];
	
	char element[namelength];
	int iele;
	char *line=NULL;
	size_t len=0;
	
	//First pass will be to count elements
	//Fist line contains number of atoms
	getline(&line,&len,infile);
	sscanf(line,"%d",&(crys->totalAtoms));
	//Second line contains lattice vectors
	getline(&line,&len,infile);
	sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
	(crys->latticeVectors+0),(crys->latticeVectors+1),(crys->latticeVectors+2),
	(crys->latticeVectors+3),(crys->latticeVectors+4),(crys->latticeVectors+5),
	(crys->latticeVectors+6),(crys->latticeVectors+7),(crys->latticeVectors+8));
	
	double x,y,z;
	crys->numElements=0;
	while(getline(&line,&len,infile) != -1)
	{
		sscanf(line,"%s %lf %lf %lf",element,&x,&y,&z);
		//printf("%s, %g %g %g\n",element,x,y,z);
		iele = crys_elementInString(elements,crys->numElements,element);
		if(iele == -1)
		{
			strcpy(elements+crys->numElements*namelength,element);
			(crys->numElements)++;
			
		}
		else
		(crys->totalAtoms)++;
	}
	crys_allocate(crys);
	
	//Second pass will be to fill in atoms
	//These counters will be re-incremented as ions are read in
	for(iele=0;iele<crys->numElements;iele++)
	{
	*(crys->totalEachSpecies+iele)=0;
	//This seems a bit redundant but it prevents crys->elements being overwritten during allocation
	strcpy(crys->elements+iele*namelength,elements+iele*namelength);
	}
	crys->totalAtoms=0;
	
	rewind(infile);
	//First line contains number of atoms
	getline(&line,&len,infile);
	//Second line contains lattice vectors
	getline(&line,&len,infile);
	
	while(getline(&line,&len,infile) != -1)
	{
		sscanf(line,"%s %lf %lf %lf",element,&x,&y,&z);
		//printf("%s, %g %g %g\n",element,x,y,z);
		
		*(crys->positions+3*crys->totalAtoms) = x;
		*(crys->positions+3*crys->totalAtoms+1) = y;
		*(crys->positions+3*crys->totalAtoms+2) = z;
		iele = crys_elementInString(crys->elements,crys->numElements,element);
		strcpy(crys->species+crys->totalAtoms*namelength,element);
		
		(*(crys->totalEachSpecies+iele))++;
		(crys->totalAtoms)++;
	}
	
	fclose(infile);
	return crys;
}


//This function is currently only used in VASPUI.c. It will be left here to allow for future recombining.
void crys_makeVASP(struct crystal *crys, char *name, int cart)
{
	printf("Making POSCAR for %s\n",name);
	char *buffer = (char *)calloc(1000,sizeof(char));
	sprintf(buffer,"%s",name);
	//printf("Outfile will be %s\n",buffer);
	FILE *outfile = fopen(buffer,"w");
	if(outfile == NULL)
	printf("Outfile is null in makePOSCAR\n");
	
	
	char *basefile =
	"%s\n"
 "1.0          ! universal scaling parameters\n"
 "%s"
 "%s\n"
 "%s           ! number of atoms\n"
 "%s         \n"
 "%s";

	char vectormatrix[10000];
	makematrix(vectormatrix,crys->latticeVectors,3,3,"\t");
	
	//printf("%s\n",vectormatrix);
	char positionmatrix[100000];
	char *keyword;
	if(cart)
	{
	makematrix(positionmatrix,crys->positions,crys->totalAtoms,3,"\t");
	keyword = "cart ! positions in cartesian coordinates";
	}
	else
	{
		keyword = "direct ! positions in fractional coordinates";
		double *inverse = inversematrix(crys->latticeVectors,3);
		double *fractionalcoordinates = (double *)malloc(crys->totalAtoms*3*sizeof(double));
		matrixmultiplication(fractionalcoordinates,crys->positions,inverse,crys->totalAtoms,3,3);
		makematrix(positionmatrix,fractionalcoordinates,crys->totalAtoms,3,"\t");
	}
	
	//printf("%s\n",positionmatrix);
	
	char atomnames[10000];
	sprintf(atomnames,"");

	
	char atomnumbers[10000];
	sprintf(atomnumbers,"");
	int index;
	for(index = 0;index<crys->numElements;index++)
	{
		if(index == 0)
		sprintf(buffer,"");
		else
		sprintf(buffer,"%s ",atomnumbers);
		sprintf(atomnumbers,"%s%d",buffer,*(crys->totalEachSpecies+index));
		
		strcpy(buffer,atomnames);
		sprintf(atomnames,"%s %s",buffer,(crys->elements+index*namelength));
	}
	
	fprintf(outfile,basefile,name,vectormatrix,atomnames,atomnumbers,keyword,positionmatrix);
	fclose(outfile);
}

void jmol(char *filename)
{
	char *command = (char *)calloc(1000,sizeof(char));
	sprintf(command,"jmol %s &",filename);
	system(command);
	free(command);
}


void vmdRender(char *dir, int id)
{
	char *basefile = 
	"mol new traj/%04d.xyz\n"
	"rotate y by 180\n"
	"mol rep CPK\n"
	"mol addrep 0\n"
	"axes location off\n"
	"color Name H black\n"
	"render TachyonInternal traj/%04d.TGA\n"
	"exit\n";
	FILE *outfile = fopen(strcat2(dir,"vmdscript"),"w");
	fprintf(outfile,basefile,id,id);
	fclose(outfile);
	system("vmd -dispdev text -e vmdscript");
}


//This is only used by VASPUI.c
struct crystal * crys_readVASP(char *filename)
{
	struct crystal *crys = (crystal *)malloc(sizeof(struct crystal));
	
	char buffer[1000];
	char buffer2[1000];
	char garbage[1000];
	//sprintf(buffer,"%s/CONTCAR",name);
	FILE *infile = fopen(filename,"r");
	if(infile == NULL)
	{
		printf("crystal %s not found \n",filename);
		return NULL;
	}
	printf("Extracting Crystal %s\n",filename);
	int index;
	
	
	fgets(garbage,1000,infile);
	 fgets(garbage,1000,infile);
	 for(index=0;index<3;index++)
	 {
		 fscanf(infile,"%lf %lf %lf%[^\n]%[\n]",crys->latticeVectors+3*index,crys->latticeVectors+3*index+1,crys->latticeVectors+3*index+2,garbage,garbage);
	 }

	 fgets(garbage,1000,infile);
	 fgets(buffer,1000,infile);
	 fgets(buffer2,1000,infile);
	 
	 char species[50];
	 char numeachelement[50];
	 int numspecies;
	 int speciesindex;
	 for(index=0,speciesindex = 0,numspecies=0;*(buffer+index) != '\0';index++,speciesindex = 0)
	 {
		 while(*(buffer+index) != ' ' && *(buffer+index) != '\t' && *(buffer+index) != '\n' && *(buffer+index) != '\0')
		 {
			 *(species+namelength*numspecies+speciesindex) = *(buffer+index);
			 index++;
			 speciesindex++;
		}
			if(speciesindex != 0)
			{
			 *(species+namelength*numspecies+speciesindex) = '\0';
			 //printf("Found species %s with length %d\n",species + namelength*numspecies,speciesindex);
			 numspecies++;
		 }
	 }
	
	//printf("buffer 2 %s \n",buffer2);
	 for(index=0,speciesindex = 0,numspecies=0,crys->totalAtoms = 0;*(buffer2+index) != '\0';index++,speciesindex = 0)
	 {
		 //startingindex = index;
		 while(*(buffer2+index) != ' ' && *(buffer2+index) != '\t' && *(buffer2+index) != '\n' && *(buffer2+index) != '\0')
		 {
			 *(numeachelement+namelength*numspecies+speciesindex) = *(buffer2+index);
			 index++;
			 speciesindex++;
		}
		if(speciesindex != 0)
			{
			 *(numeachelement+namelength*numspecies+speciesindex) = '\0';
			 //printf("Found numeachelement %s  with length %d\n",numeachelement + namelength*numspecies,speciesindex);
			 crys->totalAtoms += atoi(numeachelement + namelength*numspecies);
			 numspecies++;
		 }
	 }
	crys->numElements = numspecies;
	crys_allocate(crys);
	 //printf("Total atoms %d\n",crys->totalAtoms);
	 for(index=0;index<numspecies;index++)
	 {
		 *(crys->totalEachSpecies+index) = atoi(numeachelement + namelength*index);
		 sprintf(crys->elements + namelength*index,species+namelength*index);
		// printf("Element %s has %d atoms\n",crys->elements + namelength*index,*(crys->totalEachSpecies+index));
	 }
	 
	 
	fgets(garbage,1000,infile);
	printf("The cell representation is:  %s\n",garbage);
	 int totalAtoms = 0;
	 int elementindex = 0;
	 for(index=0;index<crys->totalAtoms;index++)
	{
		fscanf(infile,"%lf %lf %lf%[^\n]%[\n]",crys->positions+3*index,crys->positions+3*index+1,crys->positions+3*index+2,buffer,buffer);
		if(!(index<*(crys->totalEachSpecies+elementindex)+totalAtoms))
		{
				totalAtoms += *(crys->totalEachSpecies+elementindex);
				elementindex++;
		}
		sprintf(crys->species+namelength*index,crys->elements+namelength*elementindex);
		//printf("Element %d is %s\n",index,crys->species+namelength*index);
	}
	//printf("Done reading positions\n");
	//printmatrix(crys->positions,crys->totalAtoms,3);
	
	fclose(infile);
	
	printf("Comparison String is %s\n",garbage);
	if(strstr(garbage,"D") || strstr(garbage,"d"))
	{
	matrixmultiplication(crys->positions, crys->positions, crys->latticeVectors, crys->totalAtoms, 3, 3);
	printf("Multiplying positions by lattice vectors because this in direct representation\n");
	}
	else
	printf("Positions read exactly as is because this is cartesian\n");
	//printf("Converting to direct format\n");
	//printmatrix(crys->positions,crys->totalAtoms,3);
	
	//printf("Done extracting crystal\n");
	return crys;
}

//The only use of this function is in vestigialcode
void crys_interatomicDistances(crystal *crys, double *intercrys_atomDistances)
{
	int i,j;
	for(i=0;i<crys->totalAtoms;i++)
	for(j=0;j<crys->totalAtoms;j++)
	*(intercrys_atomDistances + i*(crys->totalAtoms) + j) = crys_atomDistance(crys,i,j);
}

//Only used by VASPUI.c
double minsetdistances(double *intercrys_atomDistances, int numtoremove,int *atomnumbers,int totalAtoms)
{
	int i,j;
	double mindistance = 1e10;
	double distance;
	for(i=0;i<numtoremove-1;i++)
	for(j=i+1;j<numtoremove;j++)
	{
	if(*(atomnumbers+i) >= totalAtoms || *(atomnumbers+j) >= totalAtoms || *(atomnumbers+j)<0 || *(atomnumbers+i)<0)
	{
	//printf("Atom numbers of %d and %d for i,j %d, %d\n",*(atomnumbers+i),*(atomnumbers+j),i,j);
	continue;
	}	
	distance = *(intercrys_atomDistances + *(atomnumbers+i)*totalAtoms+*(atomnumbers+j));
	if(distance<mindistance)
	{
	mindistance = distance;
	//printf("Mindistance reduced to %g from atoms %d and %d\n",distance,*(atomnumbers+i),*(atomnumbers+j));
	}
	if(mindistance == 0)
	printf("Somehow distance is zero\n");
	}
	//printf("Min distance of set is %g\n",mindistance);
	return mindistance;
}

//Only used by VASPUI.c
crystal * crys_maxDefectSpacing(crystal *inputcrys, char **atoms,int numremove)
{
	int i,j;
	printf("Searching for removal of %d atoms\n",numremove);
	for(i=0;i<numremove;i++)
	printf("%s\t",*(atoms+i));
	printf("\n");
	
	int numElements = inputcrys->numElements;
	int numtoremove[numElements];
	int crys_atomsOfElement[numElements];
	int crys_elementOffsets[numElements];
	int maxcombinations[numElements];
	int totalAtoms = inputcrys->totalAtoms;
	
	for(i=0;i<numElements;i++)
	{
		numtoremove[i] = 0;
		crys_atomsOfElement[i] = *(inputcrys->totalEachSpecies + i);
		crys_elementOffsets[i] = crys_elementOffset(inputcrys,inputcrys->elements + namelength*i);
	for(j=0;j<numremove;j++)
	if(!strcmp(inputcrys->elements + namelength*i,*(atoms+j)))
	{
		numtoremove[i]++;
	}
	}
	
	double intercrys_atomDistances[totalAtoms*totalAtoms];
	crys_interatomicDistances(inputcrys,intercrys_atomDistances);
	
	int iterations[numElements];
	int *combinations[numElements];
	int totalcombinations = 1;

	for(i=0;i<numElements;i++)
	{
	maxcombinations[i] = choose(crys_atomsOfElement[i],numtoremove[i]);	
	totalcombinations *= maxcombinations[i];
	iterations[i] = 0;
	combinations[i] = (int *)malloc(numtoremove[i]*sizeof(int));
	combination(combinations[i],crys_atomsOfElement[i],numtoremove[i],0);
	printf("There are %d combinations for removing %d of %d atoms of element %s\n",maxcombinations[i],numtoremove[i],crys_atomsOfElement[i],inputcrys->elements+i*namelength);
	}
	printf("Searching among %d combinations\n",totalcombinations);
	
	double minmaxdistance = 0;
	double newdistance;
	int atomnumbers[numremove];
	int removenumbers[numremove];
	int atomnumberfill;
	while(iterations[numElements-1]<maxcombinations[numElements-1])
	{
		atomnumberfill = 0;
		//printf("combination, atom numbers are:\t");
		for(i=0;i<numElements;i++)
		for(j=0;j<numtoremove[i];j++)
		{
		atomnumbers[atomnumberfill] = *(combinations[i]+j) + crys_elementOffsets[i] -1;//The minus 1 is there because the author of combination made it 1 based instead of 0 based
		//atomnumbers[atomnumberfill] = combinations[i][j] + crys_elementOffset[i];
		//printf("%d,%d\t",*(*(combinations+i)+j),atomnumbers[atomnumberfill]);
		atomnumberfill++;
		}
		//printf("\n");
		
		
		newdistance = minsetdistances(intercrys_atomDistances, numremove, atomnumbers,totalAtoms);
		
		if(newdistance > minmaxdistance)
		{
			//store the combinations
			minmaxdistance = newdistance;
			for(i=0;i<numremove;i++)
			removenumbers[i] = atomnumbers[i];
		}
		iterations[0]++;
		combination(combinations[0],crys_atomsOfElement[0],numtoremove[0],iterations[0]);
		i = 0;
		while(iterations[i] >= maxcombinations[i] && iterations[numElements-1]<maxcombinations[numElements-1])
		{
			iterations[i] = 0;
			iterations[i+1]++;
			combination(combinations[i],crys_atomsOfElement[i],numtoremove[i],iterations[i]);
			if(iterations[i+1] < maxcombinations[i+1])
			combination(combinations[i+1],crys_atomsOfElement[i+1],numtoremove[i+1],iterations[i+1]);
			i++;
		}
		//for(i=0;i<numElements;i++)
		//printf("Iteration %d of %d for element %s\n",iterations[i],maxcombinations[i],inputcrys->elements+i*namelength);
	}
	printf("The maximum distance is %g\n",minmaxdistance);
	
	crystal *crys = crys_multiply(inputcrys,1,1,1);
	for(i=numremove-1;i>=0;i--)//This loop is ran in reverse because if you start with the highest atom numbers it wont run into problems.
	{
		printf("Performing action on atom %d %s\n", removenumbers[i],crys->species+removenumbers[i]*namelength);
		//atomaction(crys,removenumbers[i],data);
		crys_removeAtom(crys,removenumbers[i]);//If this is methylammonium and it is not the very last thing to be removed were gonna have a bad time

	}
	return crys;
}

void cn_free(crystalnetwork *network)
{
	free(network->adjacencyList);
	free(network->numAdjacent);
	//free(network->vacancyList);
	free(network);
}

crystalnetwork *crys_duplicateNetwork(crystal *crys)
{

	crystalnetwork *newcn = (crystalnetwork *)malloc(sizeof(crystalnetwork));
	//newcn->numVacancies = crys->network->numVacancies;
	newcn->maxconnections = crys->network->maxconnections;
	newcn->adjacencyList = (int *)malloc(newcn->maxconnections*crys->totalAtoms*sizeof(int));
	newcn->numAdjacent = (int *)malloc(crys->totalAtoms*sizeof(int));
	//newcn->vacancyList = malloc(crys->network->numVacancies*sizeof(int));
	
	int i,j;
	for(i=0;i<crys->totalAtoms;i++)
	{
	*(newcn->numAdjacent+i) = *(crys->network->numAdjacent+i);
	for(j=0;j<*(newcn->numAdjacent+i);j++)
	*(newcn->adjacencyList+i*newcn->maxconnections+j) = *(crys->network->adjacencyList+i*newcn->maxconnections+j);
	}
		
	//for(i=0;i<crys->network->numVacancies;i++)
	//*(newcn->vacancyList+i) = *(crys->network->vacancyList+i);
	
	return newcn;
}

void crys_replaceRandomAtoms(crystal *crys, int numReplace, char *replaceElements, int numReplaceElements, char *newElement)
{
    printf("Replacing %d atoms from amongst %s with %s\n",numReplace,crys_formatElementString(replaceElements,numReplaceElements),newElement);
	//crys_printElements(crys);
	//printf("%s\n",crys_formatElementString(replaceElements,numReplaceElements));
	
	int ireplace,testAtom,removeAtoms[numReplace];
	double x[numReplace],y[numReplace],z[numReplace];
	for(ireplace=0;ireplace<numReplace;ireplace++)
	{
		testAtom = (crys->totalAtoms)*generaterandom();//The atom to change to vacancy
		while((crys_elementInString(replaceElements,numReplaceElements,crys->species+testAtom*namelength) == -1) || ((insertionSort(removeAtoms,ireplace,testAtom)) == -1))
		testAtom = (crys->totalAtoms)*generaterandom();
		//printf("Changing atom %d which is element %s to a vacancy\n",j,crys->species+j*namelength);
		
		//x[ireplace]=*(crys->positions+3*removeAtoms[ireplace]);
		//y[ireplace]=*(crys->positions+3*removeAtoms[ireplace]+1);
		//z[ireplace]=*(crys->positions+3*removeAtoms[ireplace]+2);
		//crys_removeAtom(crys,iatom);
		//crys_addAtom(crys,newElement,x,y,z);
	}
	for(ireplace=0;ireplace<numReplace;ireplace++)
	{
		x[ireplace]=*(crys->positions+3*removeAtoms[ireplace]);
		y[ireplace]=*(crys->positions+3*removeAtoms[ireplace]+1);
		z[ireplace]=*(crys->positions+3*removeAtoms[ireplace]+2);
	}
	//printlist(removeAtoms,numReplace);
	//crys_printAllAtoms(crys);
	crys_removeAtoms(crys,removeAtoms,numReplace);
	//crys_printAllAtoms(crys);
	crys_addAtoms(crys,newElement,numReplace,x,y,z);
	//crys_printAllAtoms(crys);
}

void crys_printElements(crystal *crys)
{
	printf("%s\n",crys_formatElementString(crys->elements,crys->numElements));
	int iele;
	for(iele=0;iele<crys->numElements;iele++)
	printf("%d\t",*(crys->totalEachSpecies+iele));
	printf("\n");
}

void crys_setupMolecules(char *newMolecules, int newNumMolecules, crystal *newMolecularSubstitutions)
{
	molecules=newMolecules;
	numMolecules=newNumMolecules;
	molecularSubstitutions = newMolecularSubstitutions;
}

crystal *crys_simpleCubic(char *element, double latticeConstant)
{
	crystal *outcrys = (crystal *)malloc(sizeof(crystal));
	outcrys->numElements=1;
	outcrys->totalAtoms=1;
	crys_allocate(outcrys);
	*(outcrys->totalEachSpecies)=1;
	
	double *lv = outcrys->latticeVectors;
	*(lv+0) = latticeConstant, *(lv+1) = 0,*(lv+2) = 0;
	*(lv+3) = 0,*(lv+4) = latticeConstant,*(lv+5) = 0;
	*(lv+6) = 0, *(lv+7) = 0, *(lv+8)=latticeConstant;
	
	sprintf(outcrys->elements,"%s",element);
	sprintf(outcrys->species,"%s",element);
	double *bp = outcrys->positions;
	*(bp+0) = 0,*(bp+1) = 0,*(bp+2) = 0;
	return outcrys;
}
