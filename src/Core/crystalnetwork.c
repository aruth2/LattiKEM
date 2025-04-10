/* crystalnetwork.h - Created By Anthony Ruth on June 17, 2017
 * Provides functionality for treating a crystal as a network of bonds and atoms.
 * 
 * */

#include "crystalnetwork.h"

char *vacancyatom = "V";

void cn_allocateCrystalNetwork(crystal *crys)
{
	printf("total Atoms %d, max connections %d\n",crys->totalAtoms,crys->network->maxconnections);
	crys->network->adjacencyList = (int *)calloc(crys->network->maxconnections*crys->totalAtoms,sizeof(int));
	crys->network->numAdjacent = (int *)calloc(crys->totalAtoms,sizeof(int));
	//int i;
	//for(i=0;i<crys->totalAtoms;i++)
	//*(crys->network->numAdjacent+i) = 0;
	
}

void cn_allocatedSize(crystal *crys)
{
	printf("Expected Size of crystal network is %d\n",crys->network->maxconnections*crys->totalAtoms*sizeof(int)+crys->totalAtoms*sizeof(int));
}

void cn_addLink(crystalnetwork *cn, int i, int j)
{
	//int *adjacencyList = crys->network->adjacencyList;
	//int *numAdjacent = crys->network->numAdjacent;
	
	int *adjacencyList = cn->adjacencyList;
	int *numAdjacent = cn->numAdjacent;
	
	*(adjacencyList+i*cn->maxconnections + *(numAdjacent+i)) = j;
    (*(numAdjacent+i))++;
    *(adjacencyList+j*cn->maxconnections + *(numAdjacent+j)) = i;
    (*(numAdjacent+j))++;
}

int cn_areConnected(crystal *crys, int i, int j)
{
	int *adjacencyList = crys->network->adjacencyList;
	int *numAdjacent = crys->network->numAdjacent;

	int iNeighbor;
	//printf("Checking if %d and %d are connected using maxconnections %d\n",i,j,crys->network->maxconnections);
	for(iNeighbor = 0;iNeighbor<*(numAdjacent+i);iNeighbor++)
		if (*(adjacencyList+i*crys->network->maxconnections + iNeighbor) == j)
			return 1;
	
	return 0;
}

void cn_bucketFillNetwork(crystal *crys, double nndistance, char *elementList, int numele, int bucketDirection)
{
	/*
	 * This function uses a bucket method to fill the adjacency list of the crystal.
	 * Only atoms which are within the bucket are compared.
	 * Therefore the scaling is nbucket * bucketsize^2
	 * Filling the buckets costs nbucket * totalAtoms
	 * The buckets are twice the width of the nndistance, 
	 * and the buckets are evenly spaced by the nndistance.
	 * A bucket exists which spans around the periodic boundary conditions.
	 * */
	//This assumes an orthorhombic cell with the lattice vectors oriented along the x, y, and z directions respectively.
	double maxExtent = *(crys->latticeVectors+ 3*bucketDirection + bucketDirection);
	//printf("maxExtent %g\n",maxExtent);
	int numBuckets = ceil(maxExtent/nndistance);
	int totalAtoms = crys->totalAtoms;
	//printf("Allocating %d buckets with size %d each\n",numBuckets,totalAtoms);
	int * buckets = (int *)calloc(numBuckets * totalAtoms,sizeof(int));
	int bucketSize[numBuckets];
	double bucketStart[numBuckets];
	int iBucket;
	int iAtom1, iAtom2;
	int iAtomBucket1, iAtomBucket2;
	double position;
	double bucketEnd;
	double distance;
	
	//cn_printAdjacencyList(crys);
	//Create Buckets
	for(iBucket = 0 ; iBucket < numBuckets ; iBucket++)
	{
		bucketStart[iBucket] = iBucket * nndistance;
		bucketSize[iBucket] = 0;
	}
	
	//Fill Buckets
	for(iAtom1 = 0 ; iAtom1 < crys->totalAtoms ; iAtom1++)
	{
		if(crys_elementInString(elementList,numele,crys->species+iAtom1*namelength)==-1)
		continue;
		position = *(crys->positions + 3 * iAtom1 + bucketDirection);
		
		for(iBucket = 0 ; iBucket < numBuckets ; iBucket++)
		{
			bucketEnd = bucketStart[iBucket] + 2 * nndistance;
			
			if( (position >= bucketStart[iBucket] && position <= bucketEnd) 
			|| ((bucketEnd >= maxExtent) 
			&& (position >= bucketStart[iBucket]) //Position is high and is in wrapping bucket
			|| (position <= bucketEnd-maxExtent))) //Position is low and is in wrapping bucket
			{
				buckets[totalAtoms*iBucket + bucketSize[iBucket]] = iAtom1;
				bucketSize[iBucket]++;
			} 
		}
	}
	
	//Compare Atoms within buckets
	for(iBucket = 0; iBucket < numBuckets; iBucket++)
	{
		for(iAtomBucket1 = 0 ; iAtomBucket1 < bucketSize[iBucket] - 1 ; iAtomBucket1++)
		{	
			iAtom1 = buckets[totalAtoms*iBucket + iAtomBucket1];
			
			for(iAtomBucket2 = iAtomBucket1+1 ; iAtomBucket2 < bucketSize[iBucket] ; iAtomBucket2++)
			{
				iAtom2 = buckets[totalAtoms*iBucket + iAtomBucket2];
				distance = crys_atomDistance(crys,iAtom1,iAtom2);
				if(distance<nndistance && !cn_areConnected(crys,iAtom1,iAtom2))
					cn_addLink(crys->network,iAtom1,iAtom2);
			}
		}
	}
	free(buckets);

}

void cn_fillNetwork(crystal *crys, double nndistance, char *elementList, int numele)
{
	//printf("Filling adjacency list\n");
		

	double distance;
	int i,j;

	//The number of attempted pairs could be reduced by using element pairs
	for(i=0;i<crys->totalAtoms;i++)
	{
		if(crys_elementInString(elementList,numele,crys->species+i*namelength)==-1)
		continue;
        for(j=i+1;j<crys->totalAtoms;j++)
        {
		if(crys_elementInString(elementList,numele,crys->species+j*namelength)==-1)
		continue;
		
        distance = crys_atomDistance(crys,i,j);
		
            if(distance<nndistance)
            {
                cn_addLink(crys->network,i,j);
                //printf("Found distance from atoms %g %g %g and %g %g %g to be %g\n",*(crys->positions+3*i),*(crys->positions+3*i+1),*(crys->positions+3*i+2),*(crys->positions+3*j),*(crys->positions+3*j+1),*(crys->positions+3*j+2),distance);
            }
        }
	}
		

}

void cn_fillFromnnd(crystal *crys, NearestNeighborDescriptor *nnd)
{
	//cn_fillNetwork(crys,nnd->nndistance,nnd->elementList,nnd->numEle);
	cn_bucketFillNetwork(crys,nnd->nndistance,nnd->elementList,nnd->numEle,DIRECTION_Z);
}

void cn_printnnd(NearestNeighborDescriptor *nnd)
{
	printf("NNdistance %g\n",nnd->nndistance);
	printf("%s\n",crys_formatElementString(nnd->elementList,nnd->numEle));
}
void cn_printAdjacencyList(crystal *crys,crystalnetwork *cn)
{
	int i,j;
	for(i=0;i<crys->totalAtoms;i++)
	{
        printf("\nAtom %d %s %d neighbors\n",i,crys->species+i*namelength,*(cn->numAdjacent+i));
        for(j=0;j<*(cn->numAdjacent+i);j++)
            printf("%d %s\t",*(cn->adjacencyList+i*cn->maxconnections+j),crys->species+*(cn->adjacencyList+i*cn->maxconnections+j) *namelength);
	}
	printf("\n");
}

void cn_nearestNeighbors(crystalnetwork *cn, int *neighbors, int *numneighbors, int index)
{
	//This function should be updated to allow restricting which atoms are considered.
	int *numAdjacent = cn->numAdjacent;
	
	*(numneighbors) = *(cn->numAdjacent+index);
	int i;
	for(i=0;i<*(numAdjacent+index);i++)
        *(neighbors+i) = *(cn->adjacencyList+index*cn->maxconnections+i);
	
}

void cn_nthNearestNeighbors(crystal *crys, int *nneighbors, int *numnneighbors, int index, int n)
{
	//This function finds the all the nth nearest neighbors of an atom.
	//The nth nearest neighbors exclude atoms which were neighbors of any previous shell.
	int usedatoms[crys->totalAtoms];
	int numusedatoms = 1;
	usedatoms[0] = index;
	
	int neighbors[crys->network->maxconnections];
	int numneighbors;
	
	int *shellatoms = (int *)malloc(crys->totalAtoms*sizeof(int));
	int *lastshellatoms= (int *)malloc(crys->totalAtoms*sizeof(int));
	lastshellatoms[0] = index;
	int numshellatoms;
	int numlastshellatoms = 1;
	int *shellswap;
	
	int shellnumber;
	int i,j;
	int insertionspot;
	for(shellnumber = 1;shellnumber<=n;shellnumber++)
	{
		//printf("In shell %d\n",shellnumber);
		numshellatoms=0;
		for(i=0;i<numlastshellatoms;i++)//Identify all atoms that are neighbors of atoms in last shell
		{
			//printf("Getting the nearest neighbors of atom %d\n",*(lastshellatoms+i));
            cn_nearestNeighbors(crys->network,neighbors,&numneighbors,*(lastshellatoms+i));
            //printf("Found %d neighbors of atom %d\n",numneighbors,*(lastshellatoms+i));
            for(j=0;j<numneighbors;j++)
            {
                //if(!intcontains(usedatoms,*(neighbors+j),numusedatoms)) // Using insert sort instead now
                //Atom has never been seen before, therefore it is part of shell
                if((insertionspot = binarysearch(usedatoms,*(neighbors+j),numusedatoms)) >= 0)
                {
                    *(shellatoms+numshellatoms) = *(neighbors+j);
                    numshellatoms++;
                    listinsert(usedatoms,numusedatoms,insertionspot,*(neighbors+j));
                    numusedatoms++;
                }
			//printf("insertionspot is %d result of listcontains is %d\n",insertionspot,intcontains(usedatoms,*neighbors+j,numusedatoms));
            }
		}
		shellswap = shellatoms; //Swap pointers so that lastshellatoms is populated and shellatoms has space to work
		shellatoms = lastshellatoms;
		lastshellatoms = shellswap;
		
		numlastshellatoms=numshellatoms;
	}
	for(i=0;i<numshellatoms;i++)
	*(nneighbors+i) = *(lastshellatoms+i); //Use lastshell instead of shell because of the swap
	*numnneighbors = numshellatoms;
	
	free(shellatoms);
	free(lastshellatoms);
}

void cn_integratednthNearestNeighbors(crystal *crys, int *integratednneighbors, int *numintegratednneighbors, int index, int nmax)
{
	int nneighbors[crys->totalAtoms];
	int numnneighbors;
	int n;
	
	*numintegratednneighbors = 0;
	int i;
	for(n=1;n<=nmax;n++)
	{
		cn_nthNearestNeighbors(crys,nneighbors,&numnneighbors,index,n);
	//	printf("Found %d cn_nthNearestNeighbors in shell %d\n",numcn_nearestNeighborseighbors,n);
		for(i=0;i<numnneighbors;i++)
		{
			*(integratednneighbors+*numintegratednneighbors) = *(nneighbors+i);
			(*numintegratednneighbors)++;
		}
	}
}


//This function is only used to determine the amount of data to allocate;
void numberOfCoordNeighbors(crystal *crys, char *coordElement, char *countedElements, int numCountedElements)
{
    int estart = crys_elementOffset(crys,coordElement);
    int n;
    int *nneighbors = (int *)malloc(crys->totalAtoms*sizeof(int));
    int numnneighbors;
    int ineighbor;
    int ielement;
    int numcounted;
    //would be nice to have more general answer than up to n=10
    for(n=1;n<10;n++)
    {
        cn_integratednthNearestNeighbors(crys,nneighbors,&numnneighbors,estart,n);
        numcounted = 0;
		for(ineighbor=0;ineighbor<numnneighbors;ineighbor++)
		{
			for(ielement=0;ielement<numCountedElements;ielement++)
			if(!strcmp(crys->species+*(nneighbors+ineighbor)*namelength,countedElements+ielement*namelength))
			{
				numcounted++;
			}
		}
        printf("Shell %d has %d coordinated atoms\n",n,numcounted);
    }
}


void cn_shellComposition(crystal *crys, int index, double *shells, char *shellelements, int numshellelements, double *distances, int *numshells)
{
	//When given an integer which points to a specific atom, this function considers the nearest neighbors, next nearest neighbors, next next nearest neighbors, etc
	//And finds the % of each element in each shell.
	//printf("Shell decomposition\n");
	
	int shellatoms[crys->totalAtoms];
	int shellnumber, crys_elementCountber, i;
	int numshellatoms=1; //Initialize to make it through loop first time
	double sum;
	for(shellnumber=1;numshellatoms>0;shellnumber++)
	{	
		printf("On shell number %d\n",shellnumber);
		cn_nthNearestNeighbors(crys,shellatoms,&numshellatoms,index,shellnumber);
		printf("Found %d atoms in this shell\n",numshellatoms);
		*(distances + shellnumber) = crys_atomDistance(crys,index,*(shellatoms));
		
		for(crys_elementCountber=0;crys_elementCountber<numshellelements;crys_elementCountber++)
		*(shells + shellnumber*numshellelements+crys_elementCountber) = 0;

		
		for(i=0;i<numshellatoms;i++)
		{
			//printf("Shellnumber %d i %d atom %d\n",shellnumber,i,*(shellatoms+i));
			for(crys_elementCountber=0;crys_elementCountber<numshellelements;crys_elementCountber++)
			if(!strcmp(crys->species+*(shellatoms+i) * namelength,shellelements+crys_elementCountber*namelength))
		 	//*(shells + shellnumber*numshellelements+crys_elementCountber) += 1.0/numshellatoms;
		 	*(shells + shellnumber*numshellelements+crys_elementCountber) += 1.0;

		}
		sum = 0;
		for(crys_elementCountber=0;crys_elementCountber<numshellelements;crys_elementCountber++)
		sum += *(shells + shellnumber*numshellelements+crys_elementCountber);
		
		printf("%g atoms in the shell were of the right type\n",sum);
		
		for(crys_elementCountber=0;crys_elementCountber<numshellelements;crys_elementCountber++)
		*(shells + shellnumber*numshellelements+crys_elementCountber) /= sum;
	}
	*numshells = shellnumber-1;
}

void cn_coordination(crystal *crys, char *coordElement, char *countedElements, int numCountedElements, int *coord, int coordinationNumber, int maxCoordination)
{
	//If passed a crystal, an element, and a series of other elements
	//This determines how many of each of the other elements surrounds the first element
	//The coordinationNumber can be used to find the coordination up to a further shell, for instance the next-nearest neighbors (coordinationNumber=2).
	
	int nneighbors[2*maxCoordination];//This is the wrong way to do this. Might not always work. A warning has been added below to identify this problem.
	
	int numnneighbors;
	
	int i, j,coordindex;
	
	int estart = crys_elementOffset(crys,coordElement);
	int numatoms = crys_elementCount(crys,coordElement);

	for(i=0;i<numatoms;i++)
	{
		//printf("On atom %d which is %s\n",i+estart,crys->species+(i+estart)*namelength);
		
		for(coordindex=0;coordindex<numCountedElements;coordindex++)
		*(coord+i*numCountedElements+coordindex) = 0;
		//printf("Before Species is located at %d\n",crys->species);

		if(coordinationNumber == 1)
		cn_nearestNeighbors(crys->network,nneighbors,&numnneighbors,i+estart);
		else
		cn_integratednthNearestNeighbors(crys,nneighbors,&numnneighbors,i+estart,coordinationNumber);
		
		//printf("After Species is located at %d\n",crys->species);
        if(numnneighbors > 2*maxCoordination)
        printf("In coordination, amount of space allocated for nearest neighbors has been exceeded %d allocated, %d used\n", 2*maxCoordination,numnneighbors);
		//printf("%d neighbors\n",numnneighbors);
		
		for(j=0;j<numnneighbors;j++)
		{
			//printf("On neighbor %d which is %d\n",j,*(nneighbors+j));
			
			//printf("On neighbor %d which is %d and atom %s\n",j,*(nneighbors+j),crys->species+*(nneighbors+j)*namelength);
			for(coordindex=0;coordindex<numCountedElements;coordindex++)
			if(!strcmp(crys->species+*(nneighbors+j)*namelength,countedElements+coordindex*namelength))
			{
				(*(coord+i*numCountedElements+coordindex))++;
			//printf("New coordination is %d\n",*(coord+i*numCountedElements+coordindex));
			}
		}
	}
}

void cn_bitA_integratednthNearestNeighbors(crystal *crys,crystalnetwork *cn, uint32_t * __restrict__ integratednneighbors, int *__restrict__ numintegratednneighbors, int * __restrict__ shellatoms, int * __restrict__ lastshellatoms, int index, int n)
{
	//printf("Location of pointers is %d, %d, %d, %d, %d\n",crys,integratednneighbors,numintegratednneighbors,shellatoms,lastshellatoms);
	//This function finds the all the nth nearest neighbors of an atom.
	//The nth nearest neighbors exclude atoms which were neighbors of any previous shell.
	int neighbors[crys->network->maxconnections];
	int numneighbors;
	int neighbor;

	int numshellatoms;
	int *shellswaptmp;
	
	int shellnumber;
	int i,j;

	SetBit(integratednneighbors,index);
	*numintegratednneighbors = 1;
	lastshellatoms[0] = index;
	int numlastshellatoms = 1;
	//printf("Finding neighbors up to %d shells\n",n);
	//printf("Allocating %d space for neighbors\n",maxconnections);

	for(shellnumber = 1;shellnumber<=n;shellnumber++)
	{
		//printf("In shell %d of %d\n",shellnumber,n);
		numshellatoms=0;
		for(i=0;(i<numlastshellatoms);i++)//Identify all atoms that are neighbors of atoms in last shell
		{
			//printf("Getting the nearest neighbors of atom %d\n",*(lastshellatoms+i));
			//printf("i %d of %d\n",i,numlastshellatoms);
            cn_nearestNeighbors(cn,neighbors,&numneighbors,*(lastshellatoms+i));
           // printf("Found %d neighbors of atom %d\n",numneighbors,*(lastshellatoms+i));
            for(j=0;j<numneighbors;j++)
            {
                //Atom has never been seen before, therefore it is part of shell
                neighbor = *(neighbors+j);
                //printf("Testing Neighbor %d\n",neighbor);
                if(!TestBit(integratednneighbors,neighbor))
                {
                    *(shellatoms+numshellatoms) = neighbor;
                    numshellatoms++;
                    SetBit(integratednneighbors,neighbor);                    
                    (*(numintegratednneighbors))++;
                    //printf("j %d integratednneighbors %d\n",j,*(numintegratednneighbors));
                }
            }
            //printf("i is %d of %d\n",i,numlastshellatoms);
		}
		
		shellswaptmp = shellatoms; //Swap pointers so that lastshellatoms is populated and shellatoms has space to work
		shellatoms = lastshellatoms;
		lastshellatoms = shellswaptmp;
		
		numlastshellatoms=numshellatoms;
	}
	
}

void cn_bitA_clusterApprox(crystal *__restrict__ crys, crystalnetwork *cn, char *__restrict__ coordElement, char *__restrict__ countedElements, int numCountedElements, int *__restrict__ coord, int coordinationNumber, int maxCoordination)
{
	//If passed a crystal, an element, and a series of other elements
	//This determines how many of each of the other elements surrounds the first element
	//The coordinationNumber can be used to find the coordination up to a further shell, for instance the next-nearest neighbors (coordinationNumber=2).
	//printf("Performing cluster approx for coordination number %d with maxCoordination %d\n",coordinationNumber,maxCoordination);

	int estart = crys_elementOffset(crys,coordElement);
	int numatoms = crys_elementCount(crys,coordElement);
	int totalAtoms = crys->totalAtoms;
	
	int totalElements = crys->numElements;
	int elementBoundsArray[totalElements+1];
	crys_elementBoundsArray(crys,elementBoundsArray);
	int elementIndicies[numCountedElements];
	
	//printf("Allocating arrays of size %d %d %d\n",numatoms*totalAtoms,numatoms*totalAtoms,numatoms*2*maxCoordination);
	int shellatoms[totalAtoms];
	int lastshellatoms[totalAtoms];
	uint32_t nneighbors[roundup(totalAtoms,32)/32];//Bit array of visited neighbors
	//printBits(nneighbors,sizeof(nneighbors));

	int i;
	int numnneighbors,coordindex;
	for(coordindex=0;coordindex<numCountedElements;coordindex++)
	{
	elementIndicies[coordindex] = crys_elementIndex(crys,countedElements+coordindex*namelength);
	for(i=0;i<numatoms;i++)
		*(coord+i*numCountedElements+coordindex) = 0;
	//printf("Element index for %s is %d\n",countedElements+namelength*coordindex,elementIndicies[coordindex]);
	}
	
	for(i=0;i<numatoms;i++)
	{
		memset(nneighbors,0,sizeof(nneighbors));
		//printf("On atom %d which is %s\n",i+estart,crys->species+(i+estart)*namelength);
		for(coordindex=0;coordindex<numCountedElements;coordindex++)
			*(coord+i*numCountedElements+coordindex) = 0;
		//printf("Before Species is located at %d\n",crys->species);
		cn_bitA_integratednthNearestNeighbors(crys,cn,nneighbors,&numnneighbors,shellatoms,lastshellatoms,i+estart,coordinationNumber);
		//printf("After Species is located at %d\n",crys->species);
		//printBits(nneighbors,sizeof(nneighbors));
		//printf("%d nneighbors %d\n",numnneighbors,&numnneighbors);
		for(coordindex = 0;coordindex<numCountedElements;coordindex++)
		{
			int iEle = elementIndicies[coordindex];
			int iAtom;
			int start=elementBoundsArray[iEle];
			int end=elementBoundsArray[iEle+1];
			int eleCount = 0;
			for(iAtom=start;iAtom<end;iAtom++)
				eleCount += TestBit(nneighbors,iAtom);
			*(coord+i*numCountedElements+coordindex) = eleCount;
		}

		//printf("Finished counting neighbors\n");
	}
	//printf("Finished Coordination\n");
}

void cn_collapseNetwork(crystal *__restrict__ crys, crystalnetwork * __restrict__ cn, char *__restrict__ coordElement, char *__restrict__ countedElements, int numCountedElements, int *__restrict__ coord, int coordinationNumber, int maxCoordination)
{
	//If passed a crystal, a coordinating element, an array of associated elements, and a number of hops,
	//This creates a network between the coordinating element 
	//and the associated elements up to n hops away along the original network
	
	int estart = crys_elementOffset(crys,coordElement);
	int numatoms = crys_elementCount(crys,coordElement);
	int totalAtoms = crys->totalAtoms;
	
	int totalElements = crys->numElements;
	int elementBoundsArray[totalElements+1];
	crys_elementBoundsArray(crys,elementBoundsArray);
	int elementIndicies[numCountedElements];
	
	//printf("Allocating arrays of size %d %d %d\n",numatoms*totalAtoms,numatoms*totalAtoms,numatoms*2*maxCoordination);
	int shellatoms[totalAtoms];
	int lastshellatoms[totalAtoms];
	uint32_t nneighbors[roundup(totalAtoms,32)/32];//Bit array of visited neighbors
	//printBits(nneighbors,sizeof(nneighbors));

	int i;
	int numnneighbors,coordindex;
	for(coordindex=0;coordindex<numCountedElements;coordindex++)
	{
	elementIndicies[coordindex] = crys_elementIndex(crys,countedElements+coordindex*namelength);
	//for(i=0;i<numatoms;i++)
	//	*(coord+i*numCountedElements+coordindex) = 0;
	//printf("Element index for %s is %d\n",countedElements+namelength*coordindex,elementIndicies[coordindex]);
	}
	
	for(i=0;i<numatoms;i++)
	{
		memset(nneighbors,0,sizeof(nneighbors));
		//printf("On atom %d which is %s\n",i+estart,crys->species+(i+estart)*namelength);
		//for(coordindex=0;coordindex<numCountedElements;coordindex++)
		//	*(coord+i*numCountedElements+coordindex) = 0;
		//printf("Before Species is located at %d\n",crys->species);
		cn_bitA_integratednthNearestNeighbors(crys,crys->network,nneighbors,&numnneighbors,shellatoms,lastshellatoms,i+estart,coordinationNumber);
		//printf("After Species is located at %d\n",crys->species);
		//printBits(nneighbors,sizeof(nneighbors));
		//printf("%d nneighbors %d\n",numnneighbors,&numnneighbors);
		for(coordindex = 0;coordindex<numCountedElements;coordindex++)
		{
			int iEle = elementIndicies[coordindex];
			int iAtom;
			int start=elementBoundsArray[iEle];
			int end=elementBoundsArray[iEle+1];
			//int eleCount = 0;
			
			for(iAtom=start;iAtom<end;iAtom++)
				if(TestBit(nneighbors,iAtom))
					{
					cn_addLink(cn,i+estart,iAtom);
				//	eleCount++;
					}
			//This does not work because there are 3 atoms in the network, I, Br, VX, but only 2 atoms counted for coordination.		
			//printf("Setting state %d coordindex %d to %d\n",i,coordindex,eleCount);		
			//*(coord+i*numCountedElements+coordindex) = eleCount;		
		}
	}
}
void cn_networkSwap(crystalnetwork *cn, int atom1, int atom2)
{
	int *adjacencyList = cn->adjacencyList;
	int *numAdjacent = cn->numAdjacent;
	int maxconnections = cn->maxconnections;	
		
	int atom1neighbors[maxconnections];
	int atom1ids[maxconnections];
	
	int numatom1neighbors;
	cn_nearestNeighbors(cn,atom1neighbors,&numatom1neighbors,atom1);
	
	int atom2neighbors[maxconnections];
	int atom2ids[maxconnections];
	int numatom2neighbors;
	cn_nearestNeighbors(cn,atom2neighbors,&numatom2neighbors,atom2);
	int i,j;
	//These two operations are quadratic in the number of neighbors. Looks relatively expensive for high neighbor count.
	for(i=0;i<numatom1neighbors;i++)//All neighbors of atom1 have their pointer to atom 1 identified
	{
		for(j=0;j<*(numAdjacent+*(atom1neighbors+i));j++)
		{
			if(*(adjacencyList+maxconnections* *(atom1neighbors+i)+j) == atom1)
			{
			atom1ids[i] = maxconnections* *(atom1neighbors+i)+j;
			break;
			}
		}
		//This check slows down the loop. Disable when not debugging.
		if(j==*(numAdjacent+*(atom1neighbors+i)))
		printf("Atom1 was not found in list\n");
	}
	for(i=0;i<numatom2neighbors;i++)//All neighbors of atom2 have their pointer to atom 2 identified
	{
		for(j=0;j<*(numAdjacent+*(atom2neighbors+i));j++)
		{
			if(*(adjacencyList+maxconnections* *(atom2neighbors+i)+j) == atom2)
			{
			atom2ids[i] = maxconnections* *(atom2neighbors+i)+j;
			break;
			}
		}
		//This check slows down the loop. Disable when not debugging.
		if(j==*(numAdjacent+*(atom2neighbors+i)))
		printf("Atom2 was not found in list\n");
	}
	for(i=0;i<numatom1neighbors;i++)//All neighbors of atom1 are changed to atom2
	*(adjacencyList+atom1ids[i])=atom2;
	for(i=0;i<numatom2neighbors;i++)//All neighbors of atom2 are changed to atom1
	*(adjacencyList+atom2ids[i])=atom1;

	int swap;
	for(i=0;i<maxconnections;i++)//The adjacency lists of atom1 and atom2 are swapped
	{
		swap = *(adjacencyList+atom1*maxconnections+i);
		 *(adjacencyList+atom1*maxconnections+i) = *(adjacencyList+atom2*maxconnections+i); 
		*(adjacencyList+atom2*maxconnections+i) = swap;
	}
		
		swap = *(numAdjacent+atom1);
		*(numAdjacent+atom1) = *(numAdjacent+atom2);
		*(numAdjacent+atom2) = swap;
}

void cn_swap(crystal *crys, int atom1, int atom2, int hasNetwork)
{
	//Swaps two atoms by interchanging their positions and their connections
	//If an atom and a vacancy are used, then this is a hop
	
	double distance = pow(*(crys->positions+3*atom1)-*(crys->positions+3*atom2),2)+ pow(*(crys->positions+3*atom1+1)-*(crys->positions+3*atom2+1),2) + pow(*(crys->positions+3*atom1+2)-*(crys->positions+3*atom2+2),2);
	distance = sqrt(distance);
	//printf("Swapping atoms %d and %d which are elements %s and %s distance %g\n",atom1,atom2,crys->species+namelength* atom1,crys->species+namelength* atom2,distance);
	double xswap,yswap,zswap;
	xswap = *(crys->positions+3*atom1);
	yswap = *(crys->positions+3*atom1+1);
	zswap = *(crys->positions+3*atom1+2);
	
	*(crys->positions+3*atom1) = *(crys->positions+3*atom2);
	*(crys->positions+3*atom1+1) = *(crys->positions+3*atom2+1);
	*(crys->positions+3*atom1+2) = *(crys->positions+3*atom2+2);
	
	*(crys->positions+3*atom2) = xswap;
	*(crys->positions+3*atom2+1) = yswap;
	*(crys->positions+3*atom2+2) = zswap;
	
	if(hasNetwork == CN_NO_NETWORK)
	return;
	
	cn_networkSwap(crys->network,atom1,atom2);
}
/*
void cn_swap(crystal *crys, int atom1, int atom2, int hasNetwork)
{
	//Swaps two atoms by interchanging their positions and their connections
	//If an atom and a vacancy are used, then this is a hop
	
	int *adjacencyList = crys->network->adjacencyList;
	int *numAdjacent = crys->network->numAdjacent;
	
	double distance = pow(*(crys->positions+3*atom1)-*(crys->positions+3*atom2),2)+ pow(*(crys->positions+3*atom1+1)-*(crys->positions+3*atom2+1),2) + pow(*(crys->positions+3*atom1+2)-*(crys->positions+3*atom2+2),2);
	distance = sqrt(distance);
	//printf("Swapping atoms %d and %d which are elements %s and %s distance %g\n",atom1,atom2,crys->species+namelength* atom1,crys->species+namelength* atom2,distance);
	double xswap,yswap,zswap;
	xswap = *(crys->positions+3*atom1);
	yswap = *(crys->positions+3*atom1+1);
	zswap = *(crys->positions+3*atom1+2);
	
	*(crys->positions+3*atom1) = *(crys->positions+3*atom2);
	*(crys->positions+3*atom1+1) = *(crys->positions+3*atom2+1);
	*(crys->positions+3*atom1+2) = *(crys->positions+3*atom2+2);
	
	*(crys->positions+3*atom2) = xswap;
	*(crys->positions+3*atom2+1) = yswap;
	*(crys->positions+3*atom2+2) = zswap;
	
	if(hasNetwork == CN_NO_NETWORK)
	return;
	
	int maxconnections = crys->network->maxconnections;	
		
	int atom1neighbors[maxconnections];
	int atom1ids[maxconnections];
	
	int numatom1neighbors;
	cn_nearestNeighbors(crys,atom1neighbors,&numatom1neighbors,atom1);
	
	int atom2neighbors[maxconnections];
	int atom2ids[maxconnections];
	int numatom2neighbors;
	cn_nearestNeighbors(crys,atom2neighbors,&numatom2neighbors,atom2);
	int i,j;
	for(i=0;i<numatom1neighbors;i++)//All neighbors of atom1 have their pointer to atom 1 identified
	{
		for(j=0;j<*(crys->network->numAdjacent+*(atom1neighbors+i));j++)
		{
			if(*(adjacencyList+maxconnections* *(atom1neighbors+i)+j) == atom1)
			{
			atom1ids[i] = maxconnections* *(atom1neighbors+i)+j;
			break;
			}
		}
		//This check slows down the loop. Disable when not debugging.
		if(j==*(crys->network->numAdjacent+*(atom1neighbors+i)))
		printf("Atom1 was not found in list\n");
	}
	for(i=0;i<numatom2neighbors;i++)//All neighbors of atom2 have their pointer to atom 2 identified
	{
		for(j=0;j<*(crys->network->numAdjacent+*(atom2neighbors+i));j++)
		{
			if(*(adjacencyList+maxconnections* *(atom2neighbors+i)+j) == atom2)
			{
			atom2ids[i] = maxconnections* *(atom2neighbors+i)+j;
			break;
			}
		}
		//This check slows down the loop. Disable when not debugging.
		if(j==*(crys->network->numAdjacent+*(atom2neighbors+i)))
		printf("Atom2 was not found in list\n");
	}
	for(i=0;i<numatom1neighbors;i++)//All neighbors of atom1 are changed to atom2
	*(adjacencyList+atom1ids[i])=atom2;
	for(i=0;i<numatom2neighbors;i++)//All neighbors of atom2 are changed to atom1
	*(adjacencyList+atom2ids[i])=atom1;

	
	
	int swap;
	for(i=0;i<maxconnections;i++)//The adjacency lists of atom1 and atom2 are swapped
	{
		swap = *(adjacencyList+atom1*maxconnections+i);
		 *(adjacencyList+atom1*maxconnections+i) = *(adjacencyList+atom2*maxconnections+i); 
		*(adjacencyList+atom2*maxconnections+i) = swap;
	}
		
		swap = *(numAdjacent+atom1);
		*(numAdjacent+atom1) = *(numAdjacent+atom2);
		*(numAdjacent+atom2) = swap;
}*/
