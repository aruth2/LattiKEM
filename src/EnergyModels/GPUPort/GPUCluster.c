#include "GPUCluster.h"
#include <nvToolsExt.h>


void ACC_nearestNeighbors(crystal *crys, int *neighbors, int *numneighbors, int index)
{
	//This function should be updated to allow restricting which atoms are considered.
	int *numAdjacent = crys->network->numAdjacent;
	
	*(numneighbors) = *(numAdjacent+index);
	int i;

#pragma acc loop seq
	for(i=0;i<*numneighbors;i++)
        *(neighbors+i) = *(crys->network->adjacencyList+index*maxconnections+i);
	
}

void ACC_integratednthNearestNeighbors(crystal *crys, uint32_t *restrict integratednneighbors, int *restrict numintegratednneighbors, int *restrict shellatoms, int *restrict lastshellatoms, int index, int n)
{
	//printf("Location of pointers is %d, %d, %d, %d, %d\n",crys,integratednneighbors,numintegratednneighbors,shellatoms,lastshellatoms);
	//This function finds the all the nth nearest neighbors of an atom.
	//The nth nearest neighbors exclude atoms which were neighbors of any previous shell.
	int neighbors[maxconnections];
	int numneighbors;

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
	int neighbor;
	
	for(shellnumber = 1;shellnumber<=n;shellnumber++)
	{
		//printf("In shell %d of %d\n",shellnumber,n);
		numshellatoms=0;
		for(i=0;;i++)//Identify all atoms that are neighbors of atoms in last shell
		//for(i=0;(i<numlastshellatoms);i++)//Identify all atoms that are neighbors of atoms in last shell
		{
			//printf("Getting the nearest neighbors of atom %d\n",*(lastshellatoms+i));
			//printf("i %d of %d\n",i,numlastshellatoms);
            ACC_nearestNeighbors(crys,neighbors,&numneighbors,*(lastshellatoms+i));
            //printf("Found %d neighbors of atom %d\n",numneighbors,*(lastshellatoms+i));
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
            if(i+1 == numlastshellatoms)//The normal loop condition of the for loop was not working correctly. So I added a variable that tracks whether or not to loop.
				break;
		}
		
		shellswaptmp = shellatoms; //Swap pointers so that lastshellatoms is populated and shellatoms has space to work
		shellatoms = lastshellatoms;
		lastshellatoms = shellswaptmp;
		
		numlastshellatoms=numshellatoms;
	}
	
}

int roundup(int number, int multiplicity)
{
	if(number % multiplicity == 0)
	return number;
	else
	return (number/multiplicity+1)*multiplicity;
}

void ACC_clusterApprox(int numGangs, int vectorLength, crystal *restrict crys, char *restrict coordElement, char *restrict countedElements, int numCountedElements, int *restrict coord, int coordinationNumber, int maxCoordination)
{
	//If passed a crystal, an element, and a series of other elements
	//This determines how many of each of the other elements surrounds the first element
	//The coordinationNumber can be used to find the coordination up to a further shell, for instance the next-nearest neighbors (coordinationNumber=2).

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
	
	uint32_t nneighbors[totalAtoms/32+1];
	
	int i;
	int numnneighbors,coordindex,iEle,eleCount;
	for(coordindex=0;coordindex<numCountedElements;coordindex++)
	{
	elementIndicies[coordindex] = crys_elementIndex(crys,countedElements+coordindex*namelength);
	for(i=0;i<numatoms;i++)
		*(coord+i*numCountedElements+coordindex) = 0;
	//printf("Element index for %s is %d\n",countedElements+namelength*coordindex,elementIndicies[coordindex]);
	}
	nvtxRangePushA("Cluster Approx");
	//for(iEle=0;iEle<totalElements+1;iEle++)
	//printf("Bound %d\n",elementBoundsArray[iEle]);
	//for(coordindex=0;coordindex<numCountedElements;coordindex++)
	//int eleCount;	
#pragma acc parallel loop vector num_gangs(numGangs) vector_length(vectorLength)\
private(coordindex,iEle,eleCount,numnneighbors)\
copyin(lastshellatoms[0:totalAtoms],shellatoms[0:totalAtoms],nneighbors[0:totalAtoms/32+1])\
private(lastshellatoms[0:totalAtoms],shellatoms[0:totalAtoms],nneighbors[0:totalAtoms/32+1])\
copyin(elementIndicies[0:numCountedElements],elementBoundsArray[0:totalElements+1])\
copyin(estart,maxCoordination,coordinationNumber,numCountedElements,maxconnections,totalAtoms)\
copyin(crys[0:1])\
copyin(crys->elements[0:namelength*totalElements],crys->positions[0:3*totalAtoms],crys->species[0:namelength*totalAtoms],crys->totalEachSpecies[0:totalElements],crys->network[0:1])\
copyin(crys->network->adjacencyList[0:maxconnections*totalAtoms],crys->network->numAdjacent[0:totalAtoms])\
copy(coord[0:numatoms*numCountedElements]) 
	for(i=0;i<numatoms;i++)
	{
		memset(nneighbors,0,sizeof(nneighbors));
		//printf("On atom %d which is %s\n",i+estart,crys->species+(i+estart)*namelength);
		//printf("%d %d %d\n",__pgi_gangidx(),__pgi_vectoridx(),i);
		//printf("Memory Locations %d %d %d %d\n",(uintptr_t)lastshellatoms,(uintptr_t)shellatoms,(uintptr_t)nneighbors,(uintptr_t)&numnneighbors);
		//printf("Before Species is located at %d\n",crys->species);
		ACC_integratednthNearestNeighbors(crys,nneighbors,&numnneighbors,shellatoms,lastshellatoms,i+estart,coordinationNumber);
		//printf("After Species is located at %d\n",crys->species);
		//printf("%d nneighbors\n",numnneighbors);
		for(coordindex = 0;coordindex<numCountedElements;coordindex++)
		{
			//*(coord+i*numCountedElements+coordindex) = 0;
			int iEle = elementIndicies[coordindex];
			int iAtom;
			int start=elementBoundsArray[iEle];
			int end=elementBoundsArray[iEle+1];
			int eleCount = 0;
			//#pragma acc loop seq
			//#pragma acc loop reduction(+:eleCount)
			for(iAtom=start;iAtom<end;iAtom++)
				eleCount += TestBit(nneighbors,iAtom);
			//printf("%d eleCount\n",eleCount);
			*(coord+i*numCountedElements+coordindex) = eleCount;
		}
		//printf("Finished counting neighbors\n");
	}
	//printf("Finished Coordination\n");
	nvtxRangePop();
}

void  SetBit( uint32_t A[],  int k )
{
    A[k/32] |= 1 << (k%32);  // Set the bit at the k-th position in A[i]
}

int TestBit( uint32_t A[],  int k )
{
    return ( (A[k/32] & (1 << (k%32) )) != 0 );     
}

// Assumes little endian
void printBits(void const * const ptr,size_t const size)
{
    unsigned char *b = (unsigned char*) ptr;
    unsigned char byte;
    int i, j;
    
    for (i = size-1; i >= 0; i--) {
        for (j = 7; j >= 0; j--) {
            byte = (b[i] >> j) & 1;
            printf("%u", byte);
        }
    }
    puts("");
}

void CPU_integratednthNearestNeighbors(crystal *crys, uint32_t *restrict integratednneighbors, int *restrict numintegratednneighbors, int *restrict shellatoms, int *restrict lastshellatoms, int index, int n)
{
	//printf("Location of pointers is %d, %d, %d, %d, %d\n",crys,integratednneighbors,numintegratednneighbors,shellatoms,lastshellatoms);
	//This function finds the all the nth nearest neighbors of an atom.
	//The nth nearest neighbors exclude atoms which were neighbors of any previous shell.
	int neighbors[maxconnections];
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
            cn_nearestNeighbors(crys,neighbors,&numneighbors,*(lastshellatoms+i));
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

void CPU_clusterApprox(crystal *restrict crys, char *restrict coordElement, char *restrict countedElements, int numCountedElements, int *restrict coord, int coordinationNumber, int maxCoordination)
{
	//If passed a crystal, an element, and a series of other elements
	//This determines how many of each of the other elements surrounds the first element
	//The coordinationNumber can be used to find the coordination up to a further shell, for instance the next-nearest neighbors (coordinationNumber=2).

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
		CPU_integratednthNearestNeighbors(crys,nneighbors,&numnneighbors,shellatoms,lastshellatoms,i+estart,coordinationNumber);
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

void CUDA_nearestNeighbors(crystal *crys, int *neighbors, int *numneighbors, int index)
{
	//This function should be updated to allow restricting which atoms are considered.
	int *numAdjacent = crys->network->numAdjacent;
	
	*(numneighbors) = *(numAdjacent+index);
	int i;

	for(i=0;i<*numneighbors;i++)
        *(neighbors+i) = *(crys->network->adjacencyList+index*maxconnections+i);
	
}


__global__ void CUDA_crys_printElements(crystal *crys)
{
	//printf("%s\n",crys_formatElementString(crys->elements,crys->numElements));
	int iele;
	for(iele=0;iele<crys->numElements;iele++)
	printf("%d\t",*(crys->totalEachSpecies+iele));
	printf("\n");
}

crystal * CUDA_crys_copyToDevice(crystal *crys)
{
	crystal * dev_crys;
	cudaMalloc((void **)&crys_device,sizeof(crystal));
	char * dev_species;
	cudaMalloc((void **)&(dev_species),crys->totalAtoms * namelength * sizeof(char));
	cudaMemcpy(&(dev_crys->species),dev_species,sizeof(char *));
	double *dev_positions;
	cudaMalloc((void **)&(dev_positions),3*crys->totalAtoms*sizeof(double));
	cudaMemcpy(&(dev_crys->positions),dev_positions,sizeof(double *));
	double *dev_totalEachSpecies;
	cudaMalloc((void **)&(dev_totalEachSpecies),maxelements*sizeof(int));
	cudaMemcpy(&(dev_crys->totalEachSpecies),dev_totalEachSpecies,sizeof(int *));
	char * dev_elements;
	cudaMalloc((void **)&(dev_elements),maxelements * namelength * sizeof(char));
	cudaMemcpy(&(dev_crys->elements),dev_elements,sizeof(char *));
	crystalnetwork * dev_network;
	cudaMalloc((void **)&(dev_network),sizeof(crystalnetwork));
	cudaMemcpy(&(dev_crys->network),dev_network,sizeof(crystalnetwork *));
	int * dev_adjacencyList;
	cudaMalloc((void **)&(dev_adjacencyList),maxconnections*crys->totalAtoms* sizeof(int));
	cudaMemcpy(&(dev_network->adjacencyList),dev_adjacencyList,sizeof(int *));
	int * dev_numAdjacent;
	cudaMalloc((void **)&(dev_numAdjacent),crys->totalAtoms * sizeof(int));
	cudaMemcpy(&(dev_network->numAdjacent),dev_numAdjacent,sizeof(int *));
	return dev_crys;
}

void CUDA_integratednthNearestNeighbors(crystal *crys, uint32_t *restrict integratednneighbors, int *restrict numintegratednneighbors, int *restrict shellatoms, int *restrict lastshellatoms, int index, int n)
{
	//printf("Location of pointers is %d, %d, %d, %d, %d\n",crys,integratednneighbors,numintegratednneighbors,shellatoms,lastshellatoms);
	//This function finds the all the nth nearest neighbors of an atom.
	//The nth nearest neighbors exclude atoms which were neighbors of any previous shell.
	int neighbors[maxconnections];
	int numneighbors;

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
	int neighbor;
	
	for(shellnumber = 1;shellnumber<=n;shellnumber++)
	{
		//printf("In shell %d of %d\n",shellnumber,n);
		numshellatoms=0;
		for(i=0;;i++)//Identify all atoms that are neighbors of atoms in last shell
		//for(i=0;(i<numlastshellatoms);i++)//Identify all atoms that are neighbors of atoms in last shell
		{
			//printf("Getting the nearest neighbors of atom %d\n",*(lastshellatoms+i));
			//printf("i %d of %d\n",i,numlastshellatoms);
            ACC_nearestNeighbors(crys,neighbors,&numneighbors,*(lastshellatoms+i));
            //printf("Found %d neighbors of atom %d\n",numneighbors,*(lastshellatoms+i));
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
            if(i+1 == numlastshellatoms)//The normal loop condition of the for loop was not working correctly. So I added a variable that tracks whether or not to loop.
				break;
		}
		
		shellswaptmp = shellatoms; //Swap pointers so that lastshellatoms is populated and shellatoms has space to work
		shellatoms = lastshellatoms;
		lastshellatoms = shellswaptmp;
		
		numlastshellatoms=numshellatoms;
	}
	
}

void CUDA_clusterApprox(int numGangs, int vectorLength, crystal *restrict crys, char *restrict coordElement, char *restrict countedElements, int numCountedElements, int *restrict coord, int coordinationNumber, int maxCoordination)
{
	//If passed a crystal, an element, and a series of other elements
	//This determines how many of each of the other elements surrounds the first element
	//The coordinationNumber can be used to find the coordination up to a further shell, for instance the next-nearest neighbors (coordinationNumber=2).

	dim3 dimGrid(1,1,1);
	dim3 dimBlock(1,1,1);
	crystal * dev_crys = CUDA_crys_copyToDevice(crys);
	CUDA_crys_printElements<<<dimGrid,dimBlock>>>(dev_crys);
	cudaDeviceSynchronize();
	return;
	


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
	
	uint32_t nneighbors[totalAtoms/32+1];
	
	int i;
	int numnneighbors,coordindex,iEle,eleCount;
	for(coordindex=0;coordindex<numCountedElements;coordindex++)
	{
	elementIndicies[coordindex] = crys_elementIndex(crys,countedElements+coordindex*namelength);
	for(i=0;i<numatoms;i++)
		*(coord+i*numCountedElements+coordindex) = 0;
	//printf("Element index for %s is %d\n",countedElements+namelength*coordindex,elementIndicies[coordindex]);
	}
	nvtxRangePushA("Cluster Approx");
	//for(iEle=0;iEle<totalElements+1;iEle++)
	//printf("Bound %d\n",elementBoundsArray[iEle]);
	//for(coordindex=0;coordindex<numCountedElements;coordindex++)
	//int eleCount;	
#pragma acc parallel loop gang num_gangs(numGangs) vector_length(vectorLength)\
private(coordindex,iEle,eleCount,numnneighbors)\
copyin(lastshellatoms[0:totalAtoms],shellatoms[0:totalAtoms],nneighbors[0:totalAtoms/32+1])\
private(lastshellatoms[0:totalAtoms],shellatoms[0:totalAtoms],nneighbors[0:totalAtoms/32+1])\
copyin(elementIndicies[0:numCountedElements],elementBoundsArray[0:totalElements+1])\
copyin(estart,maxCoordination,coordinationNumber,numCountedElements,maxconnections,totalAtoms)\
copyin(crys[0:1])\
copyin(crys->elements[0:namelength*totalElements],crys->positions[0:3*totalAtoms],crys->species[0:namelength*totalAtoms],crys->totalEachSpecies[0:totalElements],crys->network[0:1])\
copyin(crys->network->adjacencyList[0:maxconnections*totalAtoms],crys->network->numAdjacent[0:totalAtoms])\
copy(coord[0:numatoms*numCountedElements]) 
	for(i=0;i<numatoms;i++)
	{
		memset(nneighbors,0,sizeof(nneighbors));
		//printf("On atom %d which is %s\n",i+estart,crys->species+(i+estart)*namelength);
		//printf("%d %d %d\n",__pgi_gangidx(),__pgi_vectoridx(),i);
		//printf("Memory Locations %d %d %d %d\n",(uintptr_t)lastshellatoms,(uintptr_t)shellatoms,(uintptr_t)nneighbors,(uintptr_t)&numnneighbors);
		//printf("Before Species is located at %d\n",crys->species);
		ACC_integratednthNearestNeighbors(crys,nneighbors,&numnneighbors,shellatoms,lastshellatoms,i+estart,coordinationNumber);
		//printf("After Species is located at %d\n",crys->species);
		//printf("%d nneighbors\n",numnneighbors);
		for(coordindex = 0;coordindex<numCountedElements;coordindex++)
		{
			//*(coord+i*numCountedElements+coordindex) = 0;
			int iEle = elementIndicies[coordindex];
			int iAtom;
			int start=elementBoundsArray[iEle];
			int end=elementBoundsArray[iEle+1];
			eleCount = 0;
			//#pragma acc loop seq
			#pragma acc loop reduction(+:eleCount)
			for(iAtom=start;iAtom<end;iAtom++)
				eleCount += TestBit(nneighbors,iAtom);
			//printf("%d eleCount\n",eleCount);
			*(coord+i*numCountedElements+coordindex) = eleCount;
		}
		//printf("Finished counting neighbors\n");
	}
	//printf("Finished Coordination\n");
	nvtxRangePop();
}
