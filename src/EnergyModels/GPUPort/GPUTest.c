#include "lattikem.h"
#include "GPUCluster.h"

double latticeConstant = 3.905;
int sizex=12;
int sizey=12;
int sizez=12;
int coordinationNumber = 4;
int XLatticeShells[] = {6, 30, 84, 204, 354, 642, 936, 1464, 1950};//The number of X-site anion neighbors surrounding a B-site cation in perovskites up to i shells
double iodineRatio = 0.5;
double vacancyRatio = 0.0005;
int repeats = 100;//repeat the calculation several times to increase the time
int useGPU = 1;
enum method { CURRENT_METHOD,ACC_METHOD,BIT_ARRAY,CUDA_METHOD};
int numGangs = 14;
int vectorLength = 32;

int main(int argc, char **argv)
{
	if (argc == 5)
	{
		numGangs = atoi(argv[1]);
		vectorLength = atoi(argv[2]);
		repeats = atoi(argv[3]);
		useGPU = atoi(argv[4]);
	}
	else
	printf("Usage: GPUtest numGangs numRepeats useGPU\n Proceeding with default values of %d %d %d\n",numGangs, repeats, useGPU);
	
	int iRepeat;
	struct crystal *MAPbI = perovskite_newCrys(MA,"Pb","I","I","I",latticeConstant);
	struct crystal *crys = crys_multiply(MAPbI,sizez,sizey,sizez);
	//srand(10);//Changing the seed should alter the result of the follow command and resulting clusters.
	crys_replaceRandomAtoms(crys,3*(1-iodineRatio)*(sizex*sizey*sizez),"I",1,"Br");
	crys_replaceRandomAtoms(crys,3*(vacancyRatio)*(sizex*sizey*sizez),"I",1,"V");
	
	char *coordElement = "Pb";
	int numBandgapAlteringElements = 2;
	char *bandgapAlteringElements = crys_elementString(2,"Br","I");	
	
	struct NearestNeighborDescriptor *nnd = malloc(sizeof(NearestNeighborDescriptor));
	nnd->nndistance=0.75*latticeConstant;
	nnd->elementList=crys_elementString(4,"Pb","I","Br",VX);
	nnd->numEle=4;
	cn_allocateCrystalNetwork(crys);
	cn_fillFromnnd(crys,nnd); //fillFromnnd uses the more-efficient bucket method.
	//cn_fillNetwork(crys, nnd->nndistance, nnd->elementList, nnd->numEle);

	cn_allocatedSize(crys);
	crys_allocatedSize(crys);
    int numcoordinationatoms = crys_elementCount(crys,coordElement);
	int maxCoordination = XLatticeShells[coordinationNumber-1];
	
	int coord[numcoordinationatoms*numBandgapAlteringElements];
	int bins[maxCoordination+1];
	int iCoord,iBin;
		
	//crys_printAllAtoms(crys);
	crys_printElements(crys);
	//cn_printAdjacencyList(crys);
    printf("determining the coordination for %d coordination atoms from %d elements which are %s using coordination number %d and max coordination %d\n",numcoordinationatoms,numBandgapAlteringElements,bandgapAlteringElements,coordinationNumber,maxCoordination);
	
	if(useGPU == ACC_METHOD)
	{
		printf("OpenACC using %d gangs %d vector length and performing %d repeats\n",numGangs,vectorLength,repeats);
		for(iRepeat=0;iRepeat<repeats;iRepeat++)
		ACC_clusterApprox(numGangs,vectorLength,crys,coordElement,bandgapAlteringElements,numBandgapAlteringElements,coord,coordinationNumber,maxCoordination);
		
		
	}
	if(useGPU == CUDA_METHOD)
	{
		printf("CUDA using %d gangs %d vector length and performing %d repeats\n",numGangs,vectorLength,repeats);
		for(iRepeat=0;iRepeat<repeats;iRepeat++)
		CUDA_clusterApprox(numGangs,vectorLength,crys,coordElement,bandgapAlteringElements,numBandgapAlteringElements,coord,coordinationNumber,maxCoordination);
		
		
	}
	
	if(useGPU == CURRENT_METHOD)
	{
		for(iRepeat=0;iRepeat<repeats;iRepeat++)
		cn_coordination(crys,coordElement,bandgapAlteringElements,numBandgapAlteringElements,coord,coordinationNumber,maxCoordination);
	}
	if(useGPU == BIT_ARRAY)
	{
		for(iRepeat=0;iRepeat<repeats;iRepeat++)
		CPU_clusterApprox(crys,coordElement,bandgapAlteringElements,numBandgapAlteringElements,coord,coordinationNumber,maxCoordination);
	}
	
	printf("Done calculating coordination\n");
	printf("Counting Clusters\n");
			
	for(iBin=0;iBin<maxCoordination+1;iBin++)
		bins[iBin]=0;
	for(iCoord=0;iCoord<numcoordinationatoms;iCoord++)
		bins[*(coord+numBandgapAlteringElements*iCoord)]++;
	printf("#Br #Clusters\n");
	for(iBin=0;iBin<maxCoordination+1;iBin++)
		if(bins[iBin]!=0)
			printf("%d %d\n",iBin,bins[iBin]);
	
	return 0;
}
