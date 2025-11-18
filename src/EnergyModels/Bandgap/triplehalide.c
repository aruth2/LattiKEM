/*
 * mixedhalide.c
 * 
 * Copyright 2019 Anthony Ruth <aruth@heesolar.com>
 * mixedhalide.c simulates the vacancy-mediated movement of I and Br atoms in a perovskite lattice
 * considering the local differences in bandgap surrounding I and Br atoms. The difference in bandgap 
 * is assumed to influence ion migration.
 */

#include "perovskite.h"
#include "bandgap.h"
//#include "crystal.h"

//char dir[1000];
double latticeConstant = 6.36;//This could really use something more clever
double ClRatio, BrRatio, vacancyRatio;
double ClHopEnergy, BrHopEnergy, IHopEnergy;
double IBrRepulsiveEnergy, BrClRepulsiveEnergy, IClRepulsiveEnergy;//Energy is per halide atom. If all eight of a I atoms neighbors are Br it will be assessed the full energy, otherwise it will be proportionally lower. 
double pureIBandgap, pureBrBandgap, pureClBandgap, IBrBowingParameter, BrClBowingParameter, IClBowingParameter;

int XLatticeShells[] = {6, 30, 84, 204, 354, 642, 936, 1464, 1950,2790,3516,4740,5754,7434,8784,10992,12726,15534,17700};//The number of X-site anion neighbors surrounding a B-site cation in perovskites up to i shells
//int repulsiveCoordinationNumber;
int calculateLatticeShells;
int	sizex,sizey,sizez;
int pbcMask;//a 3 bit mask for pbc in x,y and z
double th_gapEnergy(int numBandgapAlteringElements, int *numEachElement);
void th_setup();
void th_saveSettings();
void th_loadSettings(char *filename);
crystal *th_crystal();
void th_trajectory(Trajectory *traj);

double th_gapEnergy(int numBandgapAlteringElements, int *numEachElement)
{
	//For MAPb(I1-x-yBrxCly)3
	//E(x)  = 0.39x + 0.33x^2+1.57   DOI: 10.1021/nl400349b
	//Eg(x,y) = Eg_I + x * (DeltaEg_IBr) + y * (DeltaEg_BrCl + DeltaEg_IBr)
	
	int numCl,numBr,numI;
	
	numCl = *(numEachElement);
	numBr = *(numEachElement+1);
	numI = *(numEachElement+2);
	double x = (double)numBr/(numCl+numBr+numI);
	double y = (double)numCl/(numCl+numBr+numI);
	double z = (double)numI/(numCl+numBr+numI);
	
	return x*pureBrBandgap + y*pureClBandgap + z*pureIBandgap - IBrBowingParameter*x*z - IClBowingParameter*y*z - BrClBowingParameter*x*y;
}

void th_trajectory(Trajectory *traj)
{
    /*
     * Provides a description of hopping barriers for a mixed halide crystal
     * The barriers can be element specific.  
     * */
   
	LatticeDynamics *LD = traj->LD = (LatticeDynamics *)malloc(sizeof(LatticeDynamics));
	bg_trajectory(traj);
	
   
	LD->numHopPairs = 3;
	LD->hopPairs = crys_elementString(6,VX,"Br",VX,"Cl",VX,"I");	
	LD->hopPairEnergies = (double *)malloc(3*sizeof(double));
	
	*(LD->hopPairEnergies) = BrHopEnergy; //Br
	*(LD->hopPairEnergies+1) = ClHopEnergy; //Cl
	*(LD->hopPairEnergies+2) = IHopEnergy; //I
	 
    traj->crys = th_crystal();
    traj->crys->network->maxconnections = 2*XLatticeShells[0];//Not sure if this is always enough
     
	(traj->nnd)->nndistance=0.75*latticeConstant;
	(traj->nnd)->elementList=crys_elementString(5,"Pb","I","Br","Cl",VX);
	(traj->nnd)->numEle=5;
	traj->numnnds=1;
}

crystal *th_crystal()
{
    //Setup Crystal
	struct crystal *MAPbI = perovskite_newCrys(MA,"Pb","I","I","I",latticeConstant);
	struct crystal *crys = crys_multiply(MAPbI,sizez,sizey,sizez);
	crys_replaceRandomAtoms(crys,3*(BrRatio)*(sizex*sizey*sizez),"I",1,"Br");	
	crys_replaceRandomAtoms(crys,3*(ClRatio)*(sizex*sizey*sizez),"I",1,"Cl");	
		
	int numHalideVacancies = vacancyRatio*(crys->totalAtoms)*3.0/5.0;
	if(numHalideVacancies < 2)
	numHalideVacancies = 2;
	printf("%d total atoms. Replacing %d atoms by vacancies actual concentration %f \n",crys->totalAtoms,numHalideVacancies,(double)numHalideVacancies/(crys->totalAtoms*3.0/5.0));
	crys_replaceRandomAtoms(crys,numHalideVacancies,crys_elementString(3,"I","Br","Cl"),3,VX);

    crys_printElements(crys);
	
	if(!(pbcMask % 2) )//x pbc
	{
		printf("Disabled x pbc\n");
		*(crys->latticeVectors) *=2;
	}
	if(!(pbcMask >> 1 % 2) )//y pbc
	{
		printf("Disabled y pbc\n");
		*(crys->latticeVectors) *=2;
	}
	if(!(pbcMask >> 2 % 2) )//z pbc
	{
		printf("Disabled z pbc\n");
		*(crys->latticeVectors) *=2;
	}
	
    return crys;
}

void TH_XLatticeShells()
{
	printf("Calculating size of shells of X lattice\n");
	Trajectory *traj = (Trajectory *)malloc(sizeof(Trajectory));
	th_trajectory(traj);
	
	crystal *crys = traj->crys;
	cn_allocateCrystalNetwork(traj->crys);
	for(int innd=0;innd<traj->numnnds;innd++)
	{
		cn_printnnd(traj->nnd+innd);
		cn_fillFromnnd(traj->crys,traj->nnd+innd);
	}
	int iNeighbor;
	int numneighbors;
	//uint32_t *neighbors = malloc(roundup(crys->totalAtoms,32)/32);
	uint32_t *neighbors = (uint32_t *)malloc(1000000*sizeof(uint32_t));

	int iAtom1 = crys_elementOffset(crys, "Pb");
	/*
	printf("Num adjacent to atom %d, %d\n",iAtom1,*(crys->network->numAdjacent+iAtom1));
	int iAdjacent;
	for(iAdjacent = 0;iAdjacent< *(crys->network->numAdjacent+iAtom1);iAdjacent++)
		printf("Atom %d\n",*(crys->network->adjacencyList+iAtom1*crys->network->maxconnections+iAdjacent));
	*/
	
	int shellatoms[crys->totalAtoms];
	int lastshellatoms[crys->totalAtoms];
	for(iNeighbor=1;iNeighbor<20;iNeighbor++)
	{
		//memset(neighbors,0,roundup(crys->totalAtoms,32)/32);
		memset(neighbors,0,1000000*sizeof(uint32_t));
		//numneighbors=0;
		//printf("On atom %d which is %s\n",i+estart,crys->species+(i+estart)*namelength);

		cn_bitA_integratednthNearestNeighbors(crys,crys->network,neighbors,&numneighbors,shellatoms,lastshellatoms,iAtom1,iNeighbor);		
		//printf("Finished integrated NN\n");

		int iAtom2;
		int start=crys_elementOffset(crys,"I");
		int end=start+crys_elementCount(crys,"I");
		int eleCount = 0;
		//printf("Counting Neighbors between %d and %d\n",start,end);
		for(iAtom2=start;iAtom2<end;iAtom2++)
			eleCount += TestBit(neighbors,iAtom2);
		
		start=crys_elementOffset(crys,"Xe");
		end=start+crys_elementCount(crys,"Xe");
		//printf("Counting Neighbors between %d and %d\n",start,end);
		for(iAtom2=start;iAtom2<end;iAtom2++)
			eleCount += TestBit(neighbors,iAtom2);
		
		printf("Shell %d num %d of %d\n",iNeighbor,eleCount,numneighbors);
	}
	exit(0);
}

void th_registerSettings()
{
	
	
	//registerString(dir,"dir",".");
    registerDouble(&(vacancyRatio),"vacancyRatio",0.0);
    registerDouble(&(BrRatio),"BrRatio",0.5);
    registerDouble(&(ClRatio),"ClRatio",0.25);
    registerDouble(&(ClHopEnergy),"ClHopEnergy",0.25);
    registerDouble(&(BrHopEnergy),"BrHopEnergy",0.25);
    registerDouble(&(IHopEnergy),"IHopEnergy",0.25);
    registerDouble(&(IBrRepulsiveEnergy),"IBrRepulsiveEnergy",0.0);
    registerDouble(&(BrClRepulsiveEnergy),"BrClRepulsiveEnergy",0.0);
    registerDouble(&(IClRepulsiveEnergy),"IClRepulsiveEnergy",0.0);
    //For MAPb(I1-xBrx)3
	//E(x)  = 0.39x + 0.33x^2+1.57   DOI: 10.1021/nl400349b
    registerDouble(&(pureIBandgap),"pureIBandgap",1.57);
    registerDouble(&(pureBrBandgap),"pureBrBandgap",2.29);
    registerDouble(&(pureClBandgap),"pureClBandgap",3.00);
    registerDouble(&(IBrBowingParameter),"IBrBowingParameter",0.33);
    registerDouble(&(BrClBowingParameter),"BrClBowingParameter",0.33);
    registerDouble(&(IClBowingParameter),"IClBowingParameter",0.33);
    registerInt(&(sizex),"sizex",4);
    registerInt(&(sizey),"sizey",4);
    registerInt(&(sizez),"sizez",4);
    registerInt(&(pbcMask),"pbcMask",7);
    registerInt(&(calculateLatticeShells),"calculateLatticeShells",0);

    bg_registerSettings();
}

void th_setup()
{
	//Configuration for bandgap based on halide concentration.
	char *coordElement = "Pb";
    char *bandgapAlteringElements = crys_elementString(3,"Br","Cl","I");	
	int numBandgapAlteringElements = 3;
	int numStates = sizex*sizey*sizez;
	
	//Configuration for interatomic attracttion/repulsion. 
	//The plan is to gradually expand this functionality so values can be assigned to arbitrary coordination sheels (e.g. nearest neighbor, next-nearest neighbor). 
	//The code should also allow for an arbitrary number of pairs.
	//For n repulsive elements, there should be n*(n-1)/2 repulsiveEnergies.
	//int repulsiveCoordinationNumber = 8;
	//int maxRepulsiveCoordination = 10;
	char *repulsiveElements = crys_elementString(3,"Br","Cl","I");
	double **repulsiveEnergies = (double **)calloc(3,sizeof(double *));
	*(repulsiveEnergies) = &BrClRepulsiveEnergy;
	*(repulsiveEnergies+1) = &IBrRepulsiveEnergy;
	*(repulsiveEnergies+2) = &IClRepulsiveEnergy;
	int numRepulsiveElements = 3;

	
	bg_setup(th_gapEnergy,coordElement,bandgapAlteringElements,numBandgapAlteringElements,XLatticeShells,numStates,repulsiveElements,numRepulsiveElements,repulsiveEnergies,th_trajectory);
}

int main(int argc, char **argv)
{
	allocateSettings();
	th_registerSettings();
	
		
	FILE *infile = fopen(argv[1],"r");
	loadSettings(infile);
	fclose(infile);
	mkdir2(getDir());
	 
	char outfileName[1000];
	sprintf(outfileName,"%s/settings",getDir());
	FILE *outfile = fopen(outfileName,"w");
	saveSettings(outfile);
	fclose(outfile);


	th_setup();
	if(calculateLatticeShells)
		TH_XLatticeShells();
		
    simulateTrajectories();

	
	return 0;
}

