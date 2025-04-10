#include "bandgap.h"

char *bandgapAlteringElements;	
int numBandgapAlteringElements;
char *coordElement;
int	coordinationNumber;//The number of shells used to define the local bandgap.
int maxCoordination;
int numStates;
int numOpticsEnergies;
//int useGPU;
int numGangs;
int vectorLength;

//double	temperature;//In units of eV There could be seperate electronic and ionic temperatures
double	numExcitations;
double	numExcitationsSecondStep;
double optics_llimit = 1.0;
double optics_ulimit = 2.5;
double *shoulder, *shoulderenergies;
double shoulder_llimit, shoulder_ulimit;
double (*bandgapFunction)(int numBandgapAlteringElements, int *numEachElement);
int biasTurnOn;
int biasSwitch;
int thermalDistribution;
int BGWave;
double fermiConvergence;

char *repulsiveElements;
double *repulsiveEnergies;
int	numRepulsiveElements;
int	repulsiveShellSize;
int	maxRepulsiveCoordination;
int saveChargeDensity;
double weightCutoff;

void absorptionshoulder(double broadeningenergy, double Eb, double *shoulder,double *shoulderenergies, double *shoulder_llimit, double *shoulder_ulimit)
{
    /*
     * This calculates an exciton + electron-phonon absorptivity shoulder using the Elliot model and following the general form
     * prescribed by Saba et al. DOI:10.1038/ncomms6049
     * */
    
//	printf("Calculating absorption shoulder with %g broadening and %d optics points\n",broadeningenergy,numOpticsEnergies);
	*shoulder_llimit = -3*broadeningenergy;
	*shoulder_ulimit = 4.5 + 3*broadeningenergy;
	double absorptionpeaks[2*numOpticsEnergies];
	double peakenergies[2*numOpticsEnergies];
	int j;

	absorptionpeaks[0] = 0;
	peakenergies[0] =0;
	for(j=1;j<numOpticsEnergies;j++)
	{
		peakenergies[j] = -Eb/pow(j,2);
		absorptionpeaks[j] = 4*M_PI*pow(Eb,1.5)/pow(j,3);
		//printf("Delta peak at %g,%g\n",peakenergies[j],absorptionpeaks[j]);
	}
	for(j=0;j<numOpticsEnergies;j++)
	{
		peakenergies[j+numOpticsEnergies] = (*shoulder_ulimit-*shoulder_llimit)*j/(numOpticsEnergies-1);
		absorptionpeaks[j+numOpticsEnergies] = (*shoulder_ulimit-*shoulder_llimit)/numOpticsEnergies* 2*M_PI*sqrt(Eb)/(1-exp(-2*M_PI*sqrt((Eb)/peakenergies[j+numOpticsEnergies])));
//		printf("Theta peak at %g,%g\n",peakenergies[j+numOpticsEnergies],absorptionpeaks[j+numOpticsEnergies]);
	}
	//broadening(peakenergies,absorptionpeaks,shoulderenergies,shoulder,2*numOpticsEnergies,numOpticsEnergies,*shoulder_llimit,*shoulder_ulimit,broadeningenergy,1,0);
	gaussianconvolution(shoulderenergies, shoulder, peakenergies, absorptionpeaks, numOpticsEnergies, 2*numOpticsEnergies, *shoulder_llimit,*shoulder_ulimit,  broadeningenergy,1);
	
	for(j=0;j<numOpticsEnergies;j++)
	{
		//printf("%g %g\n",shoulderenergies[j],shoulder[j]);
	}
}

void OS_save(Configuration *config)
{
    
	FILE *outfile = fopen(config->dataFileName,"w");
	OptoelectronicState *data = config->data;
	int i;
	for(i=0;i<numOpticsEnergies;i++)
	fprintf(outfile,"%g %g %g %g %g\n",*(data->photonenergies+i),*(data->absorption+i),*(data->emission+i),*(data->DOS+i),*(data->current+i));
	fclose(outfile);
	if(saveChargeDensity)
	{
		char *filename = strcat2(config->dataFileName,".situs");
		outfile = fopen(filename,"w");
		free(filename);
		int dimensions = pow(numStates+0.001,1/3.0);//Only works for cubic cells
		fprintf(outfile,"%.06f %.06f %.06f %.06f %d %d %d \n\n",data->voxelSize,data->voxelSize/2,data->voxelSize/2,data->voxelSize/2,dimensions,dimensions,dimensions);
		int counter = 0;
		
		char buf[64];
		int d,e;
		//This assumes the data (atom locations) are correctly oriented in first X, then Y, then Z
		for(i = 0; i < numStates; i++)
		{
			d = sprintf(buf, "%.06f",*(data->chargeDensity+i));
			fprintf(outfile,"%*s%.06f ",11-d,"",*(data->chargeDensity+i));
			//printf("%g %.06f\n",*(data->chargeDensity+i),*(data->chargeDensity+i));
			//fprintf(outfile,"   %.06f ",*(data->chargeDensity+i));
			counter++;
			if(counter == 10)
			{
				fprintf(outfile,"\n");
				counter = 0;
			}
		}
		fclose(outfile);
	}
}

void OS_combineWeighted(Configuration *configs, int numCombine, Configuration *outconfig, double *weights)
{
	OptoelectronicState *data = outconfig->data;
	OptoelectronicState *OS;
	int numScalars = data->numScalars;
	int iScalar, iOE, iConfig, iState;
	double sumWeight = 0;
	data->voxelSize = 0;

	for(iConfig = 0; iConfig < numCombine; iConfig++)
		sumWeight += *(weights + iConfig);
	for(iScalar = 0; iScalar < numScalars; iScalar++)
		*(*(data->scalars + iScalar)) = 0;

	for(iScalar = 0; iScalar < numScalars; iScalar++)
	{
		for(iConfig = 0; iConfig < numCombine; iConfig++)
		{
			OS = (OptoelectronicState *)((configs + iConfig)->data);
			*(*(data->scalars + iScalar)) += *(weights + iConfig)* *(*(OS->scalars + iScalar));
			
		}
		*(*(data->scalars + iScalar)) /= sumWeight;
		*(outconfig->seriesData + iScalar) = *(*(data->scalars + iScalar));
	}
	
	for(iConfig = 0; iConfig < numCombine; iConfig++)
	{
		OS = (OptoelectronicState *)((configs + iConfig)->data);
		data->voxelSize += *(weights + iConfig)* OS->voxelSize/sumWeight;	
	}
	
	for(iState = 0; iState < numStates; iState++)
	{
		*(data->chargeDensity + iState) = 0;
		for(iConfig = 0; iConfig < numCombine; iConfig++)
		{
			OS = (OptoelectronicState *)((configs + iConfig)->data);
			*(data->chargeDensity + iState) += *(weights + iConfig)* *(OS->chargeDensity + iState);
			//printf("%g\n",*(OS->chargeDensity + iState));
		}
		*(data->chargeDensity + iState) /= sumWeight;
		//printf("\n%g\n",*(data->chargeDensity + iState));
	}

	for(iOE = 0; iOE <numOpticsEnergies;iOE++)
	{
		*(data->photonenergies+iOE)=0;
		*(data->absorption+iOE)=0;
		*(data->emission+iOE)=0;
		*(data->DOS+iOE)=0;
		*(data->current+iOE)=0;
		for(iConfig = 0;iConfig < numCombine; iConfig++)
		{
			OS = (OptoelectronicState *)((configs + iConfig)->data);
			*(data->photonenergies+iOE) += *(weights + iConfig) * *(OS->photonenergies+iOE);
			*(data->absorption+iOE) += *(weights + iConfig) * *(OS->absorption+iOE);
			*(data->emission+iOE) += *(weights + iConfig) * *(OS->emission+iOE);
			*(data->DOS+iOE) += *(weights + iConfig) * *(OS->DOS+iOE);
			*(data->current+iOE) += *(weights + iConfig) * *(OS->current+iOE);
		}
		*(data->photonenergies+iOE)/= sumWeight;
		*(data->absorption+iOE)/= sumWeight;
		*(data->emission+iOE)/= sumWeight;
		*(data->DOS+iOE)/= sumWeight;
		*(data->current+iOE)/= sumWeight;
	}
	
}


void OS_allocate(OptoelectronicState *OS)
{
    OS->energystates = malloc(numStates*sizeof(double));
    OS->weights = malloc(numStates*sizeof(double));
    OS->chargeDensity = malloc(numStates*sizeof(double));
    OS->emission = malloc(numOpticsEnergies*sizeof(double));
    OS->DOS = malloc(numOpticsEnergies*sizeof(double));
    OS->current = malloc(numOpticsEnergies*sizeof(double));
    OS->absorption = malloc(numOpticsEnergies*sizeof(double));
    OS->photonenergies = malloc(numOpticsEnergies*sizeof(double));
    
   //int numcoordinationatoms = crys_elementCount(crys,coordElement);
    //double energy;
	//int coord[numcoordinationatoms*numBandgapAlteringElements];
	OS->bandgapNetwork = NULL;
	OS->coord = calloc(numStates*numBandgapAlteringElements,sizeof(int));
    
    //Scalars
	//This is somewhat redundant in that the number of scalars should match the number of series data. It should only be set once
	OS->numScalars=3;
	OS->scalars=calloc(OS->numScalars,sizeof(double *));
	*(OS->scalars) = &OS->fermiEnergy;
	*(OS->scalars+1) = &OS->photocarrierEnergy;
	*(OS->scalars+2) = &OS->interatomicEnergy;
}

void OS_free(OptoelectronicState *OS)
{
    free(OS->energystates);
    free(OS->weights);
    free(OS->chargeDensity);
    free(OS->emission);
    free(OS->absorption);
    free(OS->photonenergies);
    free(OS->DOS);
    free(OS->current);
}

void fermiSum( double *energyLevels, double fermiEnergy, double temperature, int numEnergyLevels, double *weights, double *weightSum, double *energySum)
{
 	/*
	Evaluates the fermi Sum, N=sum(fi) and E=sum(fi * Ei)
	where fi is the fermi function fi = 1/(1+exp((Ei-Ef)/kT))
	
	Inputs are:
	energyLevels (Ei)
	fermiEnergy (Ef)
	temperature (kT)
	numEnergyLevels (n)
	
	outputs are:
	weights (fi)
	weightSum (N)
	energySum (E)
	*/
	*weightSum=0;
	*energySum=0;
	for(int i = 0;i<numEnergyLevels,i++)
	{
	*(weights+i) = 1.0/(1.0+exp((*(energyLevels+i)-fermiEnergy)/temperature));
	*weightSum += *(weights+i);
	*energySum += *(weights+i) * *(energyLevels+i);
	}
}

double solveFermiEnergy(double currentExcitations, double fermiGuess, double fermiConvergence, double *energyLevels,double *weights)
{
		double fermiMin = fermiGuess - 7.5*fermiConvergence;
		double fermiMax = fermiGuess + 8.49*fermiConvergence;//This is slightly offset so that the fermi level could be almost unchanged. 
		double fermiEnergy;
		double totalCarriersMin=0;
		double totalCarriersMax=0; 
		double energy=0;
		//Try to bound the fermi level within a small box around the previous result which requires not more than 6 iterations to converge.
		fermiSum( energyLevels, fermiMin, temperature, numcoordinationatoms, weights, &totalCarriersMin, &energy);
		fermiSum( energyLevels, fermiMin, temperature, numcoordinationatoms, weights, &totalCarriersMax, &energy);
		
		if((totalCarriersMin > currentExcitations) || (totalCarriersMax < currentExcitations))//If the solution is not within the bounds expand the bounds.  	
		{
		//printf("Fermi level outside bounds, Min %g max %g, min carriers %g, max carriers %g, target carriers %g\n",fermiMin,fermiMax,totalCarriersMin,totalCarriersMax,currentExcitations);
		fermiMin = optics_llimit;
		fermiMax = optics_ulimit;
		}
		
		double totalCarriers=0;
		while(fermiMax-fermiMin > fermiConvergence)
		{
			fermiEnergy = (fermiMax+fermiMin)/2.0;	
			fermiSum( energyLevels, fermiMin, temperature, numcoordinationatoms, weights, &totalCarriersMax, &energy);
			if(totalCarriers > currentExcitations) //Lower search range
				fermiMax = fermiEnergy;
			else
				fermiMin = fermiEnergy;
		}
	return fermiEnergy;
}

double OS_energy(crystal *crys, int *coord,Configuration *config)
{
	//This function calculates the energy of n band edge excitations based on local composition. 
    //First local bandgaps are identified and then excitons are distributed to preferentially reside
    //in low bandgap regions.
	//The probability is considered via boltzmann or fermi distribution
	
	int numcoordinationatoms = crys_elementCount(crys,coordElement);	
    double currentExcitations =  *(config->externalConditions+0);
	double temperature = getTemperature();
	
    OptoelectronicState *OS = config->data;
    //printf("Dat is located at %d\n",(int *)dat);
	int i;

	double lowestenergystate = 10;
	double highestenergystate = 0;
	//This should be updated for partial sums. Only those states which change need to be calculated.
	for(i=0;i<numcoordinationatoms;i++)
	{
		*(OS->energystates+i) = bandgapFunction(numBandgapAlteringElements, coord+numBandgapAlteringElements*i);

        if(*(OS->energystates+i)<lowestenergystate)
            lowestenergystate=*(OS->energystates+i);
            
		if(*(OS->energystates+i)>highestenergystate)
            highestenergystate = *(OS->energystates+i);
	}
	double energy = 0;
	
	if(thermalDistribution == BOLTZMANN_DISTRIBUTION)
	{
		double probabilitysum = 0;
		for(i=0;i<numcoordinationatoms;i++)
		{
			//By using the middle energy we keep the weights close to 1 to reduce floating point errors.
			*(OS->weights+i) = exp(-(*(OS->energystates+i) - (highestenergystate+lowestenergystate)/2.0)/temperature);
			probabilitysum += *(OS->weights+i);
			//printf("Energy state %g weight %g probability sum %g\n",*(OS->energystates+i),*(OS->weights+i),probabilitysum);
		}
		
		//printf("Lowest energy state %g highest %g\n",lowestenergystate,highestenergystate);
		
		for(i=0;i<numcoordinationatoms;i++)
		{
			*(OS->weights+i) = *(OS->weights+i) * currentExcitations/probabilitysum;
			//energy += phaseenergy[i] * exp(-phaseenergy[i]/temperature) * numExcitations/probabilitysum;
			energy += *(OS->energystates+i) * *(OS->weights+i);
		    //printf("Energy state %g weight %g probability sum %g\n",*(OS->energystates+i),*(OS->weights+i)/numExcitations,probabilitysum);
		}
	}
	
	if(thermalDistribution == FERMI_DISTRIBUTION)
	{
		OS->fermiEnergy =  solveFermiEnergy(currentExcitations, OS->fermiEnergy, fermiConvergence, OS->energystates,OS->weights);
	}
	
	
    //Only calculate optics if they will be saved. This helps a lot on kinetic steps
    if(config->savedata)
    {
		//printf("Here\n");
		convolve(OS->photonenergies,OS->absorption, OS->energystates,OS->weights, shoulderenergies, shoulder, numOpticsEnergies, numcoordinationatoms, numOpticsEnergies, optics_llimit,optics_ulimit,0);
		gaussianconvolution(OS->photonenergies,OS->emission, OS->energystates,OS->weights, numOpticsEnergies, numcoordinationatoms, optics_llimit,optics_ulimit,temperature,1);
		gaussianconvolution(OS->photonenergies,OS->DOS, OS->energystates,OS->weights, numOpticsEnergies, numcoordinationatoms, optics_llimit,optics_ulimit,temperature,0);
		stepconvolution(OS->photonenergies,OS->current, OS->photonenergies,OS->emission, numOpticsEnergies, numOpticsEnergies, optics_llimit,optics_ulimit,1,1);
    }
    if(saveChargeDensity)
    {
		int coordOffset = crys_elementOffset(crys,coordElement);
		OS->voxelSize = crys_atomDistance(crys,coordOffset,coordOffset+1);
		double radiusCutoff = OS->voxelSize*pow(numStates*3/4/M_PI,1/3.0)+0.001;
		for(i=0;i<numcoordinationatoms;i++)
			*(OS->chargeDensity+i) = 0;
		
		for(i=0;i<numcoordinationatoms;i++)
		{
			int iAtom;
			if(*(OS->weights+i)>weightCutoff)
			for(iAtom=0;iAtom<numcoordinationatoms;iAtom++)
			{
				if(crys_atomDistance(crys,i+coordOffset,iAtom+coordOffset)<radiusCutoff)
					*(OS->chargeDensity+iAtom) += *(OS->weights+i);
			}
		}
		//for(i=0;i<numcoordinationatoms;i++)
			//printf("%g\n",*(OS->chargeDensity+i));
		normalize(OS->chargeDensity,numStates,currentExcitations);
	}
    
    int iScalar;
	//This is kind of a ridiculous ask
	if(config->seriesData!=NULL)
	for(iScalar = 0;iScalar<OS->numScalars;iScalar++)
		*(config->seriesData+iScalar) = *(*(OS->scalars+iScalar));
    //config->energy = energy;
    
	return energy;
}

void OS_printCoord(OptoelectronicState *OS)
{
	//OptoelectronicState *OS = config->data;
	int *coord = OS->coord;
	int iState, iEle, sum;
	for(iState=0;iState<numStates;iState++)
	{
		printf("State %d\n",iState);
		sum = 0;
		for(iEle=0;iEle<numBandgapAlteringElements;iEle++)
			{
			printf("%s\t%d\t\t",bandgapAlteringElements+iEle*namelength,*(coord+iState*numBandgapAlteringElements+iEle));
			sum += *(coord+iState*numBandgapAlteringElements+iEle);
			}
		printf("%d tot\n",sum);
		if(sum > OS->bandgapNetwork->maxconnections)
			exit(1);
	}
}

double interatomic_energy(Configuration *config)
{
	//Calculates the energy of n*m neaarest neighbor atomic pairs based on an energy of repulsion/attraction U12
	crystal *crys = config->crys;
	//OptoelectronicState *data = config->data;
	
    int numcoordinationatoms=0;
    int iele1,iele2;
    int repulsivePair=0;
    int numEle1,numEle2;
    int numElementNeighbors;
    int iAtom1;
    double energy=0;
    
    int atomCount;
    for(iele1=0;iele1<numRepulsiveElements;iele1++)
    {
		atomCount = crys_elementCount(crys,repulsiveElements+namelength*iele1);
		if(atomCount>numcoordinationatoms)
		numcoordinationatoms=atomCount;
	}
    int coord[numcoordinationatoms*(numRepulsiveElements-1)];
    
    //printf("computing interatomic energy from %s with %d repulsive Elements with %d max shell\n",repulsiveElements,numRepulsiveElements,maxRepulsiveCoordination);
    for(iele1=0;iele1<numRepulsiveElements-1;iele1++)
    {
		//printf("Counting element %s\n",repulsiveElements+namelength*iele1);
		numEle1 = crys_elementCount(crys,repulsiveElements+namelength*iele1);
		numElementNeighbors = numRepulsiveElements-iele1-1;
		//printf("Computing coordination of element %s which has %d element neighbors\n",repulsiveElements+namelength*iele1,numElementNeighbors);
		//cn_coordination(crys,repulsiveElements+namelength*iele1,repulsiveElements+namelength*(iele1+1),numElementNeighbors,coord,1,maxRepulsiveCoordination);//Implicitly nearest neighbors (coordination #1)
		
		//CPU_clusterApprox(crys,repulsiveElements+namelength*iele1,repulsiveElements+namelength*(iele1+1),numElementNeighbors,coord,1,maxRepulsiveCoordination);//Implicitly nearest neighbors (coordination #1)
		cn_bitA_clusterApprox(crys,crys->network,repulsiveElements+namelength*iele1,repulsiveElements+namelength*(iele1+1),numElementNeighbors,coord,1,maxRepulsiveCoordination);//Implicitly nearest neighbors (coordination #1)
		//printf("Finnished Computing coordination of element %s\n",repulsiveElements+namelength*iele1);
		for(iele2 = iele1+1;iele2<numRepulsiveElements;iele2++)
		{
			numEle2 = 0;
			//printf("Summing neighbors of element %s and %s\n",repulsiveElements+namelength*iele1,repulsiveElements+namelength*iele2);
			for(iAtom1=0;iAtom1<numEle1;iAtom1++)
				numEle2 += coord[iAtom1*numElementNeighbors+(iele2-iele1-1)];
			//printf("Found %d neighbors of type %s around atom %s which has %d members. This is repulsive pair %d\n",numEle2,repulsiveElements+namelength*iele2,repulsiveElements+namelength*iele1,numEle1,repulsivePair);
			//This implictly assumes that repulsiveShellSize is the number of nearest neighbors for each atom type. However, if some shells are larger than others, this will not work and we will need separate shell sizes for different atom types.
			//printf("Adding %g to energy\n",*(repulsiveEnergies+repulsivePair)*numEle2/repulsiveShellSize);
			energy += *(repulsiveEnergies+repulsivePair)*numEle2/repulsiveShellSize;
			repulsivePair++;
		}
	}
	
	return energy;
}
/*
double bg_energy_old(Configuration *config)
{
	//Calculates the energy of n band edge excitations
	//The excitation energy is estimated based on local bandgaps and varies from place to place
	crystal *crys = config->crys;
	OptoelectronicState *data = config->data;
	//int bins[maxCoordination*numBandgapAlteringElements];
    int numcoordinationatoms = crys_elementCount(crys,coordElement);
    double energy;
	
	
	int coord[numcoordinationatoms*numBandgapAlteringElements];
	
	cn_bitA_clusterApprox(crys,coordElement,bandgapAlteringElements,numBandgapAlteringElements,coord,coordinationNumber,maxCoordination);
	
	data->interatomicEnergy = interatomic_energy(config);
	data->photocarrierEnergy = OS_energy(crys,coord,config);
	energy = data->interatomicEnergy + data->photocarrierEnergy;

	//printf("interatomic energy is %g, photocarrier Energy is %g, total Energy is %g saving it to %d\n",data->interatomicEnergy, data->photocarrierEnergy,energy,(config->energy));
	(config->energy) = energy;
    
	return energy;
}*/
void bg_coordSwap(Configuration *config, int atom1, int atom2, int forward_reverse)
{
	//1 for forward, -1 to go in reverse and undo a swap.
	
	crystal *crys = config->crys;
	//printf("Coord swap %d %s %d %s\n",atom1,crys->species+atom1*namelength,atom2,crys->species+atom2*namelength);
	OptoelectronicState *data = config->data;
	crystalnetwork *cn = data->bandgapNetwork;
	int *coord = data->coord;
	if(data->bandgapNetwork == NULL)
		return;
	
	//printf("Swapping %d and %d\n",atom1,atom2);
	//OS_printCoord(data);	
	//cn_printAdjacencyList(crys,cn);
	
	
	int *adjacencyList = cn->adjacencyList;
	int *numAdjacent = cn->numAdjacent;
	int maxconnections = cn->maxconnections;
	
	int estart = crys_elementOffset(crys,coordElement);
	int numatoms = crys_elementCount(crys,coordElement);	
	
	int iEle1 = crys_elementInString(bandgapAlteringElements, numBandgapAlteringElements,crys->species+atom1*namelength);	
	//int iEle1 = crys_elementIndex(crys, crys->species+atom1*namelength);	
	int atom1neighbors[maxconnections];	
	int numatom1neighbors;
	cn_nearestNeighbors(cn,atom1neighbors,&numatom1neighbors,atom1);
	
	int iEle2 = crys_elementInString(bandgapAlteringElements, numBandgapAlteringElements,crys->species+atom2*namelength);	
	//int iEle2 = crys_elementIndex(crys, crys->species+atom2*namelength);	
	int atom2neighbors[maxconnections];
	int numatom2neighbors;
	cn_nearestNeighbors(cn,atom2neighbors,&numatom2neighbors,atom2);
	int i,j;
	//printf("iEle1 %d iEle2 %d\n",iEle1,iEle2);
	for(i=0;i<numatom1neighbors;i++)//subtract 1 from the counter for each coordination neighbor of atom 1 of element 1. add 1 from the counter for each coordination neighbor of atom 1 of element 2.
	{
		//printf("Checking if neighbor %d %s is a %s\n",i,crys->species+*(atom1neighbors+i)*namelength,coordElement);
		//This strcmp is inefficient. We should use ranges instead.
		//if(!strcmp(crys->species+*(atom1neighbors+i)*namelength,coordElement))
		if(*(atom1neighbors+i) >= estart && *(atom1neighbors+i) <estart+numatoms)
		{
			//printf("Here\n");
			//printf("Atom %d is a neighbor of atom 1, subtracting ele %d adding ele %d\n",*(atom1neighbors+i),iEle1,iEle2);
			if(iEle1 != -1)//Invoked when one atom is a vacancy
			(*(coord+(*(atom1neighbors+i)-estart)*numBandgapAlteringElements+iEle1)) -= forward_reverse;
			if(iEle2 != -1)
			(*(coord+(*(atom1neighbors+i)-estart)*numBandgapAlteringElements+iEle2)) += forward_reverse;
		}
	}
	for(i=0;i<numatom2neighbors;i++)//subtract 1 from the counter for each coordination neighbor of atom 2 of element 2. add 1 from the counter for each coordination neighbor of atom 2 of element 1
	{
		//This strcmp is inefficient. We should use ranges instead.
		//if(!strcmp(crys->species+*(atom2neighbors+i)*namelength,coordElement))
		if(*(atom2neighbors+i) >= estart && *(atom2neighbors+i) < estart+numatoms)
		{
			//printf("Here\n");
			//printf("Atom %d is a neighbor of atom 2, subtracting ele %d adding ele %d\n",*(atom2neighbors+i),iEle2,iEle1);
			if(iEle2 != -1)//Invoked when one atom is a vacancy
			(*(coord+(*(atom2neighbors+i)-estart)*numBandgapAlteringElements+iEle2)) -= forward_reverse;
			if(iEle1 != -1)//Invoked when one atom is a vacancy
			(*(coord+(*(atom2neighbors+i)-estart)*numBandgapAlteringElements+iEle1)) += forward_reverse;
		}
	}
	//OS_printCoord(data);
	
	//cn_printAdjacencyList(crys,cn);
}

//wrapper function. networkswap is typically much more EXPENSIVE THAN COORDSwap.
void bg_networkSwap(Configuration *config, int atom1, int atom2)
{
	OptoelectronicState *data = config->data;
	//printf("Network swap %d %d\n",atom1,atom2);
	if(data->bandgapNetwork != NULL)
	{
	bg_coordSwap(config,atom1,atom2,1);//network swap is always forward because the network is updated so performing a second forward swap will undo the first.
	cn_networkSwap(data->bandgapNetwork,atom1,atom2);
	}
}


double bg_energy(Configuration *config)
{
	//Calculates the energy of n band edge excitations
	crystal *crys = config->crys;
	OptoelectronicState *data = config->data;
    int numcoordinationatoms = crys_elementCount(crys,coordElement);
    double energy;

	if(data->bandgapNetwork == NULL)
	{
		data->bandgapNetwork = malloc(sizeof(crystalnetwork));
		data->bandgapNetwork->maxconnections=maxCoordination;
		//printf("Allocation blocks of %d and %d for bandgapNetwork\n",crys->totalAtoms,data->bandgapNetwork->maxconnections*crys->totalAtoms);
		data->bandgapNetwork->adjacencyList = calloc(data->bandgapNetwork->maxconnections*crys->totalAtoms,sizeof(int));
		data->bandgapNetwork->numAdjacent = calloc(crys->totalAtoms,sizeof(int));
		cn_collapseNetwork(crys, data->bandgapNetwork, coordElement, crys_appendElementString(bandgapAlteringElements,numBandgapAlteringElements,VX), numBandgapAlteringElements+1, data->coord, coordinationNumber, maxCoordination);
		//cn_printAdjacencyList(crys,data->bandgapNetwork);
		//printf("Memory locations %X %X %d\n",data->bandgapNetwork,data->coord,data->bandgapNetwork->maxconnections);
		cn_bitA_clusterApprox(crys, crys->network, coordElement, bandgapAlteringElements, numBandgapAlteringElements, data->coord, coordinationNumber, data->bandgapNetwork->maxconnections);
		//This call should be faster than the above, but it seems to fail.
		//cn_bitA_clusterApprox(crys, data->bandgapNetwork, coordElement, bandgapAlteringElements, numBandgapAlteringElements, data->coord, 1, data->bandgapNetwork->maxconnections);
		//OS_printCoord(data);
		//exit(0);
	}


	//cn_printAdjacencyList(crys,data->bandgapNetwork);
	//OS_printCoord(data);
	//interatomic_energy still uses the slow poor-scaling method to compute. We will need to update it for speed. 
	data->interatomicEnergy = 0;
	//data->interatomicEnergy = interatomic_energy(config);
	data->photocarrierEnergy = OS_energy(crys,data->coord,config);
	energy = data->interatomicEnergy + data->photocarrierEnergy;
	
	if(!isnormal(data->photocarrierEnergy))
	{
		printf("Energy is not normal\n");
		printf("Max connections %d\n",data->bandgapNetwork->maxconnections);
		OS_printCoord(data);
		//cn_printAdjacencyList(crys,data->bandgapNetwork);
		exit(0);
	}
	//printf("interatomic energy is %g, photocarrier Energy is %g, total Energy is %g saving it to %d\n",data->interatomicEnergy, data->photocarrierEnergy,energy,(config->energy));
	(config->energy) = energy;
    
	return energy;
}

void bg_registerSettings()
{
	lattikem_registerSettings();
	registerDouble(&numExcitations,"numExcitations",10);
	registerDouble(&numExcitationsSecondStep,"numExcitationsSecondStep",0);
	registerDouble(&fermiConvergence,"fermiConvergence",1e-5);
	registerDouble(&weightCutoff,"weightCutoff",0.1);
	registerInt(&coordinationNumber,"coordinationNumber",3);
	registerInt(&biasTurnOn,"biasTurnOn",0);
	registerInt(&biasSwitch,"biasSwitch",1000000);
	registerInt(&numOpticsEnergies,"numOpticsEnergies",1000);
	registerInt(&saveChargeDensity,"saveChargeDensity",0);
	registerInt(&numGangs,"numGangs",1);
	registerInt(&vectorLength,"vectorLength",32);
	registerInt(&BGWave,"BGWave",0);
	registerEnum(2,&thermalDistribution,"thermalDistribution",FERMI_DISTRIBUTION,"boltzmann","fermi");
	sourceControl_RegisterSettings();	
}

void bg_setup(double *newbandgapFunction(int numBandgapAlteringElements, int *numEachElement),char *newcoordElement,char *newbandgapAlteringElements,int newnumBandgapAlteringElements,int *coordinationShells,int newnumStates,
char *newrepulsiveElements,int newnumRepulsiveElements, double *newrepulsiveEnergies,int newrepulsiveShellSize, int newmaxRepulsiveCoordination,
void *traj_generator(Trajectory *traj))
{
	shoulder = malloc(numOpticsEnergies*sizeof(double));
    shoulderenergies = malloc(numOpticsEnergies*sizeof(double));
	absorptionshoulder(getTemperature(),0.025,shoulder,shoulderenergies,&shoulder_llimit,&shoulder_ulimit);
	bandgapFunction = newbandgapFunction;	
	coordElement = newcoordElement;
	bandgapAlteringElements = newbandgapAlteringElements;
	numBandgapAlteringElements = newnumBandgapAlteringElements;
	maxCoordination = coordinationShells[coordinationNumber-1];
	numStates = newnumStates;
	
	//This repulsive energy functionality should be pushed into a new module which could be loaded by the bandgap module.
	repulsiveElements = newrepulsiveElements;
	numRepulsiveElements = newnumRepulsiveElements;
	repulsiveEnergies = newrepulsiveEnergies;
	repulsiveShellSize = newrepulsiveShellSize;
	maxRepulsiveCoordination = newmaxRepulsiveCoordination;
	
	//This is currently set without a function for combining optoelectronic states. If you try to post-process it will crash
	//Combining nonweighted is now deprecated and weighted is the norm.
	LD_setup(bg_energy,OS_save,OS_allocate,OS_free,NULL,OS_combineWeighted,traj_generator,sizeof(OptoelectronicState),bg_networkSwap,bg_coordSwap);
}

void bg_trajectory(Trajectory *traj)
{
	int numSteps = getNumSteps();
	traj->numExternalConditions=1;//Number of carriers
    traj->externalConditions=malloc(1*numSteps * sizeof(double));
    
    if(BGWave)//Use SourceControl.c to set a periodic repeating excitation condition
    {
		sourceControl_Trajectory(traj);
	}
	else
	{	
    int iStep;
	for(iStep =0;iStep<numSteps;iStep++)
	if(iStep < biasTurnOn)
	*(traj->externalConditions+traj->numExternalConditions*iStep+0) = 0; //Initial step
	else
	if(iStep >= biasTurnOn && iStep < biasSwitch)
	*(traj->externalConditions+traj->numExternalConditions*iStep+0) = numExcitations; //Turn on Light
	else
	*(traj->externalConditions+traj->numExternalConditions*iStep+0) = numExcitationsSecondStep;//Switch to second intensity (e.g. turn off)
	}
	
	traj->numSeriesData=3;
	traj->seriesData=calloc(traj->numSeriesData*numSteps,sizeof(double));	
	traj->seriesDataNames = malloc(traj->numSeriesData*sizeof(char *));
	*(traj->seriesDataNames) = "Efermi";	
	*(traj->seriesDataNames+1) = "photocarrierEnergy";	
	*(traj->seriesDataNames+2) = "interatomicEnergy";	
			
}
