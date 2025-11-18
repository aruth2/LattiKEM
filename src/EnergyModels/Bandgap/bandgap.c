#ifdef __CUDACC__
#include "bandgap_GPU.h"
#else
#include "bandgap.h"
#endif

//Global settings for bandgap model
char *bandgapAlteringElements;	
int numBandgapAlteringElements;
char *coordElement;
int	coordinationNumber;//The number of shells used to define the local bandgap.
int maxCoordination;
int numStates;

double (*bandgapFunction)(int numBandgapAlteringElements, int *numEachElement);

//Spectra settings
int numOpticsEnergies;
double *shoulder, *shoulderenergies;
double shoulder_llimit, shoulder_ulimit;

//Thermodynamic Distribution
int thermalDistribution;
double fermiConvergence;
//For saving charge density to plot
int saveChargeDensity;
double weightCutoff;

//int numGangs;
//int vectorLength;

//Global settings for interatomic interactions.
int interatomicMode;  
char *repulsiveElements;
double **repulsiveEnergies;
int	numRepulsiveElements;
int	repulsiveCoordinationNumber;
int	maxRepulsiveCoordination;

//CW illumination
double	numExcitations;
double	numExcitationsSecondStep;
int biasTurnOn;
int biasSwitch;

//Pulsed illumination
int BGWave;

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
	OptoelectronicState *data = (OptoelectronicState *)config->data;
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
		int d;
		//This assumes the data (atom locations) are correctly oriented in first X, then Y, then Z
		for(i = 0; i < numStates; i++)
		{
			//This has to do with fixed character length floating point numbers. 
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
	OptoelectronicState *data = (OptoelectronicState *)outconfig->data;
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

#define maxUpdate 10000
void OS_allocate(Configuration *config)
{
	OptoelectronicState *OS = (OptoelectronicState *)config->data;
    OS->chargeDensity = (double *)malloc(numStates*sizeof(double));
    OS->emission = (double *)malloc(numOpticsEnergies*sizeof(double));
    OS->DOS = (double *)malloc(numOpticsEnergies*sizeof(double));
    OS->current = (double *)malloc(numOpticsEnergies*sizeof(double));
    OS->absorption = (double *)malloc(numOpticsEnergies*sizeof(double));
    OS->photonenergies = (double *)malloc(numOpticsEnergies*sizeof(double));
    
#ifdef __CUDACC__
	bgg_allocate(config,numStates);
#else	
    OS->energystates = (double *)malloc(numStates*sizeof(double));
    OS->weights = (double *)malloc(numStates*sizeof(double));
    OS->fermiEnergy = (double *)calloc(1,sizeof(double));//The first guess uses the initialized value.  
    OS->photocarrierEnergy = (double *)malloc(sizeof(double));
#endif    
 
 	OS->bandgapNetwork = NULL;
	OS->BGcoord = (int *)calloc(numStates*numBandgapAlteringElements,sizeof(int));
    
    OS->IANetwork = NULL;
	OS->IAcoord = (int *)calloc(numStates*numBandgapAlteringElements,sizeof(int));
    
    OS->BGUpdateList = (int *)calloc(maxUpdate,sizeof(int));
    
    //OS->fermiEnergy=0;
    //Scalars
	//This is somewhat redundant in that the number of scalars should match the number of series data. It should only be set once
	OS->numScalars=3;
	OS->scalars=(double **)calloc(OS->numScalars,sizeof(double *));
	*(OS->scalars) = OS->fermiEnergy;
	*(OS->scalars+1) = OS->photocarrierEnergy;
	*(OS->scalars+2) = &OS->interatomicEnergy;
}

void OS_free(Configuration *config)
{
	OptoelectronicState *OS = (OptoelectronicState *)config->data;
    free(OS->chargeDensity);
    free(OS->emission);
    free(OS->absorption);
    free(OS->photonenergies);
    free(OS->DOS);
    free(OS->current);
#ifdef __CUDACC__
	bgg_free(OS);
#else
    free(OS->energystates);
    free(OS->weights);
    free(OS->fermiEnergy);
	free(OS->photocarrierEnergy);
#endif
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
	for(int i = 0;i<numEnergyLevels;i++)
	{
		*(weights+i) = 1.0/(1.0+exp((*(energyLevels+i)-fermiEnergy)/temperature));
		*weightSum += *(weights+i);
		*energySum += *(energyLevels+i) * *(weights+i);
	}
}


double solveFermiEnergy(double currentExcitations, double fermiGuess, double fermiConvergence, double *energyLevels,double *weights, int numEnergyLevels, double temperature, double *energy)
{
	double fermiMin = fermiGuess - 7.5*fermiConvergence;
	double fermiMax = fermiGuess + 8.49*fermiConvergence;//This is slightly offset so that the fermi level could be almost unchanged. 
	//double fermiEnergy;
	double totalCarriersMin=0;
	double totalCarriersMax=0; 
	//Try to bound the fermi level within a small box around the previous result which requires not more than 6 iterations to converge.
	fermiSum( energyLevels, fermiMin, temperature, numEnergyLevels, weights, &totalCarriersMin, energy);
	fermiSum( energyLevels, fermiMax, temperature, numEnergyLevels, weights, &totalCarriersMax, energy);
	
	
	if((totalCarriersMin > currentExcitations) || (totalCarriersMax < currentExcitations))//If the solution is not within the bounds expand the bounds.  	
	{
		//printf("Fermi level outside bounds, Min %.10g max %.10g, min carriers %.10g, max carriers %.10g, target carriers %g\n",fermiMin,fermiMax,totalCarriersMin,totalCarriersMax,currentExcitations);
		fermiMin = optics_llimit;
		fermiMax = optics_ulimit;
	}
		
	double totalCarriers=0;
	while(fermiMax-fermiMin > fermiConvergence)
	{
		fermiGuess = (fermiMax+fermiMin)/2.0;	
		fermiSum( energyLevels, fermiGuess, temperature, numEnergyLevels, weights, &totalCarriers, energy);
		if(totalCarriers > currentExcitations) //Lower search range
			fermiMax = fermiGuess;
		else
			fermiMin = fermiGuess;
		//printf("Fermi range %g %g Carriers %g %g\n",fermiMin,fermiMax,totalCarriers,currentExcitations);
	}
	return fermiGuess;
}

void OS_energy(crystal *crys, int *BGcoord,Configuration *config)
{
	//This function calculates the energy of n band edge excitations based on local composition. 
    //First local bandgaps are identified and then excitons are distributed to preferentially reside
    //in low bandgap regions.
	//The probability is considered via boltzmann or fermi distribution
	
	
	int numcoordinationatoms = crys_elementCount(crys,coordElement);	
    double currentExcitations =  *(config->externalConditions+0);
	double temperature = getTemperature();
	
    OptoelectronicState *OS = (OptoelectronicState *)config->data;
    //printf("Data is located at %X\n",config->data);
	int i;
	
#ifdef __CUDACC__
	//nvtxRangePushA(":ELEVELS_IN_IND");
	//int copyElevels = !(OS->numBGUpdate);
	int copyElevels = 1;
#endif
	if((OS->numBGUpdate))//Recalcualte the bandgap for only those energy states that have been altered. 
	{
		for(int iUpdate=0;iUpdate<(OS->numBGUpdate);iUpdate++)
			{
				i = *(OS->BGUpdateList+iUpdate);
				*(OS->energystates+i) = bandgapFunction(numBandgapAlteringElements, BGcoord+numBandgapAlteringElements*i);
#ifdef __CUDACC__
				//cudaMemcpyAsync(OS->d_energyLevels+i, OS->energystates+i, sizeof(double), cudaMemcpyHostToDevice,0);//Copying one element at a time has been far slower than copying the block even if only a small number are copied.
#endif
			}
	}
	else
	for(i=0;i<numcoordinationatoms;i++)
	{
		*(OS->energystates+i) = bandgapFunction(numBandgapAlteringElements, BGcoord+numBandgapAlteringElements*i);
	}
	OS->numBGUpdate = 0;
	//printf("Starting thermal distribution\n");
	
	if(thermalDistribution == BOLTZMANN_DISTRIBUTION)
	{
		double probabilitysum = 0;
		for(i=0;i<numcoordinationatoms;i++)
		{
			//By using the middle energy we keep the weights close to 1 to reduce floating point errors.
			*(OS->weights+i) = exp(-(*(OS->energystates+i) - (optics_ulimit+optics_llimit)/2.0)/temperature);
			probabilitysum += *(OS->weights+i);
			//printf("Energy state %g weight %g probability sum %g\n",*(OS->energystates+i),*(OS->weights+i),probabilitysum);
		}
		
		//printf("Lowest energy state %g highest %g\n",lowestenergystate,highestenergystate);
		
		for(i=0;i<numcoordinationatoms;i++)
		{
			*(OS->weights+i) = *(OS->weights+i) * currentExcitations/probabilitysum;
			//energy += phaseenergy[i] * exp(-phaseenergy[i]/temperature) * numExcitations/probabilitysum;
			*(OS->photocarrierEnergy) += *(OS->energystates+i) * *(OS->weights+i);
		    //printf("Energy state %g weight %g probability sum %g\n",*(OS->energystates+i),*(OS->weights+i)/numExcitations,probabilitysum);
		}
	}
	
	if(thermalDistribution == FERMI_DISTRIBUTION)
	{	
#ifdef __CUDACC__
		//printf("Calculating Fermi Distribution on GPU\n");
		int copyWeights = config->savedata || saveChargeDensity;
		bgg_runKernel(OS, copyElevels, copyWeights, numStates, currentExcitations);
#else
		//printf("Calculating Fermi Distribution on CPU\n");
		*(OS->fermiEnergy) =  solveFermiEnergy(currentExcitations, *(OS->fermiEnergy), fermiConvergence, OS->energystates,OS->weights,numcoordinationatoms,temperature,OS->photocarrierEnergy);
#endif
	}
	
	//printf("Saving output data\n");
    //Only calculate optics if they will be saved. This helps a lot on kinetic steps
    if(config->savedata)
    {
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
    
	//return energy;
}

void OS_printCoord(OptoelectronicState *OS)
{
	//OptoelectronicState *OS = config->data;
	int *BGcoord = OS->BGcoord;
	int *IAcoord = OS->IAcoord;
	int iState, iEle, sum;
	for(iState=0;iState<numStates;iState++)
	{
		printf("BG State %d\n",iState);
		sum = 0;
		for(iEle=0;iEle<numBandgapAlteringElements;iEle++)
			{
			printf("%s\t%d\t\t",bandgapAlteringElements+iEle*namelength,*(BGcoord+iState*numBandgapAlteringElements+iEle));
			sum += *(BGcoord+iState*numBandgapAlteringElements+iEle);
			}
		printf("%d tot\n",sum);
		if(sum > OS->bandgapNetwork->maxconnections)
			exit(1);
		
		if(interatomicMode == IA_SEPARATE_NETWORK)
		{	
		printf("IA State %d\n",iState);
		sum = 0;
		for(iEle=0;iEle<numBandgapAlteringElements;iEle++)
			{
			printf("%s\t%d\t\t",bandgapAlteringElements+iEle*namelength,*(IAcoord+iState*numBandgapAlteringElements+iEle));
			sum += *(IAcoord+iState*numBandgapAlteringElements+iEle);
			}
		printf("%d tot\n",sum);
		}	
	}
}
/*
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
    int BGcoord[numcoordinationatoms*(numRepulsiveElements-1)];
    
    //printf("computing interatomic energy from %s with %d repulsive Elements with %d max shell\n",repulsiveElements,numRepulsiveElements,maxRepulsiveCoordination);
    for(iele1=0;iele1<numRepulsiveElements-1;iele1++)
    {
		//printf("Counting element %s\n",repulsiveElements+namelength*iele1);
		numEle1 = crys_elementCount(crys,repulsiveElements+namelength*iele1);
		numElementNeighbors = numRepulsiveElements-iele1-1;
		//printf("Computing coordination of element %s which has %d element neighbors\n",repulsiveElements+namelength*iele1,numElementNeighbors);
		//cn_coordination(crys,repulsiveElements+namelength*iele1,repulsiveElements+namelength*(iele1+1),numElementNeighbors,coord,1,maxRepulsiveCoordination);//Implicitly nearest neighbors (coordination #1)
		
		//CPU_clusterApprox(crys,repulsiveElements+namelength*iele1,repulsiveElements+namelength*(iele1+1),numElementNeighbors,coord,1,maxRepulsiveCoordination);//Implicitly nearest neighbors (coordination #1)
		cn_bitA_clusterApprox(crys,crys->network,repulsiveElements+namelength*iele1,repulsiveElements+namelength*(iele1+1),numElementNeighbors,BGcoord,1,maxRepulsiveCoordination);//Implicitly nearest neighbors (coordination #1)
		//printf("Finnished Computing coordination of element %s\n",repulsiveElements+namelength*iele1);
		for(iele2 = iele1+1;iele2<numRepulsiveElements;iele2++)
		{
			numEle2 = 0;
			//printf("Summing neighbors of element %s and %s\n",repulsiveElements+namelength*iele1,repulsiveElements+namelength*iele2);
			for(iAtom1=0;iAtom1<numEle1;iAtom1++)
				numEle2 += BGcoord[iAtom1*numElementNeighbors+(iele2-iele1-1)];
			//printf("Found %d neighbors of type %s around atom %s which has %d members. This is repulsive pair %d\n",numEle2,repulsiveElements+namelength*iele2,repulsiveElements+namelength*iele1,numEle1,repulsivePair);
			//This implictly assumes that repulsiveShellSize is the number of nearest neighbors for each atom type. However, if some shells are larger than others, this will not work and we will need separate shell sizes for different atom types.
			//printf("Adding %g to energy\n",*(repulsiveEnergies+repulsivePair)*numEle2/repulsiveCoordinationNumber);
			energy += *(repulsiveEnergies+repulsivePair)*numEle2/repulsiveCoordinationNumber;
			repulsivePair++;
		}
	}
	
	return energy;
}*/

double interatomic_energy(Configuration *config)
{
	if(interatomicMode == IA_NONE)
		return 0;
		
	//crystal *crys = config->crys;
	OptoelectronicState *data = (OptoelectronicState *)config->data;	
	int *coord;
	//int maxConnections;
	int iele1, iele2, iAtom;
	int numEle1, numEle2;
	int repulsivePair=0;
	double energy=0;
	if(interatomicMode == IA_BG_NETWORK)
	{
		coord = data->BGcoord;
		//maxConnections = maxCoordination;//Value set but never used
	}
	
	if(interatomicMode == IA_SEPARATE_NETWORK)
	{
		coord = data->IAcoord;
		//maxConnections = maxRepulsiveCoordination;//Value set but never used
	}
	
	for(iele1=0;iele1<numBandgapAlteringElements-1;iele1++)
    {

		for(iele2 = iele1+1;iele2<numBandgapAlteringElements;iele2++)
		{
			for(iAtom=0;iAtom<numStates;iAtom++)
			{
				numEle1 = coord[iAtom*numBandgapAlteringElements+iele1];
				numEle2 = coord[iAtom*numBandgapAlteringElements+iele2];
				//This is very close to U * x1 * x2, where x1 and x2 are the fraction of ele1 and ele2.
				//However, the version below mutes the impact of vacancies.
				energy += **(repulsiveEnergies+repulsivePair)*numEle1*numEle2/(numEle1+numEle2)/(numEle1+numEle2);
			}
			repulsivePair++;
		}
	}
	return energy;
		
}

void bg_coordSwap(crystal *crys, crystalnetwork *cn, int *coord, int atom1, int atom2, int forward_reverse, int *BGUpdateList, int *numBGUpdate)
{
	//1 for forward, -1 to go in reverse and undo a swap.
	
	//printf("Memory Locations in BG_coordswap %X %X %X\n",cn,coord,crys);
	
	//if(cn == NULL)
	//return;
	
	//printf("Swapping %d and %d\n",atom1,atom2);
	//OS_printCoord(data);	
	//cn_printAdjacencyList(crys,cn);
	//printf("Initializing variables\n");
	int maxconnections = cn->maxconnections;
	
	int estart = crys_elementOffset(crys,coordElement);
	int numatoms = crys_elementCount(crys,coordElement);	
	
	int iEle1 = crys_elementInString(bandgapAlteringElements, numBandgapAlteringElements,crys->species+atom1*namelength);	
		
	int atom1neighbors[maxconnections];	
	int numatom1neighbors;
	//printf("Locating Neighbors\n");
	cn_nearestNeighbors(cn,atom1neighbors,&numatom1neighbors,atom1);
	
	int iEle2 = crys_elementInString(bandgapAlteringElements, numBandgapAlteringElements,crys->species+atom2*namelength);		
	int atom2neighbors[maxconnections];
	int numatom2neighbors;
	cn_nearestNeighbors(cn,atom2neighbors,&numatom2neighbors,atom2);
	int i;
	
	//*numBGUpdate = 0;
	//printf("iEle1 %d iEle2 %d\n",iEle1,iEle2);
	for(i=0;i<numatom1neighbors;i++)//subtract 1 from the counter for each coordination neighbor of atom 1 of element 1. add 1 from the counter for each coordination neighbor of atom 1 of element 2.
	{
		//printf("Checking if neighbor %d %s is a %s\n",i,crys->species+*(atom1neighbors+i)*namelength,coordElement);
		if(*(atom1neighbors+i) >= estart && *(atom1neighbors+i) <estart+numatoms)
		{
			//printf("Atom %d is a neighbor of atom 1, subtracting ele %d adding ele %d\n",*(atom1neighbors+i),iEle1,iEle2);
			if(iEle1 != -1)//Invoked when one atom is a vacancy
			(*(coord+(*(atom1neighbors+i)-estart)*numBandgapAlteringElements+iEle1)) -= forward_reverse;
			if(iEle2 != -1)
			(*(coord+(*(atom1neighbors+i)-estart)*numBandgapAlteringElements+iEle2)) += forward_reverse;
			
			if(insertionSort(BGUpdateList,*numBGUpdate,*(atom1neighbors+i)-estart) != -1)//This insertion for tracking which BGs need updating appears to be more costly than the time saved.
				(*(numBGUpdate))++;
			//*(BGUpdateList + *numBGUpdate) = *(atom1neighbors+i)-estart;//This can double count, but it's faster building the list
			//(*(numBGUpdate))++;
		}
	}
	for(i=0;i<numatom2neighbors;i++)//subtract 1 from the counter for each coordination neighbor of atom 2 of element 2. add 1 from the counter for each coordination neighbor of atom 2 of element 1
	{
		if(*(atom2neighbors+i) >= estart && *(atom2neighbors+i) < estart+numatoms)
		{
			if(iEle2 != -1)//Invoked when one atom is a vacancy
			(*(coord+(*(atom2neighbors+i)-estart)*numBandgapAlteringElements+iEle2)) -= forward_reverse;
			if(iEle1 != -1)//Invoked when one atom is a vacancy
			(*(coord+(*(atom2neighbors+i)-estart)*numBandgapAlteringElements+iEle1)) += forward_reverse;
			
			
			if(insertionSort(BGUpdateList,*numBGUpdate,*(atom2neighbors+i)-estart) != -1)
				(*(numBGUpdate))++;
			//*(BGUpdateList + *numBGUpdate) = *(atom1neighbors+i)-estart;//This can double count, but it's faster building the list
			//(*(numBGUpdate))++;
				
		}
	}
	//OS_printCoord(data);
	
	//cn_printAdjacencyList(crys,cn);
}

//wrapper function. networkswap is typically much more EXPENSIVE THAN COORDSwap.
void bg_networkSwap(Configuration *config, int atom1, int atom2)
{
	OptoelectronicState *OS = (OptoelectronicState *)config->data;
	//printf("Network swap %d %d\n",atom1,atom2);
	//printf("Memory locations %X %X\n",OS->bandgapNetwork,OS->BGcoord);
	if(OS->bandgapNetwork != NULL)
	{
	bg_coordSwap(config->crys,OS->bandgapNetwork,OS->BGcoord,atom1,atom2,1,OS->BGUpdateList,&(OS->numBGUpdate));//network swap is always forward because the network is updated so performing a second forward swap will undo the first.
	cn_networkSwap(OS->bandgapNetwork,atom1,atom2);
	}
	if((OS->IANetwork != NULL) && (interatomicMode == IA_SEPARATE_NETWORK))//The second clause should be redundant. 
	{
	bg_coordSwap(config->crys,OS->IANetwork,OS->IAcoord,atom1,atom2,1,OS->BGUpdateList,&(OS->numBGUpdate));//network swap is always forward because the network is updated so performing a second forward swap will undo the first.
	cn_networkSwap(OS->IANetwork,atom1,atom2);	
	}
}

void bg_partialSwap(Configuration *config, int atom1, int atom2, int forward_reverse)
{
	OptoelectronicState *OS = (OptoelectronicState *)config->data;
	if(OS->bandgapNetwork != NULL)
	{
		bg_coordSwap(config->crys,OS->bandgapNetwork,OS->BGcoord,atom1,atom2,forward_reverse,OS->BGUpdateList,&(OS->numBGUpdate));
	}
	if((OS->IANetwork != NULL) && (interatomicMode == IA_SEPARATE_NETWORK))
	{
		bg_coordSwap(config->crys,OS->IANetwork,OS->IAcoord,atom1,atom2,forward_reverse,OS->BGUpdateList,&(OS->numBGUpdate));
	}
}


double bg_energy(Configuration *config)
{
	//Calculates the energy of n band edge excitations
	crystal *crys = config->crys;
	OptoelectronicState *OS = (OptoelectronicState *)config->data;
    int numcoordinationatoms = crys_elementCount(crys,coordElement);
    double energy;

	if(OS->bandgapNetwork == NULL)
	{
		OS->bandgapNetwork = (crystalnetwork *)malloc(sizeof(crystalnetwork));
		OS->bandgapNetwork->maxconnections=maxCoordination;
		//printf("Allocation blocks of %d and %d for bandgapNetwork\n",crys->totalAtoms,OS->bandgapNetwork->maxconnections*crys->totalAtoms);
		OS->bandgapNetwork->adjacencyList = (int *)calloc(OS->bandgapNetwork->maxconnections*crys->totalAtoms,sizeof(int));
		OS->bandgapNetwork->numAdjacent = (int *)calloc(crys->totalAtoms,sizeof(int));
		cn_collapseNetwork(crys, OS->bandgapNetwork, coordElement, crys_appendElementString(bandgapAlteringElements,numBandgapAlteringElements,VX), numBandgapAlteringElements+1, OS->BGcoord, coordinationNumber, maxCoordination);
		//printf("Memory locations %X %X %d\n",OS->bandgapNetwork,OS->BGcoord,OS->bandgapNetwork->maxconnections);
		cn_bitA_clusterApprox(crys, crys->network, coordElement, bandgapAlteringElements, numBandgapAlteringElements, OS->BGcoord, coordinationNumber, OS->bandgapNetwork->maxconnections);
		//OS_printCoord(data);
		OS->numBGUpdate = 0;
	}
	if((OS->IANetwork == NULL) && (interatomicMode == IA_SEPARATE_NETWORK) )
	{
		OS->IANetwork = (crystalnetwork *)malloc(sizeof(crystalnetwork));
		OS->IANetwork->maxconnections=maxRepulsiveCoordination;
		//printf("Allocation blocks of %d and %d for bandgapNetwork\n",crys->totalAtoms,data->bandgapNetwork->maxconnections*crys->totalAtoms);
		OS->IANetwork->adjacencyList = (int *)calloc(OS->IANetwork->maxconnections*crys->totalAtoms,sizeof(int));
		OS->IANetwork->numAdjacent = (int *)calloc(crys->totalAtoms,sizeof(int));
		//Currently, the same elements are used for repulsive energy as bandgap altering. Will replace this later.
		cn_collapseNetwork(crys, OS->IANetwork, coordElement, crys_appendElementString(bandgapAlteringElements,numBandgapAlteringElements,VX), numBandgapAlteringElements+1, OS->IAcoord, repulsiveCoordinationNumber, maxRepulsiveCoordination);
		//printf("Memory locations %X %X %d\n",data->bandgapNetwork,data->coord,data->bandgapNetwork->maxconnections);
		cn_bitA_clusterApprox(crys, crys->network, coordElement, bandgapAlteringElements, numBandgapAlteringElements, OS->IAcoord, repulsiveCoordinationNumber, OS->IANetwork->maxconnections);
		//OS_printCoord(data);
	}

	//cn_printAdjacencyList(crys,data->bandgapNetwork);
	//OS_printCoord(data);
	//data->interatomicEnergy = 0;
	OS_energy(crys,OS->BGcoord,config);
	OS->interatomicEnergy = interatomic_energy(config);
	energy = OS->interatomicEnergy + *(OS->photocarrierEnergy);
	
	if(!isnormal(*(OS->photocarrierEnergy)))
	{
		printf("Energy is not normal\n");
		printf("Max connections %d\n",OS->bandgapNetwork->maxconnections);
		OS_printCoord(OS);
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
	registerInt(&repulsiveCoordinationNumber,"repulsiveCoordinationNumber",1);
	registerInt(&biasTurnOn,"biasTurnOn",0);
	registerInt(&biasSwitch,"biasSwitch",1000000);
	registerInt(&numOpticsEnergies,"numOpticsEnergies",1000);
	registerInt(&saveChargeDensity,"saveChargeDensity",0);
	//registerInt(&numGangs,"numGangs",1);
	//registerInt(&vectorLength,"vectorLength",32);
	//registerInt(&BGWave,"BGWave",0);
	registerEnum(3,&BGWave, "BGWave",CONSTANT_WAVE,"CW","PULSED","CROSSOVER");
	registerEnum(2,&thermalDistribution,"thermalDistribution",FERMI_DISTRIBUTION,"boltzmann","fermi");
	registerEnum(3,&interatomicMode, "interatomicMode",IA_NONE,"NONE","IA_BG_NETWORK","IA_SEPARATE_NETWORK");
	sourceControl_RegisterSettings();
	
#ifdef __CUDACC__
	bgg_registerSettings();
#endif	
}

void bg_setup(double newbandgapFunction(int numBandgapAlteringElements, int *numEachElement),char *newcoordElement,char *newbandgapAlteringElements,int newnumBandgapAlteringElements,int *coordinationShells,int newnumStates,
char *newrepulsiveElements,int newnumRepulsiveElements, double **newrepulsiveEnergies,
void traj_generator(Trajectory *traj))
{
	shoulder = (double *)malloc(numOpticsEnergies*sizeof(double));
    shoulderenergies = (double *)malloc(numOpticsEnergies*sizeof(double));
	absorptionshoulder(getTemperature(),0.025,shoulder,shoulderenergies,&shoulder_llimit,&shoulder_ulimit);
	bandgapFunction = newbandgapFunction;	
	coordElement = newcoordElement;
	bandgapAlteringElements = newbandgapAlteringElements;
	numBandgapAlteringElements = newnumBandgapAlteringElements;
	maxCoordination = coordinationShells[coordinationNumber-1];
	numStates = newnumStates;
	
	//This repulsive energy functionality should be pushed into a new module which could be loaded by the bandgap module.
	//Currently the bandgap altering elements and the repulsive elements are one and the same.
	repulsiveElements = newrepulsiveElements;
	numRepulsiveElements = newnumRepulsiveElements;
	
	repulsiveEnergies = newrepulsiveEnergies;
	//repulsiveCoordinationNumber = newrepulsiveCoordinationNumber;
	maxRepulsiveCoordination = coordinationShells[repulsiveCoordinationNumber-1];
	//newmaxRepulsiveCoordination;

	
	//This is currently set without a function for combining optoelectronic states. If you try to post-process it will crash
	//Combining nonweighted is now deprecated and weighted is the norm.
	LD_setup(bg_energy,OS_save,OS_allocate,OS_free,NULL,OS_combineWeighted,traj_generator,sizeof(OptoelectronicState),bg_networkSwap,bg_partialSwap);
}

void bg_trajectory(Trajectory *traj)
{
	int numSteps = getNumSteps();
	traj->numExternalConditions=1;//Number of carriers
    traj->externalConditions=(double *)malloc(1*numSteps * sizeof(double));
    
    
    if(BGWave == CONSTANT_WAVE)
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
    if(BGWave == PULSED)//Use SourceControl.c to set a periodic repeating excitation condition
    {
		sourceControl_Trajectory(traj);
	}
	if(BGWave == CROSSOVER)//First use CW up to biasSwitch, then Use SourceControl.c to set a periodic repeating excitation condition
    {
		sourceControl_Trajectory(traj);
		int iStep;
		for(iStep =0;iStep<numSteps;iStep++)
			if(iStep < biasTurnOn)
				*(traj->externalConditions+traj->numExternalConditions*iStep+0) = 0; //Initial step
			else
				if(iStep >= biasTurnOn && iStep < biasSwitch)
					*(traj->externalConditions+traj->numExternalConditions*iStep+0) = numExcitations; //Turn on Light
		
	}
	
	
	traj->numSeriesData=3;
	traj->seriesData=(double *)calloc(traj->numSeriesData*numSteps,sizeof(double));	
	traj->seriesDataNames = (char **)malloc(traj->numSeriesData*sizeof(char *));
	*(traj->seriesDataNames) = "Efermi";	
	*(traj->seriesDataNames+1) = "photocarrierEnergy";	
	*(traj->seriesDataNames+2) = "interatomicEnergy";	
			
}

double getFermiConvergence()
{
	return fermiConvergence;
}
