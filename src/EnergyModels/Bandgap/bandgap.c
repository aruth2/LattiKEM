#include "bandgap.h"

char *bandgapAlteringElements;	
int numBandgapAlteringElements;
char *coordElement;
int	coordinationNumber;//The number of shells used to define the local bandgap.
int maxCoordination;
int numStates;
int numOpticsEnergies;

double	temperature;//In units of eV There could be seperate electronic and ionic temperatures
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

double OS_energy(crystal *crys,int *coord,Configuration *config)
{
	//This function calculates the energy of n band edge excitations based on local composition. 
    //First local bandgaps are identified and then excitons are distributed to preferentially reside
    //in low bandgap regions.
	//The probability is considered via boltzmann or fermi distribution
	
	int numcoordinationatoms = crys_elementCount(crys,coordElement);	
    double currentExcitations =  *(config->externalConditions+0);

    OptoelectronicState *OS = config->data;
    //printf("Dat is located at %d\n",(int *)dat);
	int i;

	double lowestenergystate = 10;
	double highestenergystate = 0;
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
		//double fermiConvergence = 1e-3;//This should be a setting
		double fermiMin = optics_llimit;
		double fermiMax = optics_ulimit;
		double totalCarriers=0;
		while(fermiMax-fermiMin > fermiConvergence)
		{
			OS->fermiEnergy = (fermiMax+fermiMin)/2.0;	
			totalCarriers=0;
			energy = 0;
			for(i=0;i<numcoordinationatoms;i++)
			{
				*(OS->weights+i) = 1.0/(1.0+exp((*(OS->energystates+i) - OS->fermiEnergy)/temperature));
				totalCarriers += *(OS->weights+i);
				energy += *(OS->energystates+i) * *(OS->weights+i);
				//printf("Energy state %g weight %g probability sum %g\n",*(OS->energystates+i),*(OS->weights+i),probabilitysum);
			}
			if(totalCarriers > currentExcitations) //Lower fermiEnergy 
				fermiMax = OS->fermiEnergy;
			else
				fermiMin = OS->fermiEnergy;
		}
		//if(currentExcitations > 0)
			//printf("Fermi energy converged to %g %g carriers expected %g calculated\n",OS->fermiEnergy,currentExcitations,totalCarriers);
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
		cn_coordination(crys,repulsiveElements+namelength*iele1,repulsiveElements+namelength*(iele1+1),numElementNeighbors,coord,1,maxRepulsiveCoordination);//Implicitly nearest neighbors (coordination #1)
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

double bg_energy(Configuration *config)
{
	//Calculates the energy of n band edge excitations
	//The excitation energy is estimated based on local bandgaps and varies from place to place
	crystal *crys = config->crys;
	OptoelectronicState *data = config->data;
	//int bins[maxCoordination*numBandgapAlteringElements];
    int numcoordinationatoms = crys_elementCount(crys,coordElement);
    double energy;
	
	//int atomsthiselement = crys_elementCount(crys,coordElement);
	//int i,binelement;
	//for(i=0;i<(maxCoordination+1)*numBandgapAlteringElements;i++)
	//*(bins+i) = 0;
	
	int coord[numcoordinationatoms*numBandgapAlteringElements];
	
	//crys_printAllAtoms(crys);
	//crys_printElements(crys);
    //printf("determining the coordination from %d elements which are %s using coordination number %d and max coordination %d\n",numBandgapAlteringElements,bandgapAlteringElements,coordinationNumber,maxCoordination);
	
	cn_coordination(crys,coordElement,bandgapAlteringElements,numBandgapAlteringElements,coord,coordinationNumber,maxCoordination);
	
	data->interatomicEnergy = interatomic_energy(config);
	data->photocarrierEnergy = OS_energy(crys,coord,config);
	energy = data->interatomicEnergy + data->photocarrierEnergy;

	//printf("interatomic energy is %g, photocarrier Energy is %g, total Energy is %g saving it to %d\n",data->interatomicEnergy, data->photocarrierEnergy,energy,(config->energy));
	(config->energy) = energy;
    
	return energy;
}

void bg_registerSettings()
{
	pmc_registerSettings();
	registerDouble(&temperature,"temperature",0.025);
	registerDouble(&numExcitations,"numExcitations",10);
	registerDouble(&numExcitationsSecondStep,"numExcitationsSecondStep",0);
	registerDouble(&fermiConvergence,"fermiConvergence",1e-5);
	registerDouble(&weightCutoff,"weightCutoff",0.1);
	registerInt(&coordinationNumber,"coordinationNumber",3);
	registerInt(&biasTurnOn,"biasTurnOn",100);
	registerInt(&biasSwitch,"biasSwitch",1000000);
	registerInt(&numOpticsEnergies,"numOpticsEnergies",1000);
	registerInt(&saveChargeDensity,"saveChargeDensity",0);
	registerEnum(2,&thermalDistribution,"thermalDistribution",BOLTZMANN_DISTRIBUTION,"boltzmann","fermi");	
}

void bg_setup(double *newbandgapFunction(int numBandgapAlteringElements, int *numEachElement),char *newcoordElement,char *newbandgapAlteringElements,int newnumBandgapAlteringElements,int *coordinationShells,int newnumStates,
char *newrepulsiveElements,int newnumRepulsiveElements, double *newrepulsiveEnergies,int newrepulsiveShellSize, int newmaxRepulsiveCoordination,
void *traj_generator(Trajectory *traj))
{
	shoulder = malloc(numOpticsEnergies*sizeof(double));
    shoulderenergies = malloc(numOpticsEnergies*sizeof(double));
	absorptionshoulder(temperature,0.025,shoulder,shoulderenergies,&shoulder_llimit,&shoulder_ulimit);
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
	LD_setup(bg_energy,OS_save,OS_allocate,OS_free,NULL,OS_combineWeighted,traj_generator,sizeof(OptoelectronicState));
}

void bg_trajectory(Trajectory *traj)
{
	int numSteps = getNumSteps();
	traj->numExternalConditions=1;//Number of carriers
    traj->externalConditions=malloc(1*numSteps * sizeof(double));
    	
    int iStep;
	for(iStep =0;iStep<numSteps;iStep++)
	if(iStep < biasTurnOn)
	*(traj->externalConditions+traj->numExternalConditions*iStep+0) = 0; //Initial step
	else
	if(iStep >= biasTurnOn && iStep < biasSwitch)
	*(traj->externalConditions+traj->numExternalConditions*iStep+0) = numExcitations; //Turn on Light
	else
	*(traj->externalConditions+traj->numExternalConditions*iStep+0) = numExcitationsSecondStep;//Switch to second intensity (e.g. turn off)
	
	
	traj->numSeriesData=3;
	traj->seriesData=calloc(traj->numSeriesData*numSteps,sizeof(double));	
	traj->seriesDataNames = malloc(traj->numSeriesData*sizeof(char *));
	*(traj->seriesDataNames) = "Efermi";	
	*(traj->seriesDataNames+1) = "photocarrierEnergy";	
	*(traj->seriesDataNames+2) = "interatomicEnergy";	
			
}
