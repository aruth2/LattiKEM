#include "stacd.h"

double carrierGenerationRate; //Carriers per unit area (sq ang) per second
double bandRecombinationRate; //angstroms cubed per second
double mobility; // angstrom squared per volt per second
double latticeConstant = 6.3620;//DOI:10.1021/acs.jpclett.5b01432
double thermalGenerationConstant;
double maxDampingParameter;
double convergenceCriterion;
int zCompression;
int lightWave;
int stepsPerLightSwitch;
int maxElectronicIterations;
int lightStepDelay;
DeviceStack *architecture;

void CP_allocate(CarrierProfile *CP)
{
	int numzValues = architecture->numzValues;
	int numIons = architecture->numIons;
	int numStates = architecture->numStates;
	int numIonTypes = architecture->numIonTypes;

	//Scalars
	//This is somewhat redundant in that the number of scalars should match the number of series data. It should only be set once
	CP->numScalars=10;
	CP->scalars=calloc(CP->numScalars,sizeof(double *));
	*(CP->scalars) = &CP->polarization;
	*(CP->scalars+1) = &CP->chargeInjection;
	*(CP->scalars+2) = &CP->ionFieldEnergy;
	*(CP->scalars+3) = &CP->ionIonEnergy;
	*(CP->scalars+4) = &CP->electronEnergy;
	*(CP->scalars+5) = &CP->holeEnergy;
	*(CP->scalars+6) = &CP->totalEnergy;
	*(CP->scalars+7) = &CP->quasiFermiLevelSplitting;
	*(CP->scalars+8) = &CP->potentialDifference;
	*(CP->scalars+9) = &CP->carrierGeneration;
	
	//States
	CP->stateCharge = calloc(numStates,sizeof(double));
	CP->quasiFermiLevels = calloc(numStates,sizeof(double));
	CP->stateCurrent = calloc(numStates,sizeof(double));
	CP->stateCurrent2 = calloc(numStates,sizeof(double));
	CP->stateRecombination = calloc(numStates,sizeof(double));
	CP->stateChange = calloc(numStates,sizeof(double));
	
	//Layers
	CP->potential = calloc(numzValues,sizeof(double));
	CP->field = calloc(numzValues,sizeof(double));
	CP->layerCharge = calloc(numzValues,sizeof(double));
	CP->ionicCharge = calloc(numzValues,sizeof(double));
	CP->electronCharge = calloc(numzValues,sizeof(double));
	CP->holeCharge = calloc(numzValues,sizeof(double));
	CP->layerCurrent = calloc(2*numzValues,sizeof(double));
	CP->layerCurrent2 = calloc(2*numzValues,sizeof(double));
	CP->layerBandCharge = calloc(2*numzValues,sizeof(double));
	CP->layerRecombination = calloc(2*numzValues,sizeof(double));
	CP->generation = calloc(numzValues,sizeof(double));
	CP->layerChange = calloc(2*numzValues,sizeof(double));
	CP->ionsPerLayer = calloc(numIonTypes*numzValues,sizeof(double));
	CP->zValues = calloc(numzValues,sizeof(double));
		
	CP->numLayerVectors=20+numIonTypes;
	CP->layerVectors=calloc(CP->numLayerVectors,sizeof(double *));
	*(CP->layerVectors+0) = CP->zValues;
	*(CP->layerVectors+1) = CP->potential;
	*(CP->layerVectors+2) = CP->field;
	*(CP->layerVectors+3) = CP->layerCharge;
	*(CP->layerVectors+4) = CP->ionicCharge;
	*(CP->layerVectors+5) = CP->electronCharge;
	*(CP->layerVectors+6) = CP->holeCharge;
	*(CP->layerVectors+7) = CP->generation;
	*(CP->layerVectors+8) = CP->layerCurrent;
	*(CP->layerVectors+9) = CP->layerCurrent+numzValues;
	*(CP->layerVectors+10) = CP->layerCurrent2;
	*(CP->layerVectors+11) = CP->layerCurrent2+numzValues;
	*(CP->layerVectors+12) = CP->layerRecombination;
	*(CP->layerVectors+13) = CP->layerRecombination+numzValues;
	*(CP->layerVectors+14) = CP->layerChange;
	*(CP->layerVectors+15) = CP->layerChange+numzValues;
	*(CP->layerVectors+16) = CP->quasiFermiLevels;
	*(CP->layerVectors+17) = CP->quasiFermiLevels+numzValues;
	*(CP->layerVectors+18) = CP->layerBandCharge;
	*(CP->layerVectors+19) = CP->layerBandCharge+numzValues;
	//Defects
	CP->effectiveCharge = calloc(numIons,sizeof(double));
	
	int iIon,iz;
	for(iIon=0;iIon<numIonTypes;iIon++)
	*(CP->layerVectors+20+iIon) = CP->ionsPerLayer+numzValues*iIon;
	//This is only done to aid in saving
	for(iz=0;iz<numzValues;iz++)
		*(CP->zValues+iz) = *(architecture->zValues+iz);
	
	CP_initialize(CP);
}

void stacd_stateToIndex(crystal *crys,int *iState, int *iBand, int *iCharge, int *iz, int *iMaterial, int *iIon, int *iEle, int *iAtom, int *iDL, int *ionDL, int *iDS)
{
	/* Given the state index return all other indicies.
	 * IF this is a band state, iIon, iEle, iAtom, iDL, ionDL, and iDS are all set to -1
	 * If this is an ion state, IBand is set to the zIndex of the corresponding band state
	 * 
	 * Description of indicies
	 * 
	 * iState - which state in the state arrays this describes
	 * iBand - electron and hole bands associated with specific z values
	 * iCharge - 0 if electron, 1 if hole
	 * iz - the layer of the material
	 * iMaterial - which material the layer belongs to
	 * iIon enumerates the ion. Forms the index of ionEles and ionAtoms
	 * iEle - which element the ion is from amongst the defects
	 * iAtom - the index of the ion within the crystal
	 * iDL - which defect level is being described by the state
	 * ionDL - starts are 0 for each element and lists the offset for that element
	 * iDS - which defect state is being described. The state is the defect level + occupancy/locality
	 * */
	//printf("iState %d, iBand %d, iCharge %d, iz %d, iMaterial %d, iIon %d, iEle %d, iAtom %d, iDL %d, ionDL %d, iDS %d\n",
	// *iState, *iBand, *iCharge, *iz, *iMaterial, *iIon, *iEle, *iAtom, *iDL, *ionDL, *iDS); 
	 
	int numzValues = architecture->numzValues;
	int numIons = architecture->numIons;
	char *ionNames = architecture->ionNames;
	int *numIonDefectLevels = architecture->numIonDefectLevels;
	double *defectLevelCharge = architecture->defectLevelCharge;
	
	if(*iState<2*numzValues)
	{
		*iBand = *iState;
		*iz = *iBand % numzValues;
		*iCharge = *iBand / numzValues;
		*iIon = *iEle = *iAtom = *iDL = *ionDL = *iDS = -1;
	}
	else
	{
		int stateOffset = *iDS = *iState - 2*numzValues;//Find the correct element
		//With how often this is going to be called, this loop could really slow things down
		int eleCount = crys_elementCount(crys,ionNames);
		int nIDL=*(numIonDefectLevels);
		for(*iEle=0, *iIon=0, *iDL=0; stateOffset - nIDL * eleCount>=0; )
		{
		//printf("eleCount %d, nIDL %d  name %s iEle %d\n",eleCount, nIDL,ionNames + namelength* *iEle,*iEle);	
		
		if(*iEle == architecture->numIonTypes)
		printf("Went past the last element. Num ion Types %d\n",architecture->numIonTypes);
		stateOffset -= nIDL * eleCount;
		*iIon += crys_elementCount(crys,ionNames + namelength * *iEle);
		*iDL += *(numIonDefectLevels + *iEle);
		
		(*iEle)++;
		eleCount = crys_elementCount(crys,ionNames + namelength * *iEle);
		nIDL=*(numIonDefectLevels+*iEle);
		}
		
		*iIon += stateOffset / *(numIonDefectLevels + *iEle);
		*iDL += stateOffset %  *(numIonDefectLevels + *iEle);
		*ionDL = stateOffset %  *(numIonDefectLevels + *iEle);
		
		if( *(defectLevelCharge + *iEle) < 0)
		*iCharge = 0;
		else
		*iCharge = 1;
		
		*iAtom = *(architecture->ionAtoms + *iIon);
		*iz = zIndex(architecture,*(crys->positions+3* *iAtom+2));
		*iBand = *iz + *iCharge * numzValues;
	}
	*iMaterial = zMaterial(architecture,*iz);
	//printf("iState %d, iBand %d, iCharge %d, iz %d, iMaterial %d, iIon %d, iEle %d, iAtom %d, iDL %d, ionDL %d, iDS %d\n",
	//*iState, *iBand, *iCharge, *iz, *iMaterial, *iIon, *iEle, *iAtom, *iDL, *ionDL, *iDS);
}

void stacd_iDSToState(crystal *crys,int *iState, int *iBand, int *iCharge, int *iz, int *iMaterial, int *iIon, int *iEle, int *iAtom, int *iDL, int *ionDL, int *iDS)
{
	/* Given  all other indicies return the state index.
	 * IF this is a band state, iIon, iEle, and iAtom are all set to -1
	 * If this is an ion state, IBand is set to the zIndex of the corresponding band state
	 * */
	int numzValues = architecture->numzValues;
	*iState = 2*numzValues + *iDS;
	if(*iDS >= architecture->numDefectStates)
	printf("Starting from iDS %d\n",*iDS);
	stacd_stateToIndex(crys,iState, iBand, iCharge, iz, iMaterial, iIon, iEle, iAtom, iDL, ionDL, iDS);	
	
}
void stacd_iIonToState(crystal *crys,int *iState, int *iBand, int *iCharge, int *iz, int *iMaterial, int *iIon, int *iEle, int *iAtom, int *iDL, int *ionDL, int *iDS)
{
	/* Given  all other indicies return the state index.
	 * IF this is a band state, iIon, iEle, and iAtom are all set to -1
	 * If this is an ion state, IBand is set to the zIndex of the corresponding band state
	 * */
	char *ionNames = architecture->ionNames;
	int *numIonDefectLevels = architecture->numIonDefectLevels;
	 
	*iEle = *(architecture->ionEles + *iIon);
	int iEleIndex, ionRemainder = *iIon;
	for(iEleIndex = 0,*iDS = 0;iEleIndex < *iEle; iEleIndex++)
		{
			*iDS += *(numIonDefectLevels + iEleIndex) * crys_elementCount(crys, ionNames + namelength * iEleIndex);
			ionRemainder -= crys_elementCount(crys, ionNames + namelength * iEleIndex);
		}
	*iDS += ionRemainder * *(numIonDefectLevels + *iEle);
	stacd_iDSToState(crys,iState, iBand, iCharge, iz, iMaterial, iIon, iEle, iAtom, iDL, ionDL, iDS);	
	
}
void stacd_iBandToState(crystal *crys,int *iState, int *iBand, int *iCharge, int *iz, int *iMaterial, int *iIon, int *iEle, int *iAtom, int *iDL, int *ionDL, int *iDS)
{	
	*iState = *iBand;
	stacd_stateToIndex(crys,iState, iBand, iCharge, iz, iMaterial, iIon, iEle, iAtom, iDL, ionDL, iDS);
}

void CP_initialize(CarrierProfile *CP)
{
	//CP->dampingParameter = dampingParameter;
	CP->dampingParameter = 1e-3*maxDampingParameter;
	
	int numzValues = architecture->numzValues;
	int numMaterials = architecture->numMaterials;
	double *bandOffsets = architecture->bandOffsets;
	double *quasiFermiLevels = CP->quasiFermiLevels;
	double leftSideFermi = *bandOffsets;
	double rightSideFermi = *(bandOffsets + 2*numMaterials-1);
	
	for(int iz = 0;iz<numzValues;iz++)
	for(int iCharge = 0;iCharge < 2;iCharge++)
		*(quasiFermiLevels+ iz + iCharge*numzValues) = leftSideFermi + (rightSideFermi-leftSideFermi) * iz/numzValues;
	return;
}

void CP_free(CarrierProfile *CP)
{
	free(CP->scalars);
	free(CP->stateCharge);
	free(CP->quasiFermiLevels);
	free(CP->stateRecombination);
	free(CP->stateCurrent);
	free(CP->stateCurrent2);
	free(CP->stateChange);
	
	free(CP->zValues);
	free(CP->potential);
	free(CP->field);
	free(CP->layerCharge);
	free(CP->ionicCharge);
	free(CP->electronCharge);
	free(CP->holeCharge);
	free(CP->layerCurrent);
	free(CP->layerCurrent2);
	free(CP->layerBandCharge);
	free(CP->layerRecombination);
	free(CP->generation);
	free(CP->layerChange);
	free(CP->quasiFermiLevels);
	free(CP->ionsPerLayer);
	
	free(CP->layerVectors);
	
	free(CP->effectiveCharge);
	
}

void CP_combineWeighted(Configuration *configs, int numCombine, Configuration *outconfig, double *weights)
{
	//printf("Combining %d carrier profiles into one\n",numCombine);
	CarrierProfile *data = outconfig->data;
	CarrierProfile *cp;
	int numScalars = data->numScalars;
	int numzValues = architecture->numzValues;
	int numLayerVectors = data->numLayerVectors;
	
	int iz, iVector, iScalar, iConfig;
	double sumWeight = 0;

	for(iConfig = 0; iConfig < numCombine; iConfig++)
		sumWeight += *(weights + iConfig);
	for(iScalar = 0; iScalar < numScalars; iScalar++)
		*(*(data->scalars + iScalar)) = 0;

	for(iVector = 0; iVector < numLayerVectors; iVector++)
	for(iz = 0; iz < numzValues; iz++)
		*(*(data->layerVectors + iVector) + iz) = 0;
	
	for(iScalar = 0; iScalar < numScalars; iScalar++)
	{
		for(iConfig = 0; iConfig < numCombine; iConfig++)
		{
			cp = (CarrierProfile *)((configs + iConfig)->data);
			*(*(data->scalars + iScalar)) += *(weights + iConfig)* *(*(cp->scalars + iScalar));
			
		}
		//if(iScalar == 8)
		//{
		//	printf("Sum of potential differences is %g and dividing by sumWeight %g\n",*(*(data->scalars + iScalar)),sumWeight);
		//}
		*(*(data->scalars + iScalar)) /= sumWeight;
		*(outconfig->seriesData + iScalar) = *(*(data->scalars + iScalar));
	}
	for(iVector = 0; iVector < numLayerVectors; iVector++)
	for(iz = 0; iz < numzValues;iz++)
	{
		for(iConfig = 0;iConfig < numCombine; iConfig++)
		{
			cp = (CarrierProfile *)((configs + iConfig)->data);
			*(*(data->layerVectors + iVector) + iz) += *(weights + iConfig)* *(*(cp->layerVectors + iVector) + iz);
		}
		*(*(data->layerVectors + iVector) + iz) /= sumWeight;
	}
	
}


void CP_save(Configuration *config)
{
	FILE *outfile = fopen(config->dataFileName,"w");
	//printf("Saving a configuration to file %s\n",config->dataFileName);
	CarrierProfile *data = config->data;
	int numLayerVectors = data->numLayerVectors;
	int iz,iVector,iStride;
	double sum;
	for(iz = 0; iz < architecture->numzValues / saveStride; iz++)
	{
		for(iVector = 0; iVector < numLayerVectors; iVector++)
		{
			sum = 0;
			for(iStride = 0; iStride < saveStride; iStride++)
				sum += *(*(data->layerVectors + iVector) + saveStride * iz + iStride);
			fprintf(outfile,"%g ",sum / saveStride);
		}
		fprintf(outfile,"\n");
	}
	fclose(outfile);
}

void stacd_Trajectory(Trajectory *traj)
{
	int numSteps = getNumSteps();
    traj->numExternalConditions=2;//first is voltage/source control, second is light intensity
    traj->externalConditions=malloc(2*numSteps * sizeof(double));

	traj->numSeriesData=10;
	traj->seriesData=calloc(traj->numSeriesData*numSteps,sizeof(double));	
	traj->seriesDataNames = malloc(traj->numSeriesData*sizeof(char *));
	*(traj->seriesDataNames) = "pol";
	*(traj->seriesDataNames+1) = "chargeInjection";
	*(traj->seriesDataNames+2) = "ionFieldEnergy";
	*(traj->seriesDataNames+3) = "ionIonEnergy";
	*(traj->seriesDataNames+4) = "electronEnergy";
	*(traj->seriesDataNames+5) = "holeEnergy";
	*(traj->seriesDataNames+6) = "totalEnergy";	
	*(traj->seriesDataNames+7) = "quasiFermiLevelSplitting";	
	*(traj->seriesDataNames+8) = "potentialDifference";	
	*(traj->seriesDataNames+9) = "carrierGeneration";	
	
	sourceControl_Trajectory(traj);
	traj_wave(traj,stepsPerLightSwitch,1,lightWave,0,carrierGenerationRate,lightStepDelay);
}

void stacd_registerSettings()
{
	registerDouble(&carrierGenerationRate,"carrierGenerationRate",0);
	registerDouble(&bandRecombinationRate,"bandRecombinationRate",1e15);
	registerDouble(&mobility,"mobility",1e15);	
	registerDouble(&maxDampingParameter,"maxDampingParameter",4e-2);	
	registerDouble(&convergenceCriterion,"convergenceCriterion",1e-3);	
	registerInt(&maxElectronicIterations,"maxElectronicIterations",10000);
	registerInt(&zCompression,"zCompression",1);
	registerInt(&lightStepDelay,"lightStepDelay",0);
	registerInt(&stepsPerLightSwitch,"stepsPerLightSwitch",100);
	
	registerEnum(3,&lightWave,"lightWave",SQUARE_WAVE,"square","triangle","sine");
	
	sourceControl_RegisterSettings();
}

void stacd_setup(DeviceStack *newArchitecture, void (*traj_generator)(Trajectory *))
{
	architecture = newArchitecture;
	thermalGenerationConstant = 0.00635 * pow(getTemperature(),1.5);
	printf("%d z values, %d ions, of %d charged elements having %d states\n",architecture->numzValues,architecture->numIons,architecture->numIonTypes,architecture->numStates);
	LD_setup(stacd_Energy,CP_save,CP_allocate,CP_free,NULL, CP_combineWeighted, traj_generator, sizeof(CarrierProfile) );

}

void printDeviceStack(DeviceStack *stack)
{
	int iMaterial;
	for(iMaterial=0;iMaterial<stack->numMaterials;iMaterial++)
	{
		printf("Material %s occupies %g to %g ang which is indicies %d to %d, has a dielectric constant of %g, band offsets of %g and %g\n"
		,*(stack->materialNames+iMaterial),*(stack->zValues+*(stack->materialTransitions+iMaterial))
		,iMaterial == stack->numMaterials-1 ? *(stack->zValues+stack->numzValues-1) : *(stack->zValues+*(stack->materialTransitions+iMaterial+1))
		,*(stack->materialTransitions+iMaterial)
		,iMaterial == stack->numMaterials-1 ? stack->numzValues-1 : *(stack->materialTransitions+iMaterial+1)-1
		, *(stack->dielectricConstants+iMaterial), *(stack->bandOffsets+2*iMaterial), *(stack->bandOffsets+2*iMaterial+1)
		); 
	}
	int iBand;
	for(iBand=0;iBand<stack->numzValues;iBand++)
	{
		printf("Layer %d is at height %g is part of material %d it has %g formula units\n",iBand,*(stack->zValues+iBand),zMaterial(stack,iBand),*(stack->formulaUnits+iBand));
	}
}

double bandFermiLevelToOccupancy(double fermiLevel,double stateVolume, int iCharge)
{
	//Assuming a parabolic band with a unitary effective mass, this calculates the Fermi level of the band relative to its CBM/VBM 
	//Occupancy is assumed to be expressed in charges per cubic angstrom
	
	double temp = getTemperature();
	
	double occupancy = thermalGenerationConstant * stateVolume * exp(fermiLevel * pow(-1,iCharge) /temp);
	if(occupancy <= minorityDensityCutoff)
	{
	occupancy = minorityDensityCutoff;
	printf("band fermi level of %g results in occupancy %g\n",fermiLevel,occupancy);
	}
	//if(occupancy >= 1-majorityDensityCutoff)
	//{
	//occupancy = 1-majorityDensityCutoff;
	//printf("fermi level of %g results in occupancy %g\n",fermiLevel,occupancy);
	//}
	//This check is wasting a lot of cycles, let's see if we can remove it.
	if(!isnormal(occupancy))
	{
	printf("Band fermi level %g for occupancy %g is not normal\n",fermiLevel,occupancy);
	exit(0);
	}
	//printf("Band fermi level %g for occupancy %g\n",value,occupancy);
	return occupancy;
}

double defectFermiLevelToOccupancy(double fermiLevel, int iCharge)
{
	double temperature = getTemperature();
	double occupancy = 1.0/exp((-pow(-1,iCharge) * fermiLevel/temperature)+1);
	if(occupancy <= minorityDensityCutoff)
	{
	occupancy = minorityDensityCutoff;
	printf("defect fermi level of %g results in occupancy %g\n",fermiLevel,occupancy);
	}
	/*if(occupancy >= 1-majorityDensityCutoff)
	{
	occupancy = 1-majorityDensityCutoff;
	//printf("fermi level of %g results in occupancy %g\n",fermiLevel,occupancy);
	}*/
	
	if(!isnormal(occupancy))
	printf("Defect fermi level %g for occupancy %g is not normal\n",fermiLevel,occupancy);
	return occupancy;
}

double bandDefectFermiLevelToOccupancy(double fermiLevel, double defectDensity, double volume, double temperature, double *defectLevels, int numDefectLevels, int iCharge)
{
	int iDL;
	double occupancy= thermalGenerationConstant * volume * exp(fermiLevel * pow(-1,iCharge) /temperature);
	
	for(iDL=0;iDL<numDefectLevels;iDL++)
		occupancy += defectDensity*volume / (exp(- pow(-1,iCharge)*(fermiLevel-*(defectLevels+iDL))/temperature)+1);
	
	if(occupancy <= minorityDensityCutoff)
	{
	occupancy = minorityDensityCutoff;
	printf(" band defect fermi level of %g results in occupancy %g\n",fermiLevel,occupancy);
	}
	
	//if(occupancy >= 1-majorityDensityCutoff)
	//{
	//occupancy = 1-majorityDensityCutoff;
	//printf("fermi level of %g results in occupancy %g\n",fermiLevel,occupancy);
	//}
	
	//printf("Converged in %d iterations\n",numIterations);
	return occupancy;
}

int zIndex(DeviceStack *architecture, double z)
{
	//Given the location of an object, this returns the corresponding z index of that object
	//This could be sped up by a lookup table
	return binarysearchDouble(architecture->zValues,z,0,architecture->numzValues)-1;
}

double zPotential(DeviceStack *architecture, CarrierProfile *CP, double z)
{
	/* This is for finding the potential of Ions in STACD. Due to the z Coarsening,
	 * the potential is quite coarse and stepped between the band sites. These coarse
	 * steps do not influence ion migration unless the ions are at the edge of the step.
	 * Therefore, the potential is smoothed for the ions by linear interpolation.
	 * The z index corresponding to the band the ion is within is used as the center, and the
	 * potentials above and below are also used.
	 * */
	 int iz = zIndex(architecture,z);
	 int numzValues= architecture->numzValues;
	 double *potential = CP->potential;
	 double *zValues = architecture->zValues;
	 double *thickness = architecture->thickness;;
	 
	 int left = z-*(zValues+iz)-*(thickness+iz)/2 < 0;
	 
	 if(iz==numzValues-1 || iz == 0)
	 return *(CP->potential+iz);
	 else
	 if(left)
	 return interpolate(z,*(zValues+iz-1)+*(thickness+iz-1)/2,*(zValues+iz)+*(thickness+iz)/2,*(potential+iz-1),*(potential+iz));
	 else
	 return interpolate(z,*(zValues+iz)+*(thickness+iz)/2,*(zValues+iz+1)+*(thickness+iz+1)/2,*(potential+iz),*(potential+iz+1));
}

double zFermi(DeviceStack *architecture, CarrierProfile *CP, double z, int iCharge)
{
	/* This is for finding the potential of Ions in STACD. Due to the z Coarsening,
	 * the potential is quite coarse and stepped between the band sites. These coarse
	 * steps do not influence ion migration unless the ions are at the edge of the step.
	 * Therefore, the potential is smoothed for the ions by linear interpolation.
	 * The z index corresponding to the band the ion is within is used as the center, and the
	 * potentials above and below are also used.
	 * */
	 int iz = zIndex(architecture,z);
	 int numzValues= architecture->numzValues;
	 double *potential = CP->potential;
	 double *zValues = architecture->zValues;
	 double *thickness = architecture->thickness;;
	 
	 int left = z-*(zValues+iz)-*(thickness+iz)/2 < 0;
	 
	 if(iz==numzValues-1 || iz == 0)
	 return *(CP->quasiFermiLevels+iz+iCharge*numzValues);
	 else
	 if(left)
	 return interpolate(z,*(zValues+iz-1)+*(thickness+iz-1)/2,*(zValues+iz)+*(thickness+iz)/2,*(CP->quasiFermiLevels+iz+iCharge*numzValues-1),*(CP->quasiFermiLevels+iz+iCharge*numzValues));
	 else
	 return interpolate(z,*(zValues+iz)+*(thickness+iz)/2,*(zValues+iz+1)+*(thickness+iz+1)/2,*(CP->quasiFermiLevels+iz+iCharge*numzValues),*(CP->quasiFermiLevels+iz+iCharge*numzValues+1));
}

int zMaterial(DeviceStack *architecture, int iz)
{
	/* Converts from a z index to which material it exists in
	 * */
	int iMaterial;
	//printf("Searching for what material %d belongs to\n",iz);
	iMaterial = binarysearch(architecture->materialTransitions,iz,architecture->numMaterials);
	//This corrects for the way that binarySearch reports the index which differentiates between exact matches and values lieing between the layers
	if(iMaterial < 0)
	iMaterial = -iMaterial - 1;
	else
	iMaterial--;
	if(iMaterial == architecture->numMaterials)
	iMaterial = architecture->numMaterials -1;

	return iMaterial;
}

//Create a lookup table which connects iIon to the corresponding atom number in the crystal
void fillIonAtoms(DeviceStack *architecture)
{
	int iEle;
	int eleStart, eleEnd;
	int iAtom,iIon = 0;	
	int *ionAtoms = architecture->ionAtoms;
	int *ionEles = architecture->ionEles;
	crystal *crys = architecture->initializationCrystal;
	for(iEle = 0; iEle<architecture->numIonTypes;iEle++)
	{
	eleStart = crys_elementOffset(crys,architecture->ionNames+iEle*namelength);
	eleEnd = eleStart + crys_elementCount(crys,architecture->ionNames+iEle*namelength);
	printf("Element %d %s starts at %d and ends at %d\n",iEle,architecture->ionNames+iEle*namelength,eleStart,eleEnd);
		for(iAtom=eleStart;iAtom<eleEnd;iAtom++)
		{
			*(architecture->ionAtoms+iIon) = iAtom;
			*(architecture->ionEles+iIon) = iEle;
			iIon++;
		}
	}
}

void compressZProfile(DeviceStack *architecture)
{
	/* This coarsens the z grid reducing the number of band states.
	 * */
	double *newzValues = calloc(architecture->numzValues,sizeof(double));
	double *newFormulaUnits = calloc(architecture->numzValues,sizeof(double));
	double *newThickness = calloc(architecture->numzValues,sizeof(double));
	int *newMaterialTransitions = calloc(architecture->numMaterials,sizeof(int));
	int *materialTransitions = architecture->materialTransitions;
	double *oldFormulaUnits = architecture->formulaUnits;
	double *oldThickness = architecture->thickness;
	double *oldzValues = architecture->zValues;
	double numAbsorberLayers = 0;
	int oldnumzValues=architecture->numzValues;
	int oldiz,newiz,compressionIndex;
	int iMaterial,previousMaterial;
	
	*newzValues = *oldzValues;
	*newMaterialTransitions = 0;
	previousMaterial = 0;
	for(oldiz=0,newiz=0,compressionIndex=0;oldiz<oldnumzValues;oldiz++)
	{
		iMaterial = zMaterial(architecture,oldiz);
		//At the start of a new material always start a new z value

		if(iMaterial != previousMaterial)
		{
			newiz++;
			compressionIndex=1;
			*(newzValues+newiz) = *(oldzValues+oldiz);
			*(newMaterialTransitions+iMaterial) = newiz;
			printf("New material starts at index %d which was index %d in the old material it has zvalue %g\n",newiz,oldiz,*(newzValues+newiz));
		}
		else
		{
			if(compressionIndex == zCompression)
			{
				newiz++;
				compressionIndex=1;
				*(newzValues+newiz) = *(oldzValues+oldiz);
			}
			else
			compressionIndex++;
		}
		*(newFormulaUnits+newiz) += *(oldFormulaUnits+oldiz);
		*(newThickness+newiz) += *(oldThickness+oldiz);
		
		if(*(architecture->absorbers+iMaterial))
		{
		//numAbsorberLayers += *(oldFormulaUnits+oldIndex);
		numAbsorberLayers++;
		}
		previousMaterial=iMaterial;
	}
	architecture->zValues=newzValues;
	architecture->formulaUnits = newFormulaUnits;
	architecture->thickness = newThickness;
	architecture->materialTransitions = newMaterialTransitions;
	architecture->numzValues = newiz + 1;
	architecture->numStates = 2*architecture->numzValues + architecture->numDefectStates;
	architecture->numAbsorberLayers=numAbsorberLayers;
	free(oldzValues);
	free(oldFormulaUnits);
	free(materialTransitions);
}


void stacd_generation(Configuration *config, double generationPerCubicAngstrom)
{
	CarrierProfile *data = config->data;
	crystal *crys = config->crys;
	
	int numzValues = architecture->numzValues;
	double *generation=data->generation;
	int *absorbers = architecture->absorbers;
	double *formulaUnits=architecture->formulaUnits;
	double *bandOffsets=architecture->bandOffsets;
	double temperature = getTemperature();
	//printf("Generation per cubic angstrom is %g\n",generationPerCubicAngstrom);
	int iz, iMaterial;
	double generationDensity,materialThermalRate, stateVolume;
	for(iz=0;iz<numzValues;iz++)
	{
	iMaterial = zMaterial(architecture,iz);
	stateVolume = architecture->volumePerUnitCell * *(architecture->formulaUnits+iz);
		if(*(absorbers+iMaterial) && generationPerCubicAngstrom != 0)
			generationDensity =  generationPerCubicAngstrom;
		else 
		if(*(bandOffsets+2*iMaterial) != *(bandOffsets+2*iMaterial+1))
			generationDensity = pow(thermalGenerationConstant,2) * bandRecombinationRate * exp((-*(bandOffsets+2*iMaterial)+*(bandOffsets+2*iMaterial+1))/temperature);
		else
			generationDensity = 0;//zero generation for metals
		*(generation+iz) = generationDensity* architecture->volumePerUnitCell * *(formulaUnits+iz);
		//printf("Generation for layer %d is %g\n",iz,*(generation+iz));
	}
}

void zIndexWeights(double *weight1, double *weight2, int *iNeighbor, double zValue, int iz)
{
	double *zValues = architecture->zValues;
	double *thickness = architecture->thickness;
	double distance1,distance2;
	int left = zValue < *(zValues+iz)+*(thickness+iz)/2;
	distance1 = fabs(zValue - *(zValues+iz)-*(thickness+iz)/2);
	if(left)
	{
		distance2 = fabs(zValue-*(zValues+iz-1)-*(thickness+iz-1)/2);
		*iNeighbor = -1;
	}
	else
	{
		distance2 = fabs(*(zValues+iz+1)+*(thickness+iz+1)/2-zValue);
		*iNeighbor = 1;
	}
	*weight1 = distance2/(distance1+distance2);
	*weight2 = distance1/(distance1+distance2);
}

double stacd_charge(Configuration *config)
{
	CarrierProfile *data = config->data;
	crystal *crys = config->crys;
	
	double *layerCharge=data->layerCharge; 
    double *ionicCharge=data->ionicCharge; 
    double *electronCharge=data->electronCharge; 
    double *holeCharge=data->holeCharge;
    double *stateCharge=data->stateCharge;
    double *ionCharge = architecture->ionCharge;
    double *effectiveCharge=data->effectiveCharge;
    double *defectLevelCharge = architecture->defectLevelCharge;
  
	double *formulaUnits=architecture->formulaUnits;
	int *numIonDefectLevels = architecture->numIonDefectLevels;
	int *ionAtoms = architecture->ionAtoms;
	int *ionEles = architecture->ionEles;
    int numzValues = architecture->numzValues;
	int numIons = architecture->numIons;
	int numStates = architecture->numStates;
	
	double *defectDensities = architecture->defectDensities;
	double *bandDefectCharge = architecture->bandDefectCharge;
	double volumePerUnitCell = architecture->volumePerUnitCell;
	
	double *zValues = architecture->zValues;
	double *thickness = architecture->thickness;
	double weight1, weight2;
	//double distance0,distance1;
	int iNeighbor;
	
	int iState, iBand, iCharge, iz, iMaterial, iIon, iEle, iAtom, iDL, ionDL, iDS;
	for(iz=0;iz<numzValues;iz++)
	{
		iBand = iz;
		stacd_iBandToState(crys, &iState, &iBand, &iCharge, &iz, &iMaterial, &iIon, &iEle, &iAtom, &iDL, &ionDL, &iDS);
		*(ionicCharge+iz) = *(defectDensities+iMaterial) * *(bandDefectCharge+iMaterial) * *(formulaUnits+iz) * volumePerUnitCell;
		*(electronCharge+iz) = -*(stateCharge+iz);
		*(holeCharge+iz) = *(stateCharge+iz+numzValues);
		//printf("Layer %d electron charge %g hole charge %g\n",iz,*(electronCharge+iz),*(holeCharge+iz));
	}
	
	//This loop will need to be rewritten when including arbitrary defect levels.
	for(iIon=0;iIon<numIons;iIon++)
	{
		stacd_iIonToState(crys, &iState, &iBand, &iCharge, &iz, &iMaterial, &iIon, &iEle, &iAtom, &iDL, &ionDL, &iDS);
		
		zIndexWeights(&weight1, &weight2, &iNeighbor,  *(crys->positions+3*iAtom+DIRECTION_Z), iz);
		/*
		left = *(crys->positions+3*iAtom+2) < *(zValues+iz)+*(thickness+iz)/2;
		
		distance0 = fabs(*(crys->positions+3*iAtom+2) - *(zValues+iz)-*(thickness+iz)/2);
		if(left)
		{
			distance1 = fabs(*(crys->positions+3*iAtom+2)-*(zValues+iz-1)-*(thickness+iz-1)/2);
			iNeighbor = -1;
		}
		else
		{
			distance1 = fabs(*(zValues+iz+1)+*(thickness+iz+1)/2-*(crys->positions+3*iAtom+2));
			iNeighbor = 1;
		}
		weight1 = distance1/(distance1+distance0);
		weight2 = distance0/(distance1+distance0);
		*/
		//Here we merge the charge from the defect to its corresponding layer
		*(ionicCharge + iz) += weight1 * *(ionCharge + iEle);//Note this uses the nominal charge of the ion
		*(ionicCharge + iz + iNeighbor) += weight2 * *(ionCharge + iEle);//Note this uses the nominal charge of the ion

		*(effectiveCharge+iIon) = *(ionCharge + iEle);
		//Merge the trap charge into the layers
		for(ionDL=0;ionDL<*(numIonDefectLevels+iEle);ionDL++,iDL++,iState++)
		{
			if(*(defectLevelCharge+iDL) < 0)
			{
			*(electronCharge+iz) += weight1 * *(defectLevelCharge+iDL) * *(stateCharge+iState);
			*(electronCharge+iz+iNeighbor) += weight2 * *(defectLevelCharge+iDL) * *(stateCharge+iState);
			}
			else
			{
			*(holeCharge+iz) += weight1 * *(defectLevelCharge+iDL) * *(stateCharge+iState);	 
			*(holeCharge+iz+iNeighbor) += weight2 * *(defectLevelCharge+iDL) * *(stateCharge+iState);	 
			}
		//Calculate the effective charge which is nominal charge + trap states on the defect
		*(effectiveCharge+iIon) += *(defectLevelCharge+iDL) * *(stateCharge+iState);
		}
		//printf("Effective charge of ion %d is %g\n",iIon,*(effectiveCharge+iIon));
	}
	double chargeSum=0;
	for(iz=0;iz<numzValues;iz++)
	{
		*(layerCharge+iz) = *(ionicCharge+iz) + *(electronCharge+iz) + *(holeCharge+iz);
		//printf("Charge of layer %d is %g due to %g ions %g electrons %g holes\n",iz, *(layerCharge+iz),*(ionicCharge+iz), *(electronCharge+iz), *(holeCharge+iz));
		chargeSum += *(layerCharge+iz);
		//printf("Layer %d electron charge %g hole charge %g ionic charge %g total charge %g\n",iz,*(electronCharge+iz),*(holeCharge+iz),*(ionicCharge+iz),*(layerCharge+iz));
	}
	return chargeSum;
}
double stacd_chargeCorrection(Configuration *config, double chargeSum)
{
	CarrierProfile *data = config->data;
	crystal *crys = config->crys;
	double *stateCharge=data->stateCharge;
	int numzValues = architecture->numzValues;
	//Place correction charge on the metals
	if(chargeSum > 0)
	{
	*(stateCharge) += chargeSum/2;
	*(stateCharge+numzValues-1) += chargeSum/2;
	}
	else
	{
	*(stateCharge+numzValues) += -chargeSum/2;
	*(stateCharge+2*numzValues-1) += -chargeSum/2;
	}
}


double stacd_potential(Configuration *config, double chargeSum, double potentialDifference)
{
	CarrierProfile *CP = config->data;
	crystal *crys = config->crys;
	double *layerCharge=CP->layerCharge; 
    double *stateCharge=CP->stateCharge;
    double *potential = CP->potential;
	double *field = CP->field;
    double *zValues = architecture->zValues;
    double *thickness = architecture->thickness;
    double *dielectricConstants = architecture->dielectricConstants;
    double *quasiFermiLevels=CP->quasiFermiLevels;
    double *bandOffsets = architecture->bandOffsets;
    double *defectLevels = architecture->defectLevels;
    double *formulaUnits=architecture->formulaUnits;
    int *ionAtoms = architecture->ionAtoms;
    int *ionEles = architecture->ionEles;
    int numzValues = architecture->numzValues;
	int numIons = architecture->numIons;
	int numStates = architecture->numStates;
	int numDefectStates = architecture->numDefectStates;
	double *stateCurrent = CP->stateCurrent;
	double *stateCurrent2 = CP->stateCurrent2;
	double *defectDensities = architecture->defectDensities;
	double *bandDefectLevels = architecture->bandDefectLevels;
	double *bandDefectCharge = architecture->bandDefectCharge;
	
	double area = *(crys->latticeVectors) * *(crys->latticeVectors+4);
	int iState, iBand, iCharge, iz, iMaterial, iIon, iEle, iAtom, iDL, ionDL, iDS;
	double temp = getTemperature();
	double chargeBelow, chargeAbove,deltaz;
	double stateVolume;
	double Vion;
	double oldPotential,convergence = 0;

	
	
	*(potential+numzValues-1) = -potentialDifference/2;//Setting the right contact to zero matches with SCAPS
	*(potential) = potentialDifference/2 + *(architecture->bandOffsets) - *(architecture->bandOffsets + 2*architecture->numMaterials-2);//Setting the right contact to zero matches with SCAPS
	*(quasiFermiLevels) = *(architecture->bandOffsets);//Set the left and right contacts to match the quasi fermi levels of the metal
	*(quasiFermiLevels+numzValues) = *(architecture->bandOffsets+1);//Set the left and right contacts to match the quasi fermi levels of the metal
	*(quasiFermiLevels+numzValues-1) = *(architecture->bandOffsets + 2*architecture->numMaterials-2);//Set the left and right contacts to match the quasi fermi levels of the metal
	*(quasiFermiLevels+2*numzValues-1) = *(architecture->bandOffsets + 2*architecture->numMaterials-1);//Set the left and right contacts to match the quasi fermi levels of the metal
	
	
	//*(potential+numzValues-1) = 0;//Setting the right contact to zero matches with SCAPS
	//*(potential) = 0;//Setting the right contact to zero matches with SCAPS
	//*(quasiFermiLevels) = *(architecture->bandOffsets);//Set the left and right contacts to match the quasi fermi levels of the metal
	//*(quasiFermiLevels+numzValues) = *(architecture->bandOffsets+1);//Set the left and right contacts to match the quasi fermi levels of the metal
	//*(quasiFermiLevels+numzValues-1) = *(architecture->bandOffsets + 2*architecture->numMaterials-2)+potentialDifference;//Set the left and right contacts to match the quasi fermi levels of the metal
	//*(quasiFermiLevels+2*numzValues-1) = *(architecture->bandOffsets + 2*architecture->numMaterials-1)+potentialDifference;//Set the left and right contacts to match the quasi fermi levels of the metal
	
	
	//*(potential+numzValues-1) = 0;//Setting the right contact to zero matches with SCAPS
	//*(potential) = 0;//Setting the right contact to zero matches with SCAPS
	//*(quasiFermiLevels) = *(architecture->bandOffsets)-potentialDifference/2;//Set the left and right contacts to match the quasi fermi levels of the metal
	//*(quasiFermiLevels+numzValues) = *(architecture->bandOffsets+1)-potentialDifference/2;//Set the left and right contacts to match the quasi fermi levels of the metal
	//*(quasiFermiLevels+numzValues-1) = *(architecture->bandOffsets + 2*architecture->numMaterials-2)+potentialDifference/2;//Set the left and right contacts to match the quasi fermi levels of the metal
	//*(quasiFermiLevels+2*numzValues-1) = *(architecture->bandOffsets + 2*architecture->numMaterials-1)+potentialDifference/2;//Set the left and right contacts to match the quasi fermi levels of the metal
	
	//*(potential+numzValues-1) = 0;//Setting the right contact to zero matches with SCAPS
	//*(potential) = 0;//Setting the right contact to zero matches with SCAPS
	//*(quasiFermiLevels) = *(architecture->bandOffsets)-potentialDifference/2;//Set the left and right contacts to match the quasi fermi levels of the metal
	//*(quasiFermiLevels+numzValues) = *(architecture->bandOffsets+1)-potentialDifference/2;//Set the left and right contacts to match the quasi fermi levels of the metal
	//*(quasiFermiLevels+numzValues-1) = *(architecture->bandOffsets)+potentialDifference/2;//Set the left and right contacts to match the quasi fermi levels of the metal
	//*(quasiFermiLevels+2*numzValues-1) = *(architecture->bandOffsets+1)+potentialDifference/2;//Set the left and right contacts to match the quasi fermi levels of the metal
	
	
	
	for(iz=numzValues-1,chargeAbove=0;iz>=0;iz--)
	{
		iMaterial = zMaterial(architecture,iz);
		if(iBand == numzValues-1)
			deltaz = *(zValues+numzValues-1)-*(zValues+numzValues-2)+*(thickness+numzValues-1)/2-*(thickness+numzValues-2)/2;
		else
			deltaz = *(zValues+iz+1)-*(zValues+iz)+*(thickness+iz+1)/2-*(thickness+iz)/2;
		stateVolume = architecture->volumePerUnitCell * *(formulaUnits+iz);

		if(iz != numzValues-1 && iz != 0)
		{
		oldPotential = *(potential+iz);	
		//*(potential+iz) = 
		//((*(potential+iz-1)+*(potential+iz+1))/2 + *(layerCharge+iz)/stateVolume * deltaz*deltaz/ permittivity/ *(dielectricConstants+iMaterial));
		*(potential+iz) = (1-CP->dampingParameter) * *(potential+iz) + CP->dampingParameter *
		((*(potential+iz-1)+*(potential+iz+1))/2 + *(layerCharge+iz)/stateVolume * deltaz*deltaz/ permittivity/ *(dielectricConstants+iMaterial));
		*(field + iz) = (*(potential+iz+1)-*(potential+iz-1))/2/deltaz;
		//if(iMaterial == 4)
		//printf("Charge of layer %d is %g and it has potential %g compared to its neighbors %g and %g\n",iz, *(layerCharge+iz), *(potential+iz), *(potential+iz-1), *(potential+iz+1));
		convergence += pow(oldPotential-*(potential+iz),2);
		}
		if(iz == 0)
		*(field + iz) = (*(potential+iz+1)-*(potential+iz))/deltaz;
		if(iz == numzValues-1)
		*(field + iz) = (*(potential+iz)-*(potential+iz-1))/deltaz;
		
		
		for(iCharge = 0;iCharge<2;iCharge++)
		{
			//This should be replaced by the index lookup function
			iBand = iz + iCharge*numzValues;
			iState=iBand;
			if(*(bandOffsets+2*iMaterial) == *(bandOffsets+2*iMaterial+1))// calculating the charge on the metals may needs to be redone. 
				{
				//*(quasiFermiLevels+iState) = *(bandOffsets+2*iMaterial);
				//*(stateCharge+iState) = bandFermiLevelToOccupancy(*(quasiFermiLevels+iState) + *(potential+iz) - *(bandOffsets+2*iMaterial+iCharge),stateVolume,iCharge);		
				*(stateCharge+iState) = 1e-10;		
				}
			else
			{
			
			if((*(bandDefectCharge+iMaterial) > 0 && iCharge == 0) || (*(bandDefectCharge+iMaterial) < 0 && iCharge == 1))
			*(stateCharge+iState) = bandDefectFermiLevelToOccupancy(*(quasiFermiLevels+iState) + *(potential+iz) -  *(bandOffsets+2*iMaterial+iCharge),  *(defectDensities+iMaterial), stateVolume, temp, (bandDefectLevels+iMaterial), 1, iCharge);
			else		
				*(stateCharge+iState) = bandFermiLevelToOccupancy(*(quasiFermiLevels+iState) + *(potential+iz) - *(bandOffsets+2*iMaterial+iCharge),stateVolume,iCharge);			
			//printf("Quasifermilevels of band state %d is %g based on occupancies %g potential %g band offset %g layer charge %g\n",iBand,*(quasiFermiLevels+iState),*(stateCharge+iState),*(potential+iz),*(bandOffsets+2*iMaterial+iCharge), *(layerCharge+iz) );
			//printf("Layer %d icharge %d F %g n %g\n",iz,iCharge,*(quasiFermiLevels+iState) + *(potential+iz) - *(bandOffsets+2*iMaterial+iCharge),*(stateCharge+iState));
			}
			iMaterial = zMaterial(architecture,iz);
			if(iz != 0 && iz != numzValues-1)
			{
			//J=sigma * E = mu n e (- Delta (F+V) / Delta x)
			//*(stateCurrent+iBand) = mobility * *(stateCharge+iState)/ (architecture->volumePerUnitCell* *(formulaUnits+iz)) *
			 //-pow(-1,iCharge) * ((*(potential+iz+1) + *(quasiFermiLevels+iState+1)) - (*(potential+iz-1)+*(quasiFermiLevels+iState-1)))/2/deltaz; 
			//*(stateCurrent+iBand) = mobility * *(stateCharge+iState)/ (architecture->volumePerUnitCell* *(formulaUnits+iz)) *
			 //-pow(-1,iCharge) * (*(quasiFermiLevels+iState+1) - *(quasiFermiLevels+iState-1))/2/deltaz;
			 //*(stateCurrent+iBand) = mobility / (architecture->volumePerUnitCell* *(formulaUnits+iz)) *
			 //-pow(-1,iCharge) * (*(quasiFermiLevels+iState+1) - *(quasiFermiLevels+iState-1))/2/deltaz; 
			deltaz = *(zValues+iz+1)-*(zValues+iz-1)+*(thickness+iz+1)/2-*(thickness+iz-1)/2;
			*(stateCurrent+iBand) = mobility/ (architecture->volumePerUnitCell* *(formulaUnits+iz)) *
			 -pow(-1,iCharge) * (*(quasiFermiLevels+iState+1) - *(quasiFermiLevels+iState-1))/2.0/deltaz; 
			 //*(stateCurrent2+iBand) = mobility / (architecture->volumePerUnitCell* *(formulaUnits+iz)) * *(stateCharge+iState) *
			 //-pow(-1,iCharge) * (( *(quasiFermiLevels+iState+1) - *(potential+iz+1)) - (*(quasiFermiLevels+iState-1) - *(potential+iz-1)))/2/deltaz;
			//*(stateCurrent2+iBand) = mobility / (architecture->volumePerUnitCell* *(formulaUnits+iz)) * *(stateCharge+iState) *
			// -pow(-1,iCharge) * ( *(quasiFermiLevels+iState+1) - *(quasiFermiLevels+iState-1) )/2/deltaz;
			*(stateCurrent2+iBand) = mobility / (architecture->volumePerUnitCell* *(formulaUnits+iz)) * bandFermiLevelToOccupancy(*(quasiFermiLevels+iState) + *(potential+iz) - *(bandOffsets+2*iMaterial+iCharge),stateVolume,iCharge) *
			 -pow(-1,iCharge) * ( *(quasiFermiLevels+iState+1) - *(quasiFermiLevels+iState-1) )/2/deltaz;

			}
			//printf("Quasifermilevels of band state %d is %g based on occupancies %g\n",iBand,*(quasiFermiLevels+iState),*(stateCharge+iState));
		}
		//printf("Carrier Product density of layer %d is %g based on occupancy %g and %g and quasifermilevel splitting is %g\n",iz, *(stateCharge+iz)* *(stateCharge+iz+numzValues)/pow(stateVolume,2),*(stateCharge+iz), *(stateCharge+iz+numzValues),*(quasiFermiLevels+iz)- *(quasiFermiLevels+iz+numzValues));
		
		chargeAbove += *(layerCharge+iz);
	}
	for(iDS=0;iDS<numDefectStates;iDS++)
	{
		stacd_iDSToState(crys, &iState, &iBand, &iCharge, &iz, &iMaterial, &iIon, &iEle, &iAtom, &iDL, &ionDL, &iDS);		
		Vion = zPotential(architecture,CP,*(crys->positions+3*iAtom+2));
		*(quasiFermiLevels+iState) = zFermi(architecture,CP,*(crys->positions+3*iAtom+2),iCharge);
		*(stateCharge + iState) =  defectFermiLevelToOccupancy(*(quasiFermiLevels+iState) + Vion - (*(bandOffsets+2*iMaterial+iCharge) + *(defectLevels+iDL)),iCharge);
		/*printf("Quasifermilevels of ion state %d which is ele %d iCharge %d are %g from band fermi level %g"
		" based on occupancy %g band occupancy %g Band offset %g defect level %g and potential %g band potential %g zvalue %g band iz %g\n"
		,iState,iEle,iCharge,*(quasiFermiLevels+iState),*(quasiFermiLevels+iBand),
		*(stateCharge+iState),*(stateCharge+iBand),*(bandOffsets+2*iMaterial+iCharge),*(defectLevels+iDL),Vion,*(potential + iz),*(crys->positions+3*iAtom+2),*(zValues+iz));
		*/
		/*printf("Quasifermilevels of ion state %d which is ele %d iCharge %d are %g from band fermi level %g"
		" based on occupancy %g band occupancy %g Band offset %g defect level %g and potential %g band potential %g zvalue %g band iz %g\n"
		,iState,iEle,iCharge,*(quasiFermiLevels+iState) + Vion - (*(bandOffsets+2*iMaterial+iCharge) + *(defectLevels+iDL)),*(quasiFermiLevels+iBand),
		*(stateCharge+iState),*(stateCharge+iBand),*(bandOffsets+2*iMaterial+iCharge),*(defectLevels+iDL),Vion,*(potential + iz),*(crys->positions+3*iAtom+2),*(zValues+iz));
		*/
	}
	return sqrt(convergence/numzValues)/CP->dampingParameter;
}

//#define iFromInterface 12 //This arbitrary distance from the interfacemake it so the current is counted from the ETL/HTL not the metal
void stacd_chargeExtraction(Configuration *config)
{
	CarrierProfile *CP = config->data;
    double *stateCharge=CP->stateCharge;
    double *stateCurrent2=CP->stateCurrent2;
    int numzValues = architecture->numzValues;
    
    int iz;
    double eLeft, hLeft, eRight, hRight;
    for(iz = 1;iz<numzValues-1;iz++)
    {
		if(zMaterial(architecture,iz) == 2 && zMaterial(architecture,iz-1) == 1)
		{
			eLeft = *(stateCurrent2 + iz);
			hLeft =  *(stateCurrent2+numzValues+iz);
		}
		if(zMaterial(architecture,iz) == 2 && zMaterial(architecture,iz+1) == 3)
		{
			eRight = *(stateCurrent2 + iz);
			hRight =  *(stateCurrent2+numzValues+iz);
		}
	}
   // double eLeft = *(stateCurrent + iFromInterface);
    //double eRight = *(stateCurrent+numzValues-iFromInterface-1);
    //double hLeft =  *(stateCurrent+numzValues+iFromInterface);
    //double hRight = *(stateCurrent+2*numzValues-iFromInterface-1);
    double eLhR = (-eLeft+hRight)/2;
    double eRhL = (-eRight+hLeft)/2;
    CP->chargeInjection = eLhR-eRhL;
    //printf("Charge injection is %g\n",CP->chargeInjection);
}


double stacd_fermiSolver(Configuration *config)
{
	/*This function implements the equation
	* F(i) = 1/(n(k)+n(k+1) * (n(k) F(i-1) + n(k+1) F(i+1)) 
	* + dx2/(2 mu) * ((G(k)-R(k))/n(k) +(G(k+1)-R(k+1))/n(k+1))
	* Where F is the carrier-specific quasi fermi level, n is a carrier density
	* dx2 is the square of the x displacement, mu is the mobility, G is the carrier generation
	* R(k) = R n(k) p(k) is the carrier recombination with opposite charge carrer density p(k).
	* F(k) = (F(i-1) + F(i))/2
	* F(k+1) = (F(i)+F(i+1))/2
	* The k grid is interstitial to the i grid in the above equation and is found by averaging potentials.  
	*/
	CarrierProfile *CP = config->data;
	crystal *crys = config->crys;
	double area = *(crys->latticeVectors) * *(crys->latticeVectors+4);
    int numzValues = architecture->numzValues;
	int numIons = architecture->numIons;
	int numStates = architecture->numStates;
	//printf("There are %d z Values and %d ions for %d states\n",numzValues,numIons,numStates);
	
	int *ionAtoms = architecture->ionAtoms;
	double *stateCharge=CP->stateCharge;
	double *quasiFermiLevels=CP->quasiFermiLevels;
	double *stateRecombination=CP->stateRecombination;
	double *stateChange=CP->stateChange;
	double *generation=CP->generation;
	double *potential=CP->potential;
	double *zValues = architecture->zValues;
	double *formulaUnits=architecture->formulaUnits;
	double *thickness=architecture->thickness;
	double *bandOffsets=architecture->bandOffsets;
	double *bandDefectLevels=architecture->bandDefectLevels;
	double *defectDensities=architecture->defectDensities;
    int *absorbers=architecture->absorbers;	
    
	int iState1, iBand1, iCharge1, iz1, iMaterial1, iIon1, iEle1, iAtom1, iDL1, ionDL1, iDS1;
	int iState2, iBand2, iCharge2, iz2, iMaterial2, iIon2, iEle2, iAtom2, iDL2, ionDL2, iDS2;	
	int iState3, iBand3, iCharge3, iz3, iMaterial3, iIon3, iEle3, iAtom3, iDL3, ionDL3, iDS3;	
	
	int iMaterialHasDefects;
	int oCharge; //opposite charge
	
	double temp = getTemperature();
	
	double statekVolume, statekp1Volume;
	
	double deltaz;
	
	double Fk,Fkp1,pFk,pFkp1,Vk,Vkp1,nk,nkp1,pk,pkp1,Gk,Gkp1,Rk,Rkp1,BOk,BOkp1,pBOk,pBOkp1;
	
	double oldFermi;
	double convergence = 0;
	double egenerationSum = 0;
	double erecombinationSum = 0;
	double ekSum = 0;
	double ekp1Sum = 0;
	double hgenerationSum = 0;
	double hrecombinationSum = 0;
	double hkSum = 0;
	double hkp1Sum = 0;
	int denomFix = 0;
	
	for(iState1=0;iState1<numStates;iState1++)
	{	
		*(stateRecombination + iState1) = 0;
		*(stateChange + iState1) = 0;		
	}
	
	//Fill transfer and recombination tensors
	for(iCharge2 = 0;iCharge2<2;iCharge2++)
	for(iz2=1;iz2<numzValues-1;iz2++)
	{
		iState1 = (iz2 - 1) + iCharge2*numzValues; //i-1
		iState2 = (iz2) + iCharge2*numzValues; //i
		iState3 = (iz2 + 1) + iCharge2*numzValues; //i+1
		stacd_stateToIndex(crys, &iState1, &iBand1, &iCharge1, &iz1, &iMaterial1, &iIon1, &iEle1, &iAtom1, &iDL1, &ionDL1, &iDS1);
		stacd_stateToIndex(crys, &iState2, &iBand2, &iCharge2, &iz2, &iMaterial2, &iIon2, &iEle2, &iAtom2, &iDL2, &ionDL2, &iDS2);
		stacd_stateToIndex(crys, &iState3, &iBand3, &iCharge3, &iz3, &iMaterial3, &iIon3, &iEle3, &iAtom3, &iDL3, &ionDL3, &iDS3);
		
		oCharge = 1-iCharge2;	
		//printf("iz1 %d iz2 %d iz3 %d iCharge1 %d iCharge2 %d iCharge3 %d\n",iz1,iz2,iz3,iCharge1,iCharge2,iCharge3);
			
		statekVolume = architecture->volumePerUnitCell * (*(formulaUnits+iz1)+*(formulaUnits+iz2))/2;
		statekp1Volume = architecture->volumePerUnitCell * (*(formulaUnits+iz2)+*(formulaUnits+iz3))/2;
		
		Fk = (*(quasiFermiLevels+iState1) + *(quasiFermiLevels+iState2))/2;
		Fkp1 = (*(quasiFermiLevels+iState2) + *(quasiFermiLevels+iState3))/2;
		pFk = (*(quasiFermiLevels + (iz1 - 1) + oCharge * numzValues ) + *(quasiFermiLevels+ (iz1) + oCharge *numzValues ))/2;//Fermi levels of opposite charge state
		pFkp1 = (*(quasiFermiLevels + iz1 + oCharge * numzValues ) + *(quasiFermiLevels + (iz1 + 1) + oCharge *numzValues ))/2;
		Vk = (*(potential+iz1)+*(potential+iz2))/2;
		Vkp1 = (*(potential+iz2)+*(potential+iz3))/2;
		BOk = (*(bandOffsets+2*iMaterial1+iCharge2)+*(bandOffsets+2*iMaterial2+iCharge2))/2;
		BOkp1 = (*(bandOffsets+2*iMaterial2+iCharge2)+*(bandOffsets+2*iMaterial3+iCharge2))/2;
		pBOk = (*(bandOffsets+2*iMaterial1+oCharge)+*(bandOffsets+2*iMaterial2+oCharge))/2;//Band offset of opposite charge
		pBOkp1 = (*(bandOffsets+2*iMaterial2+oCharge)+*(bandOffsets+2*iMaterial3+oCharge))/2;
		
		if(*(defectDensities+iMaterial1)>0)//Determines which material to use when calculating defect densities at interfaces. 
		iMaterialHasDefects = iMaterial1;//This has a downside that if two materials with defects are joined at an interface it will not work properly.
		else if(*(defectDensities+iMaterial2)>0) 
		iMaterialHasDefects = iMaterial2;
		else if(*(defectDensities+iMaterial3)>0) 
		iMaterialHasDefects = iMaterial3;
		else
		iMaterialHasDefects = iMaterial1;

		nk = bandDefectFermiLevelToOccupancy(Fk + Vk -  BOk,  *(defectDensities+iMaterialHasDefects), statekVolume, temp, (bandDefectLevels+iMaterialHasDefects),1, iCharge2);//This needs to be amended to consider defect density at interfaces
		nkp1 = bandDefectFermiLevelToOccupancy(Fkp1 + Vkp1 -  BOkp1,  *(defectDensities+iMaterialHasDefects), statekp1Volume, temp, (bandDefectLevels+iMaterialHasDefects), 1, iCharge2);//This needs to be amended to consider defect density at interfaces
		//nk/pk is bad nomenclature. nk refers to the charge carrier with iCharge1=iCharge2=iCharge3. pk refers to !iCharge1. This is different than the standard nomenclature where nk is electrons.
		pk = bandDefectFermiLevelToOccupancy(pFk + Vk -  pBOk,  *(defectDensities+iMaterialHasDefects), statekVolume, temp, (bandDefectLevels+iMaterialHasDefects), 1, oCharge); 
		pkp1 = bandDefectFermiLevelToOccupancy(pFkp1 + Vkp1 -  pBOkp1,  *(defectDensities+iMaterialHasDefects), statekp1Volume, temp, (bandDefectLevels+iMaterialHasDefects), 1, oCharge);
		
		/*printf("Layer %d charge %d Fk %g BOk %g nk %g Fkp1 %g nkp1 %g pFk %g pk %g pFkp1 %g pkp1 %g\n",iz2,iCharge2,
		Fk+Vk-BOk,BOk,nk,Fkp1+Vkp1-BOkp1,nkp1,
		pFk+Vk-pBOk,pk,pFkp1+Vkp1-pBOkp1,pkp1);
		*/
		Gk = (*(generation+iz1)+*(generation+iz2))/2;
		Gkp1 = (*(generation+iz2)+*(generation+iz3))/2;
		if((*(bandOffsets+2*iMaterial1) != *(bandOffsets+2*iMaterial1+1)) && (*(bandOffsets+2*iMaterial2) != *(bandOffsets+2*iMaterial1+1)))
			Rk = bandRecombinationRate * nk * pk / statekVolume;
		else
			Rk = 0;
		if((*(bandOffsets+2*iMaterial2) != *(bandOffsets+2*iMaterial2+1)) && (*(bandOffsets+2*iMaterial3) != *(bandOffsets+2*iMaterial3+1)))	
			Rkp1 = bandRecombinationRate * nkp1 * pkp1 / statekp1Volume;
		else
			Rkp1 = 0;
		deltaz = (*(thickness+iz1)+*(thickness+iz2))/2;
		
		
		//if(iMaterial1 != iMaterial2 || iMaterial2 != iMaterial3)
		/*printf("Layer %d charge %d materials %d %d %d \nn %g %g %g %g %g p %g %g\n G %g %g R %g %g\n",iz2,iCharge2,iMaterial1,iMaterial2,iMaterial3,
		*(stateCharge+iState1),nk,*(stateCharge+iState2),nkp1,*(stateCharge+iState3),pk,pkp1,Gk,Gkp1,Rk,Rkp1);
		*/
		oldFermi = *(quasiFermiLevels+iState2);
		if(iCharge2 == 0)
		{
		egenerationSum += Gk + Gkp1;
		erecombinationSum += Rk + Rkp1;
		ekSum += (Gk-Rk)/(nk+pk);
		ekp1Sum += (Gkp1-Rkp1)/(nkp1+pkp1);
		}
		else
		{
		hgenerationSum += Gk + Gkp1;
		hrecombinationSum += Rk + Rkp1;
		hkSum += (Gk-Rk)/(nk+pk);
		hkp1Sum += (Gkp1-Rkp1)/(nkp1+pkp1);	
		}

		//It is unclear if dividing by the sum of nk and pk is correct for the generation terms, but it increases the numerical stability a lot
		
		//*(quasiFermiLevels+iState2) = 1/(nk+nkp1)*(nk * *(quasiFermiLevels+iState1) + nkp1 * *(quasiFermiLevels+iState3)) 
		//+ pow(-1,iCharge2) * deltaz*deltaz/(2.0 * mobility) * ( (Gk-Rk)/(nk+pk) + (Gkp1-Rkp1)/(nkp1+pkp1)); 
		*(quasiFermiLevels+iState2) = 1/(nk+nkp1)*(nk * *(quasiFermiLevels+iState1) + nkp1 * *(quasiFermiLevels+iState3)) 
		+ pow(-1,iCharge2) * deltaz*deltaz/(2.0 * mobility) * ( (Gk-Rk)/nk + (Gkp1-Rkp1)/nkp1);
		
		if(!isnormal(*(quasiFermiLevels+iState2)) || fabs(*(quasiFermiLevels+iState2)) > 10)
		{
		//printf("QFL not normal, was %g replacing by approximation nk %g, pk %g npk1 %g pkp1 %g\n",*(quasiFermiLevels+iState2),nk,pk,nkp1,pkp1);
		//printf("Layer %d charge %d materials %d %d %d \nn %g %g %g %g %g p %g %g\n G %g %g R %g %g\n",iz2,iCharge2,iMaterial1,iMaterial2,iMaterial3,
		//*(stateCharge+iState1),nk,*(stateCharge+iState2),nkp1,*(stateCharge+iState3),pk,pkp1,Gk,Gkp1,Rk,Rkp1);
		denomFix++;
		*(quasiFermiLevels+iState2) = 1/(nk+nkp1)*(nk * *(quasiFermiLevels+iState1) + nkp1 * *(quasiFermiLevels+iState3)) 
		+ pow(-1,iCharge2) * deltaz*deltaz/(2.0 * mobility) * ( (Gk-Rk)/(nk+pk) + (Gkp1-Rkp1)/(nkp1+pkp1));
		}
		//else
		
		//printf("New QFL %g\n:",*(quasiFermiLevels+iState2));
		
		*(stateChange+iState2) = *(quasiFermiLevels+iState2)-oldFermi;
		
		*(quasiFermiLevels+iState2) = (1-CP->dampingParameter) * oldFermi + CP->dampingParameter * *(quasiFermiLevels+iState2);
		
		*(stateRecombination+iState2) = (Rk+Rkp1)/2;
		
		convergence += pow(oldFermi-*(quasiFermiLevels+iState2),2);
	}
	//printf("G %g R %g k %g kp1 %g G %g R %g k %g kp1 %g\n",egenerationSum,erecombinationSum,ekSum,ekp1Sum,hgenerationSum,hrecombinationSum,hkSum,hkp1Sum);
	return (1+denomFix)*sqrt(convergence/numzValues)/CP->dampingParameter;
}


void stacd_output(Configuration *config)
{
	CarrierProfile *data = config->data;
	crystal *crys = config->crys;
    int numzValues = architecture->numzValues;
    int numIonTypes = architecture->numIonTypes;
	int numIons = architecture->numIons;
	int numStates = architecture->numStates;
	int numDefectStates = architecture->numDefectStates;
	
	int *ionAtoms = architecture->ionAtoms;
	int *ionEles = architecture->ionEles;
	double *ionCharge = architecture->ionCharge;
	double *stateCharge=data->stateCharge;
	double *quasiFermiLevels=data->quasiFermiLevels;
	double *stateCurrent=data->stateCurrent;
	double *stateCurrent2=data->stateCurrent2;
	double *stateRecombination=data->stateRecombination;
	double *stateChange=data->stateChange;
	double *layerCurrent=data->layerCurrent;
	double *layerCurrent2=data->layerCurrent2;
	double *layerRecombination=data->layerRecombination;
	double *generation=data->generation;
	double *layerChange=data->layerChange;
	double *layerCharge=data->layerCharge;
	double *layerBandCharge=data->layerBandCharge;
	double *ionsPerLayer=data->ionsPerLayer;
	double *potential = data->potential;
	double *effectiveCharge=data->effectiveCharge;
	double *dielectricConstants = architecture->dielectricConstants;
	double *zValues = architecture->zValues;
	double *thickness = architecture->thickness;
    double *formulaUnits=architecture->formulaUnits;	
	int iState, iBand, iCharge, iz, iMaterial, iIon, iEle, iAtom, iDL, ionDL, iDS;
	int iState2, iBand2, iCharge2, iz2, iMaterial2, iIon2, iEle2, iAtom2, iDL2, ionDL2, iDS2;
	
	double weight1, weight2;
	int iNeighbor;
	double maxQFLS = 0;
	
	//Extra calculations of things to be saved
	data->polarization=0; 
	data->ionFieldEnergy=0; 
	data->ionIonEnergy=0;
	data->electronEnergy=0;
	data->holeEnergy=0;
	data->totalEnergy=0;

	for(iz=0;iz<numzValues;iz++)
	{
	if(*(quasiFermiLevels+iz)-*(quasiFermiLevels+iz+numzValues) > maxQFLS)
	maxQFLS = *(quasiFermiLevels+iz)-*(quasiFermiLevels+iz+numzValues);
	for(iEle = 0;iEle < numIonTypes;iEle++)
		*(ionsPerLayer + iEle*numzValues+ iz) = 0;
	}
	for(iBand=0;iBand<2*numzValues;iBand++)
	{
		stacd_iBandToState(crys, &iState, &iBand, &iCharge, &iz, &iMaterial, &iIon, &iEle, &iAtom, &iDL, &ionDL, &iDS);

		//*(layerDiffusion + iBand) =  *(stateDiffusion + iState);
		*(layerCurrent + iBand) =  *(stateCurrent + iState);
		*(layerCurrent2 + iBand) =  *(stateCurrent2 + iState);
		*(layerBandCharge + iBand) =  *(stateCharge + iState);
		*(layerRecombination + iBand) = *(stateRecombination + iState);		
		*(layerChange + iBand) = *(stateChange + iState);
		
		data->polarization += (*(zValues+iz)+*(thickness+iz)/2) * *(layerCharge+iz);
		
		if(iCharge == 0)
		data->electronEnergy += *(stateCharge+iState)  * (*(quasiFermiLevels+iState)+ *(potential+iz));
		else
		data->holeEnergy += -*(stateCharge+iState) * (*(quasiFermiLevels+iState)+ *(potential+iz));		
	}
	double Vion;
	for(iDS=0;iDS<numDefectStates;iDS++)
	{
		stacd_iDSToState(crys, &iState, &iBand, &iCharge, &iz, &iMaterial, &iIon, &iEle, &iAtom, &iDL, &ionDL, &iDS);

		zIndexWeights(&weight1, &weight2, &iNeighbor, *(crys->positions+3*iAtom+DIRECTION_Z), iz);
		
		//*(layerDiffusion + iBand) += *(stateDiffusion + iState);
		*(layerCurrent + iBand) += weight1 * *(stateCurrent + iState);
		*(layerCurrent2 + iBand) += weight1 * *(stateCurrent2 + iState);
		*(layerRecombination + iBand) += weight1 * *(stateRecombination + iState);
		*(layerChange + iBand) += weight1 * *(stateChange + iState);
		(*(ionsPerLayer + iEle*numzValues+ iz)) += weight1;
		*(layerCurrent + iBand+ iNeighbor) += weight2 * *(stateCurrent + iState);
		*(layerCurrent2 + iBand+ iNeighbor) += weight2 * *(stateCurrent2 + iState);
		*(layerRecombination + iBand+ iNeighbor) += weight2 * *(stateRecombination + iState);
		*(layerChange + iBand+ iNeighbor) += weight2 * *(stateChange + iState);
		(*(ionsPerLayer + iEle*numzValues+ iz + iNeighbor)) += weight2;
		//printf("Ion type %d zvalue %g is in layer %d\n",iEle,*(crys->positions+3*iAtom+DIRECTION_Z),iz);

		
		Vion=zPotential(architecture,data,*(crys->positions+3*iAtom+DIRECTION_Z));
		
		if(iCharge == 0)
		data->electronEnergy += *(stateCharge + iState) * (*(quasiFermiLevels+iState)+ Vion);
		else
		data->holeEnergy += -*(stateCharge + iState) * (*(quasiFermiLevels+iState)+ Vion);	
		
		//The remainder of this loop should only be done once for each ion
		if(ionDL != 0)
		continue;
		


		//This term uses the nominal charge
		data->ionFieldEnergy += *(ionCharge+iEle) * Vion;
		//data->ionFieldEnergy += *(effectiveCharge + iIon) * Vion;
		for(iIon2=iIon+1;iIon2<numIons;iIon2++)
		{
			stacd_iIonToState(crys, &iState2, &iBand2, &iCharge2, &iz2, &iMaterial2, &iIon2, &iEle2, &iAtom2, &iDL2, &ionDL2, &iDS2);
			//This assumes that iIon1 and iIon2 are in the same material. Maybe we should have 2 dielectrics with an effective interaction
			//printf("Q1: %g Q2: %g epsilon: %g r: %g E: %g\n", *(effectiveCharge+iIon1),*(effectiveCharge+iIon2),*(dielectricConstants+iMaterial1),*(interstateDistance + iState1*numStates+iState2),
			data->ionIonEnergy += *(effectiveCharge+iIon) * *(effectiveCharge+iIon2) * coulombconstant/ *(dielectricConstants+iMaterial) / crys_atomDistance(crys, *(ionAtoms+iIon), *(ionAtoms+iIon2));			
		}
		
	}
	
	
	//data->totalEnergy = data->ionFieldEnergy+data->ionIonEnergy+data->electronEnergy+data->holeEnergy;
	data->totalEnergy = data->ionFieldEnergy+data->ionIonEnergy;
	//data->totalEnergy = data->ionFieldEnergy;

	//printf("IonField Energy %g IonIonEnergy %g electronEnergy %g holeEnergy %g totalEnergy %g\n",data->ionFieldEnergy,data->ionIonEnergy,data->electronEnergy,data->holeEnergy,data->totalEnergy);
	//This assumes holes come out on the left, electrons on the right
	//data->quasiFermiLevelSplitting = *(quasiFermiLevels+numzValues-1)-*(quasiFermiLevels+numzValues);
	data->quasiFermiLevelSplitting = maxQFLS;
	//printf("Quasi Fermi Levels %g %g %g %g resulting in splitting %g\n",*(quasiFermiLevels),*(quasiFermiLevels+numStates),*(quasiFermiLevels+numzValues-1),*(quasiFermiLevels+numzValues-1+numStates),data->quasiFermiLevelSplitting);
	//printf("Total energy %g Quasi Fermi Level splitting %g based on generation %g\n",data->totalEnergy,data->quasiFermiLevelSplitting,*(config->externalConditions+1));
	data->potentialDifference = *(config->externalConditions+0);
	data->carrierGeneration = *(config->externalConditions+1);
	
	//printf("Carrier generation %g Potential difference %g\n",data->carrierGeneration,data->potentialDifference);
	
	int iScalar;
	//This is kind of a ridiculous ask
	if(config->seriesData!=NULL)
	for(iScalar = 0;iScalar<data->numScalars;iScalar++)
	*(config->seriesData+iScalar) = *(*(data->scalars+iScalar));
	config->energy = data->totalEnergy;
	//printf("Charge injection is %g\n",data->chargeInjection);
	//printf("Total energy %g\n",config->energy);
}

double stacd_Energy(Configuration *config)
{
	CarrierProfile *data = config->data;
	crystal *crys = config->crys;
	double generationPerCubicAngstrom = *(config->externalConditions+1) /architecture->numAbsorberLayers/latticeConstant;
	double potentialDifference = *(config->externalConditions+0);

	double fermiConvergence,potentialConvergence;
	double convergence = 1e100;
	int numSteps=0;	
	double chargeSum;
	
	//Filling the generation profile only needs to be done once
	//printf("Carrier Generation %g\n",generationPerCubicAngstrom);
	
	stacd_generation(config,generationPerCubicAngstrom);
	chargeSum = stacd_charge(config);
	stacd_potential(config,chargeSum,potentialDifference);
	//printf("Initial charge of device is %g splitting it amongst band states\n",chargeSum);

	while( convergenceCriterion < convergence)
	{
	//Calculate charges
	chargeSum = stacd_charge(config);
	//printf("Charge of device is %g\n",chargeSum);
	//stacd_chargeCorrection(config,chargeSum);
	chargeSum = stacd_charge(config);
	potentialConvergence = stacd_potential(config,chargeSum,potentialDifference);
	fermiConvergence = stacd_fermiSolver(config);
	//This should be moved to stacd_output
	stacd_chargeExtraction(config);
		
	convergence = fmax(potentialConvergence,fermiConvergence);
	
	if(convergence * data->dampingParameter < convergenceCriterion && data->dampingParameter<maxDampingParameter)
	data->dampingParameter *= 1.001;

	numSteps++;
	//if(numSteps == 10) exit(0);
	
	if(numSteps == maxElectronicIterations)
	break;

	if(numSteps % (maxElectronicIterations/10) == 0)
	printf("Convergence is %g of %g on step %d fC %g pC %g potential %g charge Injection %g generation %g damping %g\n",convergence,convergenceCriterion,numSteps,fermiConvergence,potentialConvergence,*(data->potential),data->chargeInjection,generationPerCubicAngstrom,data->dampingParameter);

	}
	stacd_output(config);
	//if(numSteps == maxElectronicIterations)
	//printf("Convergence is %g of %g on step %d fC %g pC %g potential %g charge Injection %g\n",convergence,convergenceCriterion,numSteps,fermiConvergence,potentialConvergence,*(data->potential),data->chargeInjection);
	//printf("Charge of device is %g\n",chargeSum);
	
	return config->energy;
}
