
#include "electrostatics.h"

char *ionNames;
double *ionCharge;	
int numIons;
//double fieldstrength; // V/angstrom
double voltageOffset;
double voltageModulation;
double *zValues;
int numzValues;
double freeElectronDensity;//e/angstrom^3
double freeHoleDensity;
double screeningLength[2];
//double screeningLength=1e300;
double dielectricConstant;
int vMode;
int stepsPerSwitch;//switches from forward bias to reverse bias with this frequency.
int balanceSteps;

void es_registerSettings()
{
	registerDouble(&freeElectronDensity,"freeElectronDensity",1e-300);
	registerDouble(&freeHoleDensity,"freeHoleDensity",1e-300);
	registerDouble(&dielectricConstant,"dielectricConstant",1);
	
	sourceControl_RegisterSettings();
}

/*
void es_loadSettings(FILE *settingsfile)
{
	readDouble(settingsfile,"voltageOffset",&(voltageOffset));
	readDouble(settingsfile,"voltageModulation",&(voltageModulation));
	readDouble(settingsfile,"freeElectronDensity",&(freeElectronDensity));
	readDouble(settingsfile,"freeHoleDensity",&(freeHoleDensity));
	readDouble(settingsfile,"dielectricConstant",&(dielectricConstant));
	readInt(settingsfile,"stepsPerSwitch",&(stepsPerSwitch));
	readInt(settingsfile,"balanceSteps",&(balanceSteps));
	
	if(freeElectronDensity==0)
	screeningLength[0]=1e10;
	else
	screeningLength[0] = sqrt(dielectricConstant*permittivity*traj_getSettings()->temperature/freeElectronDensity);	
	if(freeHoleDensity==0)
	screeningLength[1]=1e10;
	else
	screeningLength[1] = sqrt(dielectricConstant*permittivity*traj_getSettings()->temperature/freeHoleDensity);	
	char *vModeText = readString(settingsfile,"vMode");
	if(!strcmp(vModeText,"square"))
	{
	vMode=VOLTAGE_SQUARE_WAVE;
	printf("Square Wave\n");
	}
	if(!strcmp(vModeText,"sine"))
	{
	vMode=VOLTAGE_SINE_WAVE;
	printf("Sine Wave\n");
	}

}*/
/*
void es_saveSettings(FILE *settingsfile)
{
	fprintf(settingsfile,"voltageOffset = %g\n",voltageOffset);
	fprintf(settingsfile,"voltageModulation = %g\n",voltageModulation);
	fprintf(settingsfile,"freeElectronDensity = %g\n",freeElectronDensity);
	fprintf(settingsfile,"freeHoleDensity = %g\n",freeHoleDensity);
	fprintf(settingsfile,"dielectricConstant = %g\n",dielectricConstant);
	fprintf(settingsfile,"stepsPerSwitch = %d\n",stepsPerSwitch);
	fprintf(settingsfile,"balanceSteps = %d\n",balanceSteps);
	switch(vMode)
	{
		case VOLTAGE_SQUARE_WAVE:
		fprintf(settingsfile,"vMode = square\n");
		break;
		case VOLTAGE_SINE_WAVE:
		fprintf(settingsfile,"vMode = sine\n");
		break;
	}

}*/

void VP_allocate(VoltageProfile *VP)
{
    VP->voltage = calloc(numzValues,sizeof(double));  
    VP->field = calloc(numzValues,sizeof(double));  
    VP->charge = calloc(numzValues,sizeof(double)); 
    VP->ionicCharge = calloc(numzValues,sizeof(double)); 
    VP->electronCharge = calloc(numzValues,sizeof(double)); 
    VP->holeCharge = calloc(numzValues,sizeof(double)); 
    //The even values are used for positive charge (holes) and the odd values for negative charge (electrons)
    //I believe the factor of 2 is a typo
    VP->ionsPerLayer = calloc(2*numzValues*numIons,sizeof(double));
}

void VP_free(VoltageProfile *VP)
{
    free(VP->voltage);
    free(VP->field);
    free(VP->charge);
    free(VP->ionicCharge);
    free(VP->electronCharge);
    free(VP->holeCharge);
    free(VP->ionsPerLayer);
}

double es_Energy(Configuration *config)
{
	//Calculates the electrostatic energy of charged particles in an electric field.

	VoltageProfile *data = config->data;
	crystal *crys = config->crys;
	double energy, EEenergy=0,Eqenergy=0, qqenergy=0;

	double fieldstrength;
	int iatom1,iatom2;
	int iele1,iele2;
	int ele1start,ele1end;
	int ele2start,ele2end;
	double d;
	int iq1,iq2;
	
	double width = sqrt(*(crys->latticeVectors) * *(crys->latticeVectors+4));
	double height = *(crys->latticeVectors+8)/2;
	double deltaz;
	double convergenceCriterion = 1e-3,previousVoltage;
	double convergence = 1;
	
	double appliedVoltage=*(config->externalConditions);
	double helmholtzParameter;
	double temperature = getTemperature();
	fieldstrength = -appliedVoltage/ height/dielectricConstant;
	double debyePartitionFunction;
	
	for(iele1=0;iele1<numIons;iele1++)
	{
		ele1start = crys_elementOffset(crys,ionNames+iele1*namelength);
		ele1end = ele1start + crys_elementCount(crys,ionNames+iele1*namelength);
		//printf("Element %s from %d to %d has charge %g\n",ionNames+iele1*namelength, ele1start, ele1end, *(ionCharge+iele1));
		for(iatom1=ele1start;iatom1<ele1end;iatom1++)
		{
			Eqenergy -= fieldstrength * *(crys->positions+iatom1*3+2) * *(ionCharge+iele1);
		}
	}
	
	for(iele1=0;iele1<numIons;iele1++)
	{
		ele1start = crys_elementOffset(crys,ionNames+iele1*namelength);
		ele1end = ele1start + crys_elementCount(crys,ionNames+iele1*namelength);
		if(*(ionCharge+iele1)>0)
		iq1 = 0;
		else
		iq1 = 1;
		for(iele2=0;iele2<numIons;iele2++)
		{
			ele2start = crys_elementOffset(crys,ionNames+iele2*namelength);
			ele2end = ele2start + crys_elementCount(crys,ionNames+iele2*namelength);
			if(*(ionCharge+iele2)>0)
			iq2 = 0;
			else
			iq2 = 1;
			for(iatom1=ele1start;iatom1<ele1end;iatom1++)
				for(iatom2=ele2start;iatom2<ele2end;iatom2++)
				{
				if(iatom1 == iatom2)
				continue;	
				d = crys_atomDistance(crys,iatom1,iatom2);
				//Include Thomas-Fermi Screening into charge-charge interaction
				qqenergy += coulombconstant* *(ionCharge+iele1) * *(ionCharge + iele2)*exp(-d/2.0/screeningLength[iq1])*exp(-d/2.0/screeningLength[iq2])/(d*dielectricConstant);	
				}
		}
	}
			
	if(config->savedata)
	{
		double area = *(crys->latticeVectors) * *(crys->latticeVectors+4);
		//printf("Area is %g sq angstrom\n",area);
		double fieldterm,voltageterm;
		//printf("Eq energy %g, qqenergy %g\n",Eqenergy, qqenergy);	
		int iz,iz2,iq,iStep=0;;
		for(iz=0;iz<numzValues;iz++)
		{
			*(data->charge+iz) = 0;
			*(data->electronCharge+iz) = 0;
			*(data->holeCharge+iz) = 0;
			*(data->ionicCharge+iz) = 0;
		}
		
		*(data->voltage)=-appliedVoltage/2.0;
		*(data->voltage+numzValues-1)=appliedVoltage/2.0;
		
		for(iele1=0;iele1<numIons;iele1++)
		for(iz=0;iz<numzValues;iz++)
		*(data->ionsPerLayer+iele1*numzValues+iz)=0;

		for(iele1=0;iele1<numIons;iele1++)
		{
			ele1start = crys_elementOffset(crys,ionNames+iele1*namelength);
			ele1end = ele1start + crys_elementCount(crys,ionNames+iele1*namelength);
				for(iatom1=ele1start;iatom1<ele1end;iatom1++)
				{
				iz = binarysearchDouble(zValues,*(crys->positions+3*iatom1+2),0,numzValues);
				if(iz<0 || iz >= numzValues)
				printf("Bad z index %d from value %g\n",iz,*(crys->positions+3*iatom1+2));
				*(data->ionicCharge+iz) += *(ionCharge+iele1); //The charge is at a lower z value than this index
				(*(data->ionsPerLayer+iele1*numzValues+iz))++;
				}
		}

		for(iz=0;iz<numzValues;iz++)
		*(data->charge+iz) = *(data->ionicCharge+iz);
		
		double debyeQ,extremeVoltage[2],freeCharge[2];
		freeCharge[1]=freeElectronDensity*area*height;
		freeCharge[0]=freeHoleDensity*area*height;
		
		while(convergence > convergenceCriterion)
		{
			
		convergence = 0;
		for(iz=1;iz<numzValues-1;iz++)
		{
			previousVoltage = *(data->voltage+iz);
			*(data->voltage+iz) = 0.5 * (*(data->voltage+iz-1)+*(data->voltage+iz+1)) + pow(0.5*(*(zValues+iz+1)-*(zValues+iz-1)),2) * *(data->charge+iz)/area /(permittivity*dielectricConstant);
			 convergence += fabs((previousVoltage-*(data->voltage+iz))/(fabs(previousVoltage+*(data->voltage+iz))+fabs(appliedVoltage)))/(numzValues-2);
		}	
			
		for(iz=0;iz<numzValues;iz++)
		{
		*(data->charge+iz)=*(data->ionicCharge+iz);
		*(data->electronCharge+iz)=0;
		*(data->holeCharge+iz)=0;
		}
	
			

		iStep++;
		//printf("Convergence%g\n",convergence);
		if((iStep != 0) && (iStep %250 == 0))
		printf("On step %d to convergence, convergence %g\n",iStep,convergence);
		}/*
		while(convergence > convergenceCriterion)
		{
			
		convergence = 0;
		for(iz=1;iz<numzValues-1;iz++)
		{
			previousVoltage = *(data->voltage+iz);
			*(data->voltage+iz) = 0.5 * (*(data->voltage+iz-1)+*(data->voltage+iz+1)) + pow(0.5*(*(zValues+iz+1)-*(zValues+iz-1)),2) * *(data->charge+iz)/area /(permittivity*dielectricConstant);
			 convergence += fabs((previousVoltage-*(data->voltage+iz))/(fabs(previousVoltage+*(data->voltage+iz))+fabs(appliedVoltage)))/(numzValues-2);
		}	
			
		//This is a correction to the charge due to the free ions according to Debye-Huckel theory.
		
		//ExtremeVoltage matches iq to avoid an overflow of the exponential. For positive charge, it is the min Voltage, for negative charge it is the maxVoltage.
		extremeVoltage[0] = extremeVoltage[1] = *(data->voltage);
		for(iz=0;iz<numzValues;iz++)
		{
		*(data->charge+iz)=*(data->ionicCharge+iz);
		*(data->electronCharge+iz)=0;
		*(data->holeCharge+iz)=0;
		if(*(data->voltage+iz)<extremeVoltage[0])
			extremeVoltage[0]=*(data->voltage+iz);
		if(*(data->voltage+iz)>extremeVoltage[1])
			extremeVoltage[1]=*(data->voltage+iz);
		}
		//not allowing the electronic charge to be placed within 2 layers of the endpoints is starting to 
		//reach a simulation with well-defined materials. However, right now it is kind of kludgy
		for(iq=0;iq<2;iq++)
		{
			debyePartitionFunction=0;
			for(iz=2;iz<numzValues-2;iz++)
			//for(iz=0;iz<numzValues;iz++)
			debyePartitionFunction += exp(-pow(-1,iq) * (*(data->voltage+iz)-extremeVoltage[iq])/temperature);
			for(iz=2;iz<numzValues-2;iz++)
			//for(iz=0;iz<numzValues;iz++)
			{
			debyeQ = pow(-1,iq) * freeCharge[iq] * exp(-pow(-1,iq) * (*(data->voltage+iz)-extremeVoltage[iq])/temperature)/debyePartitionFunction;
			*(data->charge+iz) +=  debyeQ;
			if(iq == 0)
			*(data->holeCharge+iz) += debyeQ;
			else
			*(data->electronCharge+iz) += debyeQ;
			//printf("Exponential: %g Debye charge %g\n",exp(-pow(-1,iq) *  (*(data->voltage+iz)-extremeVoltage[iq])/temperature),debyeQ);
			}
		}	
			

		iStep++;
		//printf("Convergence%g\n",convergence);
		if((iStep != 0) && (iStep %250 == 0))
		printf("On step %d to convergence, convergence %g\n",iStep,convergence);
		}*/
		
		data->polarization=0;
		data->ionicPolarization=0;
		data->debyePolarization=0;
		for(iz=0;iz<numzValues;iz++)
		{
			if(iz==0)
			*(data->field) = -(*(data->voltage+1)-*(data->voltage))/(*(zValues+1)-*(zValues));	
			else if(iz ==numzValues-1)
			*(data->field+numzValues-1) = -(*(data->voltage+numzValues-1)-*(data->voltage+numzValues-2))/(*(zValues+numzValues-1)-*(zValues+numzValues-2));	
			else
			*(data->field+iz) = -(*(data->voltage+iz+1)-*(data->voltage+iz-1))/(*(zValues+iz+1)-*(zValues+iz-1));	
			
			data->polarization += *(data->charge+iz) * *(zValues+iz);
			data->ionicPolarization += *(data->ionicCharge+iz) * *(zValues+iz);
			data->debyePolarization += (*(data->holeCharge+iz)+*(data->electronCharge+iz)) * *(zValues+iz);
		}
		
	}
	
	if(!(config->seriesData==NULL))
	{
	*(config->seriesData)=data->polarization;
	*(config->seriesData+1)=data->ionicPolarization;
	*(config->seriesData+2)=data->debyePolarization;
	}
	energy = Eqenergy+qqenergy;

	//printf("Energy is %g saving it to %d\n",energy,(int *)(config->energy));
	(config->energy) = energy;
    
	return energy;
}

double chargeBalance(crystal *crys, double ionicCharge, double interfaceCharge)
{
	double width = sqrt(*(crys->latticeVectors) * *(crys->latticeVectors+4));
	double height = *(crys->latticeVectors+8)/2;
	double freeCharge[2];
	freeCharge[1]=-freeElectronDensity*width*width*height;
	freeCharge[0]=freeHoleDensity*width*width*height;
	double totalCharge = freeCharge[0] + freeCharge[1] + ionicCharge + interfaceCharge;
	return totalCharge;
}

void VP_save(Configuration *config)
{
    
	FILE *outfile = fopen(config->dataFileName,"w");
	VoltageProfile *data = config->data;

    //printf("VP_save %s Config is at %d data is at %d\n",config->dataFileName,(int)config,(int)data);
	int iz,iele;
	fprintf(outfile,"%g\n",data->polarization);
	for(iz=0;iz<numzValues;iz++)
	{
	fprintf(outfile,"%g %g %g %g %g",*(zValues+iz),*(data->charge+iz),*(data->ionicCharge+iz),*(data->electronCharge+iz),*(data->holeCharge+iz),*(data->field+iz),*(data->voltage+iz));
		for(iele=0;iele<numIons;iele++)
			fprintf(outfile," %g",*(data->ionsPerLayer+iele*numzValues+iz));
	
	fprintf(outfile,"\n");
	}
	fclose(outfile);
}

void VP_saveEveryOther(Configuration *config)
{
    
	FILE *outfile = fopen(config->dataFileName,"w");
	VoltageProfile *data = config->data;

    //printf("VP_save %s Config is at %d data is at %d\n",config->dataFileName,(int)config,(int)data);
	int iz,iele;
	fprintf(outfile,"%g\n",data->polarization);
	for(iz=0;iz<numzValues/2;iz++)
	{
	fprintf(outfile,"%g %g %g %g %g %g %g",(*(zValues+2*iz)+*(zValues+2*iz+1))/2,(*(data->charge+2*iz)+*(data->charge+2*iz+1))/2,
	(*(data->ionicCharge+2*iz)+*(data->ionicCharge+2*iz+1))/2,(*(data->electronCharge+2*iz)+*(data->electronCharge+2*iz+1))/2,(*(data->holeCharge+2*iz)+*(data->holeCharge+2*iz+1))/2,
	(*(data->field+2*iz)+*(data->field+2*iz+1))/2,(*(data->voltage+2*iz)+*(data->voltage+2*iz+1))/2);
		for(iele=0;iele<numIons;iele++)
			fprintf(outfile," %g",(*(data->ionsPerLayer+iele*numzValues+2*iz)+*(data->ionsPerLayer+iele*numzValues+2*iz+1))/2.0);
	
	fprintf(outfile,"\n");
	}
	fclose(outfile);
}

void VP_combineWeighted(Configuration *configs, int numCombine,Configuration *outconfig, double *weights)
{
	//printf("Combining %d voltage profiles into one\n",numCombine);
	VoltageProfile *data = outconfig->data;
	VoltageProfile *vp;
	int iz,iConfig,iele,iq;
	double sumWeight=0;
	for(iConfig=0;iConfig<numCombine;iConfig++)
		sumWeight+= *(weights+iConfig);
	for(iz=0;iz<numzValues;iz++)
	{
		 *(data->voltage+iz)=0;
		 *(data->field+iz)=0;
		 *(data->charge+iz)=0;
		 *(data->ionicCharge+iz)=0;
		 *(data->electronCharge+iz)=0;
		 *(data->holeCharge+iz)=0;
		for(iele=0;iele<numIons;iele++)
			*(data->ionsPerLayer+iele*numzValues+iz)=0;

		for(iConfig=0;iConfig<numCombine;iConfig++)
		{
			vp = (configs+iConfig)->data;
			*(data->voltage+iz)+=*(weights+iConfig)**(vp->voltage+iz);
			*(data->field+iz)+=*(weights+iConfig)**(vp->field+iz);
			*(data->charge+iz)+=*(weights+iConfig)**(vp->charge+iz);
			*(data->ionicCharge+iz)+=*(weights+iConfig)**(vp->ionicCharge+iz);
			*(data->electronCharge+iz)+=*(weights+iConfig)**(vp->electronCharge+iz);
			*(data->holeCharge+iz)+=*(weights+iConfig)**(vp->holeCharge+iz);
			for(iele=0;iele<numIons;iele++)
				*(data->ionsPerLayer+iele*numzValues+iz)+=*(weights+iConfig)**(vp->ionsPerLayer+iele*numzValues+iz);

		}
		 *(data->voltage+iz)/=sumWeight;
		 *(data->field+iz)/=sumWeight;
		 *(data->charge+iz)/=sumWeight;
		 *(data->ionicCharge+iz)/=sumWeight;
		 *(data->electronCharge+iz)/=sumWeight;
		 *(data->holeCharge+iz)/=sumWeight;
		// printf("Average charge at height %d is %g\n",iz,*(data->charge+iz)); 
		for(iele=0;iele<numIons;iele++)
			*(data->ionsPerLayer+iele*numzValues+iz)/=sumWeight;
	}
	data->polarization=0;
	data->ionicPolarization=0;
	data->debyePolarization=0;
	for(iConfig=0;iConfig<numCombine;iConfig++)
	{
		vp = (configs+iConfig)->data;
		data->polarization+=*(weights+iConfig)*vp->polarization;
		data->ionicPolarization+=*(weights+iConfig)*vp->ionicPolarization;
		data->debyePolarization+=*(weights+iConfig)*vp->debyePolarization;
	}
	data->polarization/=sumWeight;
	data->ionicPolarization/=sumWeight;
	data->debyePolarization/=sumWeight;
	*(outconfig->seriesData) = data->polarization;	
	*(outconfig->seriesData+1) = data->ionicPolarization;	
	*(outconfig->seriesData+2) = data->debyePolarization;	
	
}

void VP_combine(Configuration *configs, int numcombine,Configuration *outconfig)
{
	//printf("Combining %d voltage profiles into one\n",numcombine);
	VoltageProfile *data = outconfig->data;
	VoltageProfile *vp;
	int iz,iconfig,iele,iq;
	for(iz=0;iz<numzValues;iz++)
	{
		 *(data->voltage+iz)=0;
		 *(data->field+iz)=0;
		 *(data->charge+iz)=0;
		 *(data->ionicCharge+iz)=0;
		 *(data->electronCharge+iz)=0;
		 *(data->holeCharge+iz)=0;
		for(iele=0;iele<numIons;iele++)
			*(data->ionsPerLayer+iele*numzValues+iz)=0;

		for(iconfig=0;iconfig<numcombine;iconfig++)
		{
			vp = (configs+iconfig)->data;
			*(data->voltage+iz)+=*(vp->voltage+iz);
			*(data->field+iz)+=*(vp->field+iz);
			*(data->charge+iz)+=*(vp->charge+iz);
			*(data->ionicCharge+iz)+=*(vp->ionicCharge+iz);
			*(data->electronCharge+iz)+=*(vp->electronCharge+iz);
			*(data->holeCharge+iz)+=*(vp->holeCharge+iz);
			for(iele=0;iele<numIons;iele++)
				*(data->ionsPerLayer+iele*numzValues+iz)+=*(vp->ionsPerLayer+iele*numzValues+iz);

		}
		 *(data->voltage+iz)/=numcombine;
		 *(data->field+iz)/=numcombine;
		 *(data->charge+iz)/=numcombine;
		 *(data->ionicCharge+iz)/=numcombine;
		 *(data->electronCharge+iz)/=numcombine;
		 *(data->holeCharge+iz)/=numcombine;
		// printf("Average charge at height %d is %g\n",iz,*(data->charge+iz)); 
		for(iele=0;iele<numIons;iele++)
			*(data->ionsPerLayer+iele*numzValues+iz)/=numcombine;
	}
	data->polarization=0;
	data->ionicPolarization=0;
	data->debyePolarization=0;
	for(iconfig=0;iconfig<numcombine;iconfig++)
	{
		vp = (configs+iconfig)->data;
		data->polarization+=vp->polarization;
		data->ionicPolarization+=vp->ionicPolarization;
		data->debyePolarization+=vp->debyePolarization;
	}
	data->polarization/=numcombine;
	data->ionicPolarization/=numcombine;
	data->debyePolarization/=numcombine;
	*(outconfig->seriesData) = data->polarization;	
	*(outconfig->seriesData+1) = data->ionicPolarization;	
	*(outconfig->seriesData+2) = data->debyePolarization;	
	
}

void es_setup(double *newzValues, int newnumzValues, char *newionNames,int newnumIons,double *newionCharge, void (*traj_generator)(Trajectory *))
{
	zValues = newzValues;
	numzValues = newnumzValues;
	ionNames = newionNames;
	numIons = newnumIons;
	ionCharge = newionCharge;
	
	printf("%d z values, %d charged elements\n",numzValues,numIons);
	LD_setup(es_Energy,VP_saveEveryOther,VP_allocate,VP_free,VP_combine, VP_combineWeighted, traj_generator, sizeof(VoltageProfile) );
	
	if(freeElectronDensity==0)
	screeningLength[0]=1e10;
	else
	screeningLength[0] = sqrt(dielectricConstant*permittivity*getTemperature()/freeElectronDensity);	
	if(freeHoleDensity==0)
	screeningLength[1]=1e10;
	else
	screeningLength[1] = sqrt(dielectricConstant*permittivity*getTemperature()/freeHoleDensity);	
	
	printf("Electron Screening length %g Hole Screening Length %g\n",screeningLength[0],screeningLength[1]);
	
}

void es_Trajectory(Trajectory *traj)
{
	int numSteps = getNumSteps();
    traj->numExternalConditions=1;//voltage control
    traj->externalConditions=malloc(numSteps * sizeof(double));
    //traj->enthalpyScalars=malloc(traj->numSteps * sizeof(double));
    sourceControl_Trajectory(traj);
    /*
    int iState,iStep;
    int usedStepsPerSwitch;
    if(stepsPerSwitch > traj->numSteps-balanceSteps)
	usedStepsPerSwitch = traj->numSteps-balanceSteps;
	else
	usedStepsPerSwitch=stepsPerSwitch;
    
    //This counts the first step which is ions only, and the final extra step which is a buffer
	traj->numEnthalpyBreakpoints = (traj->numSteps-balanceSteps)/usedStepsPerSwitch+2;
	printf("This trajectory has %d enthalpy breakpoints \n",traj->numEnthalpyBreakpoints);
	traj->enthalpyStates=malloc(traj->numEnthalpyBreakpoints * sizeof(int));
	traj->enthalpyStateBreakpoints=malloc(traj->numEnthalpyBreakpoints * sizeof(int));		
	traj->currentEnthalpyState=ION_ONLY;
	*(traj->enthalpyStates)=ION_ONLY;
	*(traj->enthalpyStateBreakpoints) = 0;
	*/

	traj->seriesData=malloc(3*numSteps*sizeof(double));	
	traj->numSeriesData=3;
	traj->seriesDataNames = malloc(3*sizeof(char *));
	*(traj->seriesDataNames) = "pol.csv";
	*(traj->seriesDataNames+1) = "ionPol.csv";
	*(traj->seriesDataNames+2) = "debyePol.csv";
		
	/*	
	for(iState=1;iState<traj->numEnthalpyBreakpoints;iState++)
	{
			
		if(iState%2 == 1)
			*(traj->enthalpyStates+iState)=ION_FORWARD_FIELD;
		if(iState%2 == 0)
			*(traj->enthalpyStates+iState)=ION_REVERSE_FIELD;

		*(traj->enthalpyStateBreakpoints+iState) = usedStepsPerSwitch*(iState-1)+balanceSteps;
		}
	
	double phase;
	for(iState=0;iState<traj->numEnthalpyBreakpoints-1;iState++)
	{
		for(iStep=*(traj->enthalpyStateBreakpoints+iState);iStep<*(traj->enthalpyStateBreakpoints+iState+1);iStep++)
		{
			 switch(vMode)
			{
			case VOLTAGE_SINE_WAVE:
			phase = (double)(iStep-*(traj->enthalpyStateBreakpoints+iState))/(*(traj->enthalpyStateBreakpoints+iState+1)-*(traj->enthalpyStateBreakpoints+iState))*M_PI;		
			break;
		
			case VOLTAGE_SQUARE_WAVE:
			phase=M_PI/2.0;
			break;
			}

			switch(*(traj->enthalpyStates+iState))
			{
				case ION_ONLY:
				*(traj->enthalpyScalars+iStep)=0;
				break;
		
				case ION_FORWARD_FIELD:
				*(traj->enthalpyScalars+iStep) = voltageOffset-voltageModulation*sin(phase);
				break;
		
				case ION_REVERSE_FIELD:
				*(traj->enthalpyScalars+iStep) = voltageOffset+voltageModulation*sin(phase);
				break;
			}
			//printf("Voltage on step %d is %g\n",iStep,*(traj->enthalpyScalars+iStep));
		}
	}*/
}
