

#include "polarization.h"

char *ionNames;
double *ionCharge;	
int numIons;

char *dipoleNames;
double *dipolePolarizability;	
int numDipoles;

//double fieldstrength; // V/angstrom
double voltageOffset;
double voltageModulation;
double *zvalues;
int numzValues;
double screeningLength=1e300;
double dielectricConstant=1;
int vMode;
int stepsPerSwitch;//switches from forward bias to reverse bias with this frequency.
int balanceSteps;
double	temperature;//In units of eV There could be seperate electronic and ionic temperatures
//These are used by pz_Energy
int numID, numI, numD;
void pz_loadSettings(FILE *settingsfile)
{
	//readDouble(settingsfile,"fieldstrength",&(fieldstrength));
	readDouble(settingsfile,"voltageOffset",&(voltageOffset));
	readDouble(settingsfile,"voltageModulation",&(voltageModulation));
	readDouble(settingsfile,"screeningLength",&(screeningLength));
	readDouble(settingsfile,"dielectricConstant",&(dielectricConstant));
	readDouble(settingsfile,"temperature",&(temperature));
	readInt(settingsfile,"stepsPerSwitch",&(stepsPerSwitch));
	readInt(settingsfile,"balanceSteps",&(balanceSteps));
	
	
		
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

}

void pz_saveSettings(FILE *settingsfile)
{
	//fprintf(settingsfile,"fieldstrength = %g\n",fieldstrength);
	fprintf(settingsfile,"voltageOffset = %g\n",voltageOffset);
	fprintf(settingsfile,"voltageModulation = %g\n",voltageModulation);
	fprintf(settingsfile,"screeningLength = %g\n",screeningLength);
	fprintf(settingsfile,"dielectricConstant = %g\n",dielectricConstant);
	fprintf(settingsfile,"stepsPerSwitch = %d\n",stepsPerSwitch);
	fprintf(settingsfile,"balanceSteps = %d\n",balanceSteps);
	fprintf(settingsfile,"temperature = %g\n",temperature);
	switch(vMode)
	{
		case VOLTAGE_SQUARE_WAVE:
		fprintf(settingsfile,"vMode = square\n");
		break;
		case VOLTAGE_SINE_WAVE:
		fprintf(settingsfile,"vMode = sine\n");
		break;
	}

}

void FP_allocate(FieldProfile *FP)
{
	Trajectory *traj = traj_getSettings();
	crystal *crys = traj->crys;
	
    FP->layerCharge = calloc(numzValues,sizeof(double));  
    FP->layerVoltage = calloc(numzValues,sizeof(double));  
    FP->ionsPerLayer = calloc(numzValues*numIons,sizeof(double));  
    FP->dipolesPerLayer = calloc(numzValues*numDipoles,sizeof(double));  
    FP->layerMoment = calloc(numzValues * 3,sizeof(double));  
    FP->layerField = calloc(numzValues * 3,sizeof(double));  
    
    int iz, iIon,iDipole;
	
	numI=0;
	for(iIon=0;iIon<numIons;iIon++)
	numI += crys_elementCount(crys,ionNames+iIon*namelength);
	
	numD=0;
	for(iDipole=0;iDipole<numDipoles;iDipole++)
	numD += crys_elementCount(crys,dipoleNames+iDipole*namelength);

	numID=numI+numD;
	//printf("There are %d ions and %d dipoles total\n",numI,numD);

	FP->field=calloc(numID*3,sizeof(double));
	FP->moment=calloc(numD*3,sizeof(double));
	FP->voltage=calloc(numID,sizeof(double));

}

void FP_free(FieldProfile *FP)
{
    free(FP->layerCharge);  
    free(FP->layerMoment);  
    free(FP->layerField);  
    free(FP->layerVoltage);  
    free(FP->ionsPerLayer);
    free(FP->dipolesPerLayer);
    free(FP->field);
	free(FP->moment);
	free(FP->voltage);
}


//This calculates the potential, field, and expected moment of every ion and dipole in the system.
//(1) the external field is applied and the potential due the the external field is calculated
//(2) The field and potential due to the ions is added
//(3) the field and potential due to the dipoles is added
//(4) the total energy is calculated
//(5) the dipoles are given an expected polarization based on the field at their location.  
//(6) the process repeats at (1) until the energy difference between two iterations is below some criterion
//(7) properties of the layers are calculated
double pz_Energy(Configuration *config)
{
	//Calculates the electrostatic energy of charged particles in an electric field.

	FieldProfile *FP = config->data;
	crystal *crys = config->crys;

	double *voltage = FP->voltage;
	double *field = FP->field;
	double *moment = FP->moment;
	
	double *layerCharge=FP->layerCharge;
	double *layerVoltage=FP->layerVoltage;
	double *layerField=FP->layerField;
	double *layerMoment=FP->layerMoment;
	
	double *ionsPerLayer=FP->ionsPerLayer;
	double *dipolesPerLayer=FP->dipolesPerLayer;

	int iIon1, iIon2;
	int iDipole1, iDipole2;
	int iDim;
	int iatom1,iatom2;
	int ele1start,ele1end;
	int ele2start,ele2end;
	int ele1count,ele2count;
	int iIonDipole1, iIonDipole2;
	//These record how far the particular species is from the start
	int dipoleOffset,ionDipoleOffset;
	
	double r[3],dist,dp,magField;//r is the vector between two points, d is the distance.
	
	double width = sqrt(*(crys->latticeVectors) * *(crys->latticeVectors+4));
	double height = *(crys->latticeVectors+8)/2;
	
		
	double expectedPolarization;	
	double convergenceCriterion = 1e-5,previousEnergy=config->energy,energy;
	double convergence = 1;
	int numIterations=0;
	
	double appliedVoltage=config->enthalpyScalar;
	int iz;
		
	double fieldStrength = -appliedVoltage/ height/dielectricConstant;
	
	while(convergence>convergenceCriterion)
	{
	//printf("(1) the external field is applied and the potential due the the external field is calculated\n");
	ionDipoleOffset=0;
	for(iIon1=0;iIon1<numIons;iIon1++)
	{
		ele1start = crys_elementOffset(crys,ionNames+iIon1*namelength);
		ele1count = crys_elementCount(crys,ionNames+iIon1*namelength);
		
		ele1end = ele1start + ele1count;
		
		for(iatom1=ele1start;iatom1<ele1end;iatom1++)
		{
			iIonDipole1 = iatom1-ele1start+ionDipoleOffset;
			//printf("iIonDipole1 %d\n", iIonDipole1);
			*(field + 3*iIonDipole1+0) = 0;
			*(field + 3*iIonDipole1+1) = 0;
			*(field + 3*iIonDipole1+2) = fieldStrength;
			//*(voltage + iIonDipole1) = -fieldStrength * dielectricConstant * *(crys->positions+iatom1*3+2);
			*(voltage + iIonDipole1) = -fieldStrength * *(crys->positions+iatom1*3+2);
		}
		ionDipoleOffset += ele1count;
	}
	
	for(iDipole1=0;iDipole1<numDipoles;iDipole1++)
	{
		ele1start = crys_elementOffset(crys,dipoleNames+iDipole1*namelength);
		ele1count = crys_elementCount(crys,dipoleNames+iDipole1*namelength);
		ele1end = ele1start + ele1count;
		//printf("iDipole1 %d, name %s\n",iDipole1, dipoleNames+iDipole1*namelength);
		//printf("ele1start %d, ele1count %d, ele1end %d\n",ele1start,ele1count,ele1end);
		for(iatom1=ele1start;iatom1<ele1end;iatom1++)
		{
			iIonDipole1 = iatom1-ele1start+ionDipoleOffset;
			*(field + 3*iIonDipole1+0) = 0;
			*(field + 3*iIonDipole1+1) = 0;
			*(field + 3*iIonDipole1+2) = fieldStrength;
			//*(voltage + iIonDipole1) = -fieldStrength * dielectricConstant * *(crys->positions+iatom1*3+2);
			*(voltage + iIonDipole1) = -fieldStrength * *(crys->positions+iatom1*3+2);
		}
		ionDipoleOffset += ele1count;
	}
	//(2) The field and potential due to the ions is added
	//printf("(2) The field and potential due to the ions is added\n");
	for(iIon1=0;iIon1<numIons;iIon1++)
	{
		ele1start = crys_elementOffset(crys,ionNames+iIon1*namelength);
		ele1count = crys_elementCount(crys,ionNames+iIon1*namelength);
		ele1end = ele1start + ele1count;
		for(iatom1=ele1start;iatom1<ele1end;iatom1++)
		{
			ionDipoleOffset=0;
			for(iIon2=0;iIon2<numIons;iIon2++)
			{
			//printf("Applying ion %s to ion %s\n",(ionNames+namelength*iIon1),(ionNames+namelength*iIon2));
			ele2start = crys_elementOffset(crys,ionNames+iIon2*namelength);
			ele2count = crys_elementCount(crys,ionNames+iIon2*namelength);
			ele2end = ele2start + ele2count;
				for(iatom2=ele2start;iatom2<ele2end;iatom2++)
				{
				if(iatom1==iatom2)
				continue;
				dist=crys_atomVector(crys,iatom1,iatom2,r);
				iIonDipole2 = iatom2-ele2start+ionDipoleOffset;
				*(voltage+iIonDipole2) += coulombconstant* *(ionCharge+iIon1) * exp(-dist/screeningLength)/(dist*dielectricConstant);
					for(iDim=0;iDim<3;iDim++)
					*(field+iIonDipole2*3+iDim) += *(r+iDim)/pow(dist,3) *coulombconstant* *(ionCharge+iIon1) * exp(-dist/screeningLength)/(dielectricConstant);
				}
			ionDipoleOffset += ele2count;
			//printf("IonDipoleOffset %d\n",ionDipoleOffset);	
			}
			
			for(iDipole2=0;iDipole2<numDipoles;iDipole2++)
			{
			//printf("Applying ion %s to dipole %s\n",(ionNames+namelength*iIon1),(dipoleNames+namelength*iDipole2));
			ele2start = crys_elementOffset(crys,dipoleNames+iDipole2*namelength);
			ele2count = crys_elementCount(crys,dipoleNames+iDipole2*namelength);
			ele2end = ele2start + ele2count;
				for(iatom2=ele2start;iatom2<ele2end;iatom2++)
				{
				dist=crys_atomVector(crys,iatom1,iatom2,r);
				iIonDipole2 = iatom2-ele2start+ionDipoleOffset;
				*(voltage+iIonDipole2) += coulombconstant* *(ionCharge+iIon1) * exp(-dist/screeningLength)/(dist*dielectricConstant);
					for(iDim=0;iDim<3;iDim++)
					*(field+iIonDipole2*3+iDim) += *(r+iDim)/pow(dist,3) * coulombconstant* *(ionCharge+iIon1) * exp(-dist/screeningLength)/(dielectricConstant);
				}
			ionDipoleOffset += ele2count;
			//printf("IonDipoleOffset %d\n",ionDipoleOffset);	
			}
		}
		
	}
	//(3) the field and potential due to the dipoles is added
	//printf("(3) the field and potential due to the dipoles is added\n");
	dipoleOffset = 0;
	for(iDipole1=0;iDipole1<numDipoles;iDipole1++)
	{
		ele1start = crys_elementOffset(crys,dipoleNames+iDipole1*namelength);
		ele1count = crys_elementCount(crys,dipoleNames+iDipole1*namelength);
		ele1end = ele1start + ele1count;
		for(iatom1=ele1start;iatom1<ele1end;iatom1++)
		{
			iIonDipole1 = iatom1-ele1start+dipoleOffset;
			ionDipoleOffset=0;
			for(iIon2=0;iIon2<numIons;iIon2++)
			{
			ele2start = crys_elementOffset(crys,ionNames+iIon2*namelength);
			ele2count = crys_elementCount(crys,ionNames+iIon2*namelength);
			ele2end = ele2start + ele2count;
				for(iatom2=ele2start;iatom2<ele2end;iatom2++)
				{
				dist=crys_atomVector(crys,iatom1,iatom2,r);
				iIonDipole2 = iatom2-ele2start+ionDipoleOffset;
				dp=dotProduct(r,(moment+iIonDipole1*3),3);
				//*(voltage+iIonDipole2) += dp/pow(dist,3) * coulombconstant* *(dipoleMoment+iDipole1) * exp(-dist/screeningLength)/dielectricConstant;
				*(voltage+iIonDipole2) += dp/pow(dist,3) * coulombconstant * exp(-dist/screeningLength)/dielectricConstant;
					for(iDim=0;iDim<3;iDim++)
					*(field+iIonDipole2*3+iDim) += *(moment+iIonDipole1*3+iDim)/pow(dist,3) *coulombconstant* exp(-dist/screeningLength)/dielectricConstant;
				}
			ionDipoleOffset += ele2count;	
			}
			
			for(iDipole2=0;iDipole2<numDipoles;iDipole2++)
			{
			ele2start = crys_elementOffset(crys,dipoleNames+iDipole2*namelength);
			ele2count = crys_elementCount(crys,dipoleNames+iDipole2*namelength);
			ele2end = ele2start + ele2count;
				for(iatom2=ele2start;iatom2<ele2end;iatom2++)
				{
				if(iatom1==iatom2)
				continue;	
				dist=crys_atomVector(crys,iatom1,iatom2,r);
				iIonDipole2 = iatom2-ele2start+ionDipoleOffset;
				dp=dotProduct(r,(moment+iIonDipole1*3),3);
				//*(voltage+iIonDipole2) += dp/pow(dist,3) * coulombconstant* *(dipoleMoment+iDipole1) * exp(-dist/screeningLength)/dielectricConstant;
				*(voltage+iIonDipole2) += dp/pow(dist,3) * coulombconstant * exp(-dist/screeningLength)/dielectricConstant;
					for(iDim=0;iDim<3;iDim++)
					*(field+iIonDipole2*3+iDim) += *(moment+iIonDipole1*3+iDim)/pow(dist,3) *coulombconstant* exp(-dist/screeningLength)/(dielectricConstant);
				}
			ionDipoleOffset += ele2count;	
			}
		}
		dipoleOffset += ele1end-ele1start;	
	}
	//printf("(4) the total energy is calculated\n");
	//printf("(5) the dipoles are given an expected polarization based on the field at their location. \n");
	energy=0;
	ionDipoleOffset=0;
	for(iIon1=0;iIon1<numIons;iIon1++)
	{
		ele1start = crys_elementOffset(crys,ionNames+iIon1*namelength);
		ele1count = crys_elementCount(crys,ionNames+iIon1*namelength);
		ele1end = ele1start + ele1count;
		for(iatom1=ele1start;iatom1<ele1end;iatom1++)
		{
			iIonDipole1 = iatom1-ele1start+ionDipoleOffset;
			energy += *(voltage+iIonDipole1) * *(ionCharge+iIon1);
		}
		ionDipoleOffset += ele1count;
	}
	dipoleOffset=0;
	
	for(iDipole1=0;iDipole1<numDipoles;iDipole1++)
	{
		ele1start = crys_elementOffset(crys,dipoleNames+iDipole1*namelength);
		ele1count = crys_elementCount(crys,dipoleNames+iDipole1*namelength);
		ele1end = ele1start + ele1count;
		//printf("ele1start %d, ele1count %d, ele1end %d\n",ele1start,ele1count,ele1end);
		for(iatom1=ele1start;iatom1<ele1end;iatom1++)
		{
			iIonDipole1 = iatom1-ele1start+ionDipoleOffset;
			iIonDipole2 = iatom1-ele1start+dipoleOffset;
			//printf("Field is %g %g %g\n",*(field+3*iIonDipole1),*(field+3*iIonDipole1+1),*(field+3*iIonDipole1+2));
			
			magField = magnitude(field+3*iIonDipole1,3);
			
			//expectedPolarization = -temperature/magField+*(dipoleMoment+iDipole1)*1/tanh(*(dipoleMoment+iDipole1)*magField/temperature);
			expectedPolarization = magField * *(dipolePolarizability+iDipole1);
			for(iDim=0;iDim<3;iDim++)
				*(moment+3*iIonDipole2+iDim) = *(field+3*iIonDipole1+iDim)/dist * expectedPolarization;
				
			energy += dotProduct((field+3*iIonDipole1),(moment+3*iIonDipole2),3);
		}
		ionDipoleOffset += ele1count;
		dipoleOffset += ele1count;
	}
	//without dipoles there is no need to solve this iteratively.
	if(numDipoles==0)
	convergence = convergenceCriterion;
	else
	convergence = fabs((previousEnergy-energy)/energy);
	//printf("iteration %d old energy %g new energy %g convergence %g\n",numIterations,previousEnergy,energy,convergence);
	previousEnergy=energy;
	numIterations++;
	if(numIterations==ITERATION_LIMIT)
	{
		printf("Reached %d iterations without converging\n",ITERATION_LIMIT);
		break;
	}
	}
	//printf("Convergence was reached after %d iterations\n",numIterations);
	//printf("(7) properties of the layers are calculated\n");
	int numIonsEachZ[numzValues],numDipolesEachZ[numzValues];
	for(iz=0;iz<numzValues;iz++)
	{
		*(layerCharge+iz) = 0;
		*(layerVoltage+iz) = 0;
		numIonsEachZ[iz]=0;
		numDipolesEachZ[iz]=0;
		for(iDim=0;iDim<3;iDim++)
		{
			*(layerField+3*iz+iDim) = 0;
			*(layerMoment+3*iz+iDim) = 0;
		}
		for(iIon1=0;iIon1<numIons;iIon1++)
			*(ionsPerLayer+iz*numIons+iIon1)=0;
		for(iDipole1=0;iDipole1<numDipoles;iDipole1++)
			*(dipolesPerLayer+iz*numDipoles+iDipole1)=0;
	}
	ionDipoleOffset=0;
	for(iIon1=0;iIon1<numIons;iIon1++)
	{
		ele1start = crys_elementOffset(crys,ionNames+iIon1*namelength);
		ele1count = crys_elementCount(crys,ionNames+iIon1*namelength);
		ele1end = ele1start + ele1count;
		for(iatom1=ele1start;iatom1<ele1end;iatom1++)
		{
			iIonDipole1 = iatom1-ele1start+ionDipoleOffset;
			//printf("iIonDipole1 %d\n", iIonDipole1);
			iz = binarysearchDouble(zvalues,*(crys->positions+3*iatom1+2),0,numzValues);
			*(layerCharge+iz) += *(ionCharge+iIon1);
			//printf("Voltage of iIonDipole1 %d is %g\n",iIonDipole1,*(voltage+iIonDipole1));
			*(layerVoltage+iz) += *(voltage+iIonDipole1);
			for(iDim=0;iDim<3;iDim++)
			{
			//printf("Field of iIonDipole1 %d in direction %d is %g\n",iIonDipole1,iDim,*(field+3*iIonDipole1+iDim));
			*(layerField+3*iz+iDim) += *(field+3*iIonDipole1+iDim);
			}
			(*(ionsPerLayer+iz*numIons+iIon1))++;
			(*(numIonsEachZ+iz))++;
		}
		ionDipoleOffset += ele1count;
	}
	dipoleOffset=0;
	for(iDipole1=0;iDipole1<numDipoles;iDipole1++)
	{
		ele1start = crys_elementOffset(crys,dipoleNames+iDipole1*namelength);
		ele1count = crys_elementCount(crys,dipoleNames+iDipole1*namelength);
		ele1end = ele1start + ele1count;
		for(iatom1=ele1start;iatom1<ele1end;iatom1++)
		{
			iIonDipole1 = iatom1-ele1start+ionDipoleOffset;
			iIonDipole2 = iatom1-ele1start+dipoleOffset;
			//printf("iIonDipole1 %d\n", iIonDipole1);
			//printf("iIonDipole2 %d\n", iIonDipole2);
			iz = binarysearchDouble(zvalues,*(crys->positions+3*iatom1+2),0,numzValues);
			//printf("Voltage of iIonDipole1 %d is %g\n",iIonDipole1,*(voltage+iIonDipole1));
			*(layerVoltage+iz) += *(voltage+iIonDipole1);
			//printf("Voltage of layer %d is %g\n",iz,*(layerVoltage+iz));
			for(iDim=0;iDim<3;iDim++)
			{
			*(layerField+3*iz+iDim) += *(field+3*iIonDipole1+iDim);
			*(layerMoment+3*iz+iDim) += *(moment+3*iIonDipole2+iDim);
			}
			(*(dipolesPerLayer+iz*numDipoles+iDipole1))++;
			(*(numDipolesEachZ+iz))++;
		}
		ionDipoleOffset += ele1count;
		dipoleOffset += ele1count;
	}
	FP->polarization =0;
	for(iz=0;iz<numzValues;iz++)
	{
		if(*(numIonsEachZ+iz)+*(numDipolesEachZ+iz) != 0)
			*(layerVoltage+iz) /= *(numIonsEachZ+iz) + *(numDipolesEachZ+iz);
		else
			*(layerVoltage+iz) =NAN;
			
		for(iDim=0;iDim<3;iDim++)
		{
			if(*(numIonsEachZ+iz)+*(numDipolesEachZ+iz) != 0)
				*(layerField+3*iz+iDim) /= *(numIonsEachZ+iz) + *(numDipolesEachZ+iz);
			else
				*(layerField+3*iz+iDim) = NAN;
		}
		FP->polarization += *(layerCharge+iz) * *(zvalues+iz);
		FP->polarization += *(layerMoment+3*iz);
		
	}
	if(!(config->seriesData==NULL))
	*(config->seriesData)=FP->polarization;
	
	//printf("Energy is %g saving it to %d\n",energy,(int *)(config->energy));
	(config->energy) = energy;
	return energy;
}

void FP_save(Configuration *config)
{
    
	FILE *outfile = fopen(config->dataFileName,"w");
	FieldProfile *FP = config->data;

    //printf("VP_save %s Config is at %d data is at %d\n",config->dataFileName,(int)config,(int)data);
	int iz,iele;
	fprintf(outfile,"%g\n",FP->polarization);
	for(iz=0;iz<numzValues;iz++)
	{
	fprintf(outfile,"%g %g %g %g %g %g %g %g %g",*(zvalues+iz),*(FP->layerCharge+iz),*(FP->layerVoltage+iz),
	*(FP->layerField+3*iz),*(FP->layerField+3*iz+1),*(FP->layerField+3*iz+2),
	*(FP->layerMoment+3*iz),*(FP->layerMoment+3*iz+1),*(FP->layerMoment+3*iz+2));
		for(iele=0;iele<numIons;iele++)
			fprintf(outfile," %g",*(FP->ionsPerLayer+iz*numIons+iele));
		for(iele=0;iele<numDipoles;iele++)
			fprintf(outfile," %g",*(FP->dipolesPerLayer+iz*numDipoles+iele));
	
	fprintf(outfile,"\n");
	}
	fclose(outfile);
}
void FP_saveEveryOther(Configuration *config)
{
    
	FILE *outfile = fopen(config->dataFileName,"w");
	FieldProfile *FP = config->data;

    //printf("VP_save %s Config is at %d data is at %d\n",config->dataFileName,(int)config,(int)data);
	int iz,iele,idim;
	double zvalue,charge,voltage,field,moment;
	fprintf(outfile,"%g\n",FP->polarization);
	for(iz=0;iz<numzValues/2;iz++)
	{
	zvalue = avgWithNAN(*(zvalues+2*iz),*(zvalues+2*iz+1));
	voltage = avgWithNAN(*(FP->layerVoltage+2*iz),*(FP->layerVoltage+2*iz+1));
	charge = avgWithNAN(*(FP->layerCharge+2*iz),*(FP->layerCharge+2*iz+1));
	fprintf(outfile,"%g %g %g",zvalue,charge,voltage);
	//fprintf(outfile,"%g %g %g %g %g %g %g %g %g",(*(zvalues+2*iz)+*(zvalues+2*iz+1))/2,(*(FP->layerCharge+2*iz)+*(FP->layerCharge+2*iz+1))/2,(*(FP->layerVoltage+2*iz)+*(FP->layerVoltage+2*iz+1))/2,
	//(*(FP->layerField+3*2*iz)+*(FP->layerField+3*(2*iz+1)))/2,(*(FP->layerField+3*2*iz+1)+*(FP->layerField+3*(2*iz+1)+1))/2,(*(FP->layerField+3*2*iz+2)+*(FP->layerField+3*(2*iz+1)+2))/2,
	//(*(FP->layerMoment+3*2*iz)+*(FP->layerMoment+3*(2*iz+1)))/2,(*(FP->layerMoment+3*2*iz+1)+*(FP->layerMoment+3*(2*iz+1)+1))/2,(*(FP->layerMoment+3*2*iz+2)+*(FP->layerMoment+3*(2*iz+1)+2))/2);
		for(idim=0;idim<3;idim++)
		{
			field = avgWithNAN(*(FP->layerField+3*2*iz+idim),*(FP->layerField+3*(2*iz+1)+idim));
			fprintf(outfile," %g",field);
		}
		for(idim=0;idim<3;idim++)
		{
			moment = avgWithNAN(*(FP->layerMoment+3*2*iz+idim),*(FP->layerMoment+3*(2*iz+1)+idim));
			fprintf(outfile," %g",moment);
		}
		for(iele=0;iele<numIons;iele++)
			fprintf(outfile," %g",avgWithNAN(*(FP->ionsPerLayer+2*iz*numIons+iele),*(FP->ionsPerLayer+(2*iz+1)*numIons+iele)));
		for(iele=0;iele<numDipoles;iele++)
			fprintf(outfile," %g",avgWithNAN(*(FP->dipolesPerLayer+2*iz*numDipoles+iele),*(FP->dipolesPerLayer+(2*iz+1)*numDipoles+iele)));
	
	fprintf(outfile,"\n");
	}
	fclose(outfile);
}


void FP_combine(Configuration *configs, int numcombine,Configuration *outconfig)
{
	//printf("Combining %d voltage profiles into one\n",numcombine);
	FieldProfile *FP = outconfig->data;
	FieldProfile *fp;
	int iz,iconfig,iele, iDim;
	//Field and Voltage is only calculated where there are charges or dipoles. This is used to ignore empty entries
	int numField,numVoltage;
	for(iz=0;iz<numzValues;iz++)
	{
		 *(FP->layerVoltage+iz)=0;
		 *(FP->layerCharge+iz)=0;
		 for(iDim=0;iDim<3;iDim++)
		 {
			 *(FP->layerField+3*iz+iDim)=0;
			 *(FP->layerMoment+3*iz+iDim)=0;
		 }
		for(iele=0;iele<numIons;iele++)
			*(FP->ionsPerLayer+iz*numIons+iele)=0;
		for(iele=0;iele<numDipoles;iele++)
			*(FP->dipolesPerLayer+iz*numDipoles+iele)=0;
		
		numField=numVoltage=0;	
		for(iconfig=0;iconfig<numcombine;iconfig++)
		{
			fp = (configs+iconfig)->data;
			if(!isnan(*(fp->layerVoltage+iz)))
			{
			numVoltage++;
			*(FP->layerVoltage+iz)+=*(fp->layerVoltage+iz);
			}
			*(FP->layerCharge+iz)+=*(fp->layerCharge+iz);
			for(iDim=0;iDim<3;iDim++)
			{
				if(!isnan(*(fp->layerField+3*iz+iDim)))
				{
				*(FP->layerField+3*iz+iDim)+=*(fp->layerField+3*iz+iDim);
				numField++;
				}
				*(FP->layerMoment+3*iz+iDim)+=*(fp->layerMoment+3*iz+iDim);
			}
			for(iele=0;iele<numIons;iele++)
				*(FP->ionsPerLayer+iz*numIons+iele)+=*(fp->ionsPerLayer+iz*numIons+iele);
			for(iele=0;iele<numDipoles;iele++)
				*(FP->dipolesPerLayer+iz*numDipoles+iele)+=*(FP->dipolesPerLayer+iz*numDipoles+iele);

		}
		*(FP->layerVoltage+iz)/=numVoltage;
		*(FP->layerCharge+iz)/=numcombine;
		for(iDim=0;iDim<3;iDim++)
		{
			*(FP->layerField+3*iz+iDim)/=(numField/3);
			*(FP->layerMoment+3*iz+iDim)/=numcombine;
		}
		for(iele=0;iele<numIons;iele++)
			*(FP->ionsPerLayer+iz*numIons+iele)/=numcombine;
		for(iele=0;iele<numDipoles;iele++)
			*(FP->dipolesPerLayer+iz*numDipoles+iele)/=numcombine;
	}
	
	FP->polarization=0;
	for(iconfig=0;iconfig<numcombine;iconfig++)
	{
		fp = (configs+iconfig)->data;
		FP->polarization+=fp->polarization;
	}
	FP->polarization/=numcombine;
	*(outconfig->seriesData) = FP->polarization;	
	
}

void pz_setup(double *newzvalues, int newnumzvalues, char *newionNames,int newnumIons,double *newionCharge, char *newdipoleNames,int newnumDipoles,double *newdipolePolarizability, void (*traj_generator)(Trajectory *))
{
	zvalues = newzvalues;
	numzValues = newnumzvalues;
	ionNames = newionNames;
	numIons = newnumIons;
	
	ionCharge = newionCharge;
	dipoleNames = newdipoleNames;
	numDipoles = newnumDipoles;
	dipolePolarizability = newdipolePolarizability;
	
	printf("%d z values, %d charged elements %d dipoles\n",numzValues,numIons,numDipoles);
	int i;
	for(i=0;i<numIons;i++)
	printf("Ion %s charge %g\n",(ionNames+i*namelength),*(ionCharge+i));
	
	for(i=0;i<numDipoles;i++)
	printf("Dipole %s polarizability %g\n",(dipoleNames+i*namelength),*(dipolePolarizability+i));
	
	//LD_setup(pz_Energy,FP_save,FP_allocate,FP_free,FP_combine, traj_generator, sizeof(FieldProfile) );
	LD_setup(pz_Energy,FP_saveEveryOther,FP_allocate,FP_free,FP_combine, traj_generator, sizeof(FieldProfile) );

}

void pz_Trajectory(Trajectory *traj)
{
    
    traj->enthalpyScalars=malloc(traj->numSteps * sizeof(double));
    
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


	traj->seriesData=malloc(1*traj->numSteps*sizeof(double));	
	traj->numSeriesData=1;
	traj->seriesDataNames = malloc(1*sizeof(char *));
	*(traj->seriesDataNames) = "pol.csv";
		
		
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
	}
}
