/*
 * chargedcationhalidevacancies.c
 * 
 * Copyright 2019 Anthony Ruth <aruth@heesolar.com>
 * Simulates the movement of FA and I vacancies in a perovskite lattice
 * considering ionic drift and diffusion. 
 * 
 */

#include "perovskite.h"
#include "electrostatics.h"
char *interfaceElement = "Cu";
char dir[1000];
double FAVacancyRatio,IVacancyRatio;
double IHopEnergy, FAHopEnergy;
double latticeConstant = 6.3620;//DOI:10.1021/acs.jpclett.5b01432
int	sizex,sizey,sizez;
double defectCharge;

void ch_trajectory(Trajectory *traj);
crystal *ch_crystal();
//void ch_loadSettings(char *filename);
//void ch_saveSettings();
void ch_registerSettings();
void setup();


void ch_trajectory(Trajectory *traj)
{
     traj->crys = ch_crystal(traj);
     
     es_Trajectory(traj);
     
    LatticeDynamics *LD = traj->LD = malloc(sizeof(LatticeDynamics));
	
	if(FAVacancyRatio !=0 && IVacancyRatio !=0)
	{
	LD->numHopPairs = 2;
	LD->hopPairs = crys_elementString(4,VX,"I",VA,FA);	
	LD->hopPairEnergies = malloc(2*sizeof(double));
	
	*(LD->hopPairEnergies) = IHopEnergy; //I
	*(LD->hopPairEnergies+1) = FAHopEnergy; //FA
	}
	if(FAVacancyRatio ==0)
	{
	LD->numHopPairs = 1;
	LD->hopPairs = crys_elementString(2,VX,"I");	
	LD->hopPairEnergies = malloc(1*sizeof(double));
	
	*(LD->hopPairEnergies) = IHopEnergy; //I
	}
	if(IVacancyRatio ==0)
	{
	LD->numHopPairs = 1;
	LD->hopPairs = crys_elementString(2,VA,FA);	
	LD->hopPairEnergies = malloc(1*sizeof(double));
	
	*(LD->hopPairEnergies) = FAHopEnergy; //I
	}

	//This could be filled in from the hoppairs
	(traj->nnd)->nndistance=0.75*latticeConstant;
	(traj->nnd)->elementList=crys_elementString(5,"Pb","I",VA,VX,FA);
	(traj->nnd)->numEle=5;
	(traj->nnd+1)->nndistance=1.01*latticeConstant;
	(traj->nnd+1)->elementList=crys_elementString(2,VA,FA);
	(traj->nnd+1)->numEle=2;
	traj->numnnds=2;

}


crystal *ch_crystal()
{
	struct crystal *perovcrys,*crys;


    //Setup Crystal
	struct crystal *FAPbI = perovskite_newCrys(FA,"Pb","I","I","I",latticeConstant);
	perovcrys = crys_multiply(FAPbI,sizex,sizey,sizez);
	
	double aluminumLatConst = 4.05;
	int numHalideVacancies = IVacancyRatio*(perovcrys->totalAtoms)*3.0/5.0;
	int numCationVacancies = FAVacancyRatio*(perovcrys->totalAtoms)*1.0/5.0;
	double perovskiteCharge = chargeBalance(perovcrys,defectCharge*(numHalideVacancies-numCationVacancies),0);
	
	struct crystal *interface = crys_simpleCubic(interfaceElement,aluminumLatConst);
	int numxInterfaceIons = 1 + sizex*latticeConstant/aluminumLatConst;
	int numyInterfaceIons = 1 + sizey*latticeConstant/aluminumLatConst;
	crystal *interfacePlane = crys_multiply(interface,numxInterfaceIons,numyInterfaceIons,1);
	double alvector[3];
	alvector[0] = 0;
	alvector[1] = 0;
	
	//This convention for where the interface is located matches the direction of ION_FORWARD_FIELD in electrostatics.c
	if(perovskiteCharge < 0)
	alvector[2] = latticeConstant*sizez;
	else
	alvector[2] = latticeConstant*(-0.5);
	crys=crys_combine(perovcrys,interfacePlane,alvector,0);
	crys_free(perovcrys);
	crys_free(interfacePlane);
	
	printf("%d total atoms. Replacing %d atoms by vacancies actual concentration %f \n",crys->totalAtoms,numHalideVacancies,(double)numHalideVacancies/(crys->totalAtoms*3.0/5.0));
	crys_replaceRandomAtoms(crys,numHalideVacancies,crys_elementString(1,"I"),1,VX);

	printf("%d total atoms. Replacing %d atoms by vacancies actual concentration %f \n",crys->totalAtoms,numCationVacancies,(double)numCationVacancies/(crys->totalAtoms*1.0/5.0));
	crys_replaceRandomAtoms(crys,numCationVacancies,crys_elementString(1,FA),1,VA);

	//This eliminates periodic boundary conditions in the z direction.
	*(crys->latticeVectors+8) *= 2;

    crys_printElements(crys);

    
    return crys;
}

void ch_registerSettings()
{
	es_registerSettings();
	registerDouble(&FAVacancyRatio,"FAVacancyRatio",0);
	registerDouble(&IVacancyRatio,"IVacancyRatio",0);
	registerDouble(&FAHopEnergy,"FAHopEnergy",0.25);
	registerDouble(&IHopEnergy,"IHopEnergy",0.25);
	registerDouble(&defectCharge,"defectCharge",1);
	registerInt(&sizex,"sizex",4);
	registerInt(&sizey,"sizey",4);
	registerInt(&sizez,"sizez",4);
	
}
/*
void ch_loadSettings(char *filename)
{
    printf("Loading settings from file %s\n",filename);
	
	FILE *settingsfile = fopen(filename,"r");
    loadpmcSettings(settingsfile);
    es_loadSettings(settingsfile);
    
	dir = readString(settingsfile,"dir");    
    readDouble(settingsfile,"FAVacancyRatio",&(FAVacancyRatio));
    readDouble(settingsfile,"IVacancyRatio",&(IVacancyRatio));
	readDouble(settingsfile,"FAHopEnergy",&(FAHopEnergy));
	readDouble(settingsfile,"IHopEnergy",&(IHopEnergy));
	readDouble(settingsfile,"defectCharge",&(defectCharge));
    readInt(settingsfile,"sizex",&(sizex));
	readInt(settingsfile,"sizey",&(sizey));
	readInt(settingsfile,"sizez",&(sizez));	

    fclose(settingsfile);
}

void ch_saveSettings()
{
    printf("Dir is %s\n",dir);
    mkdir2(dir);
    char *filename = strcat2(dir,"settings");
    printf("Saving settings to %s\n",filename);
    FILE *settingsfile = fopen(filename,"w");

    savepmcSettings(settingsfile);
    es_saveSettings(settingsfile);
    
	fprintf(settingsfile,"IHopEnergy = %g\n",IHopEnergy);
	fprintf(settingsfile,"FAHopEnergy = %g\n",FAHopEnergy);
	fprintf(settingsfile,"FAVacancyRatio = %g\n",FAVacancyRatio);
	fprintf(settingsfile,"IVacancyRatio = %g\n",IVacancyRatio);
	fprintf(settingsfile,"defectCharge = %g\n",defectCharge);
	fprintf(settingsfile,"sizex = %d\n",sizex);
	fprintf(settingsfile,"sizey = %d\n",sizey);
	fprintf(settingsfile,"sizez = %d\n",sizez);
	    
    fclose(settingsfile);
	
	char command[1000];
	sprintf(command,"cat %s\n",strcat2(dir,"settings"));
	system(command);
}*/

void setup()
{
	
    char *ionNames; 	
	int numIons;
	double *ionCharge = malloc(3*sizeof(double));

	int iInterface;
	if(FAVacancyRatio !=0 && IVacancyRatio !=0)
	{
	numIons = 3;
	ionNames = crys_elementString(3,VX,VA,interfaceElement);
	*ionCharge = defectCharge; *(ionCharge+1) = -defectCharge;
	iInterface=2;
	}
	if(FAVacancyRatio ==0)
	{
	numIons = 2;
	ionNames = crys_elementString(2,VX,interfaceElement);
	*ionCharge = defectCharge;
	iInterface=1;
	}
	if(IVacancyRatio ==0)
	{
	numIons = 2;
	ionNames = crys_elementString(2,VA,interfaceElement);
	*ionCharge = -defectCharge;
	iInterface=1;
	}
	
	
	crystal *crys = ch_crystal();
	double ionicCharge = defectCharge * (crys_elementCount(crys,VX)-crys_elementCount(crys,VA));
	double perovskiteCharge = chargeBalance(crys,ionicCharge,0);
	//double interfaceIonCharge = -perovskiteCharge/crys_elementCount(crys,interfaceElement);
	double interfaceIonCharge = 0;
	*(ionCharge+iInterface) = interfaceIonCharge;
	printf("There is an ionic charge of %g a perovskite charge of %g, and each %s will have a charge of %g\n",ionicCharge,perovskiteCharge,interfaceElement,interfaceIonCharge);
	
	double *zValues = malloc(crys->totalAtoms*sizeof(double));
	int numzValues;
	crys_zValues(crys,zValues,&numzValues);	
	//printmatrix(zValues,numzValues,1,"");
	
	es_setup(zValues, numzValues, ionNames,numIons,ionCharge,ch_trajectory);
}

int main(int argc, char **argv)
{
	allocateSettings();
	ch_registerSettings();
	
	
	FILE *infile = fopen(argv[1],"r");
	loadSettings(infile);
	fclose(infile);
	mkdir2(dir);
	 
	char outfileName[1000];
	sprintf(outfileName,"%s/settings",dir);
	FILE *outfile = fopen(outfileName,"w");
	saveSettings(outfile);
	fclose(outfile);
		    
    setup();  
    simulateTrajectories();
		
	return 0;
}


