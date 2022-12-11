/*
 * chargedvacanciespolarFA.c
 * 
 * Copyright 2019 Anthony Ruth <aruth@heesolar.com>
 * Simulates the movement of FA and I vacancies in a perovskite lattice
 * considering ionic drift and diffusion and the static dipole of FA. 
 * 
 */

#include "perovskite.h"
#include "polarization.h"

char *dir;
//Currently FAVacancyRatio is ignored and an equal number of FA and I vacancies are used
double FAVacancyRatio, IVacancyRatio, CsRatio;
double IHopEnergy, FAHopEnergy, CsHopEnergy;
double latticeConstant = 6.3620;//DOI:10.1021/acs.jpclett.5b01432
int	sizex,sizey,sizez;
double defectCharge;
double FAPolarizability,IPolarizability;

void cvp_trajectory(Trajectory *traj);
crystal *cvp_crystal();
void cvp_loadSettings(char *filename);
void cvp_saveSettings();
void setup();


void cvp_trajectory(Trajectory *traj)
{
     traj->crys = cvp_crystal(traj);
     
     pz_Trajectory(traj);
     
    LatticeDynamics *LD = traj->LD = malloc(sizeof(LatticeDynamics));
	
	LD->numHopPairs = 3;
	LD->hopPairs = crys_elementString(6,VX,"I",VA,FA,VA,"Cs");	
	LD->hopPairEnergies = malloc(3*sizeof(double));
	
	*(LD->hopPairEnergies) = IHopEnergy; //I
	*(LD->hopPairEnergies+1) = FAHopEnergy; //FA
	*(LD->hopPairEnergies+2) = CsHopEnergy; //Cs
	
	//This could be filled in from the hoppairs
	(traj->nnd)->nndistance=0.75*latticeConstant;
	(traj->nnd)->elementList=crys_elementString(6,"Pb","I",VA,VX,FA,"Cs");
	(traj->nnd)->numEle=6;
	(traj->nnd+1)->nndistance=1.01*latticeConstant;
	(traj->nnd+1)->elementList=crys_elementString(3,VA,FA,"Cs");
	(traj->nnd+1)->numEle=2;
	traj->numnnds=2;

}


crystal *cvp_crystal()
{
	struct crystal *crys;


    //Setup Crystal
	struct crystal *FAPbI = perovskite_newCrys(FA,"Pb","I","I","I",latticeConstant);
	crys = crys_multiply(FAPbI,sizex,sizey,sizez);
	crys_replaceRandomAtoms(crys,CsRatio*(sizex*sizey*sizez),FA,1,"Cs");
		
	//int numHalideVacancies = IVacancyRatio*(crys->totalAtoms)*3.0/5.0;
	int numHalideVacancies = IVacancyRatio*(crys->totalAtoms)*3.0/5.0;
	int numCationVacancies = FAVacancyRatio*(crys->totalAtoms)*1.0/5.0;
	//if(numHalideVacancies < 1)
	//numHalideVacancies = 1;
	printf("%d total atoms. Replacing %d atoms by vacancies actual concentration %f \n",crys->totalAtoms,numHalideVacancies,(double)numHalideVacancies/(crys->totalAtoms*3.0/5.0));
	crys_replaceRandomAtoms(crys,numHalideVacancies,crys_elementString(1,"I"),1,VX);

	//int numCationVacancies = FAVacancyRatio*(crys->totalAtoms)*1.0/5.0;
	//int numCationVacancies = numHalideVacancies;
	//if(numCationVacancies < 1)
	//numCationVacancies = 1;
	printf("%d total atoms. Replacing %d atoms by vacancies actual concentration %f \n",crys->totalAtoms,numCationVacancies,(double)numCationVacancies/(crys->totalAtoms*1.0/5.0));
	//This command originally would also allow Cs to randomly be replaced by cation vacancies. However, the FieldProfile allocation assumes that every crystal has the same number
	//of dipoles and ions and thus was unable to accomodate having differing numbers of FA in each run.
	crys_replaceRandomAtoms(crys,numCationVacancies,crys_elementString(1,FA),1,VA);

	//This eliminates periodic boundary conditions in the z direction.
	*(crys->latticeVectors+8) *= 2;

    crys_printElements(crys);

    
    return crys;
}

void cvp_loadSettings(char *filename)
{
    printf("Loading settings from file %s\n",filename);
	
	FILE *settingsfile = fopen(filename,"r");
    loadpmcSettings(settingsfile);
    pz_loadSettings(settingsfile);

    
	dir = readString(settingsfile,"dir");    
    readDouble(settingsfile,"FAVacancyRatio",&(FAVacancyRatio));
    readDouble(settingsfile,"IVacancyRatio",&(IVacancyRatio));
    readDouble(settingsfile,"CsRatio",&(CsRatio));
	readDouble(settingsfile,"FAHopEnergy",&(FAHopEnergy));
	readDouble(settingsfile,"CsHopEnergy",&(CsHopEnergy));
	readDouble(settingsfile,"IHopEnergy",&(IHopEnergy));
	readDouble(settingsfile,"defectCharge",&(defectCharge));
	readDouble(settingsfile,"FAPolarizability",&(FAPolarizability));
	readDouble(settingsfile,"IPolarizability",&(IPolarizability));
    readInt(settingsfile,"sizex",&(sizex));
	readInt(settingsfile,"sizey",&(sizey));
	readInt(settingsfile,"sizez",&(sizez));	
	//readInt(settingsfile,"continuationJob",&(continuationJob));	


    fclose(settingsfile);
}

void cvp_saveSettings()
{
    printf("Dir is %s\n",dir);
    mkdir2(dir);
    char *filename = strcat2(dir,"settings");
    printf("Saving settings to %s\n",filename);
    FILE *settingsfile = fopen(filename,"w");

    savepmcSettings(settingsfile);
    pz_saveSettings(settingsfile);
    
	fprintf(settingsfile,"IHopEnergy = %g\n",IHopEnergy);
	fprintf(settingsfile,"FAHopEnergy = %g\n",FAHopEnergy);
	fprintf(settingsfile,"CsHopEnergy = %g\n",CsHopEnergy);
	fprintf(settingsfile,"FAVacancyRatio = %g\n",FAVacancyRatio);
	fprintf(settingsfile,"IVacancyRatio = %g\n",IVacancyRatio);
	fprintf(settingsfile,"CsRatio = %g\n",CsRatio);
	fprintf(settingsfile,"defectCharge = %g\n",defectCharge);
	fprintf(settingsfile,"FAPolarizability = %g\n",FAPolarizability);
	fprintf(settingsfile,"IPolarizability = %g\n",IPolarizability);
	fprintf(settingsfile,"sizex = %d\n",sizex);
	fprintf(settingsfile,"sizey = %d\n",sizey);
	fprintf(settingsfile,"sizez = %d\n",sizez);
	    
    fclose(settingsfile);
	
	char command[1000];
	sprintf(command,"cat %s\n",strcat2(dir,"settings"));
	system(command);
}

void setup()
{
	
    char *ionNames = crys_elementString(2,VX,VA);	
	int numIons = 2;
	double *ionCharge = malloc(2*sizeof(double));
	*ionCharge = defectCharge; *(ionCharge+1) = -defectCharge;
	
	//This should allow the same code to be used regardless of whether there are dipoles or not and reduce code maintainence. 
	char *dipoleNames = NULL;
	int numDipoles = 0;
	double *dipolePolarizability = NULL;
	
	if(FAPolarizability != 0 || IPolarizability != 0)//We could add an extra case where FA is there or I is there. 
	{
	dipoleNames = crys_elementString(2,FA,"I");	
	numDipoles = 2;
	dipolePolarizability = malloc(2*sizeof(double));
	*dipolePolarizability = FAPolarizability;
	*(dipolePolarizability+1) = IPolarizability;
	}
	
	crystal *crys = cvp_crystal();
	double *zvalues = malloc(crys->totalAtoms*sizeof(double));
	int numzvalues;
	crys_zvalues(crys,zvalues,&numzvalues);	
	
	pz_setup(zvalues, numzvalues, ionNames,numIons,ionCharge,dipoleNames,numDipoles,dipolePolarizability,cvp_trajectory);
}

int main(int argc, char **argv)
{
	cvp_loadSettings(argv[1]);
	setup();
    cvp_saveSettings();   
    calculateTrajectories();
		
	return 0;
}



