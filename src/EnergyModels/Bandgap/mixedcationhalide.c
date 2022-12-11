/*
 * pmc.c
 * 
 * Copyright 2017 Anthony Ruth <aruth2@nd.edu>
 * pmc.c simulates the movement of ions in a perovskite lattice
 * considering energy differences. 
 * 
 */

#include "perovskite.h"
#include "bandgap.h"

char dir[1000];
double iodineRatio, FARatio, FAvacancyRatio, IvacancyRatio;
double clHopEnergy, iHopEnergy, faHopEnergy, csHopEnergy;
int	coordinationNumber;//The number of shells used to define the local bandgap.
//This should be resolved.
int shells[] = {12, 60, 168, 408, 708, 1284, 1872, 2928, 3900};//The number of neighbors surrounding a B-site cation in perovskites up to i shells
double latticeConstant = 6.3620;//DOI:10.1021/acs.jpclett.5b01432
int	sizex,sizey,sizez;

double IBrRepulsiveEnergy,FACsRepulsiveEnergy;
double mch_gapEnergy(int numBandgapAlteringElements, int *numEachElement);
void mch_setup();
void mch_registerSettings();
//void mch_saveSettings();
//void mch_loadSettings(char *filename);
crystal *mch_crystal();
void mch_bgSettings();
void mch_trajectory(Trajectory *traj);

double mch_gapEnergy(int numBandgapAlteringElements, int *numEachElement)
{
	//For FAxCs1-xPb(IyCl1-y)3
	//Then apply a correction for A-site and X-site composition
	//CsPbCl3 - 3.04 ev
	//FAPbCl3 - 3.0 eV 
	//CsPbI3 - 1.75 eV
	//FAPbI3 - 1.53 eV
	//CsPbBr3 - 2.43 eV
	//FAPbBr3 - 2.29 eV
	//for x=y=0, the bandgap is 3.04
	//Then there is an x coefficient of 0.2
	//and a y coefficient of 1.3
	
	int numCl,numI,numFA,numCs;
	double x,y;
	numCl = *(numEachElement);
	numI = *(numEachElement+1);
	numFA = *(numEachElement+2);
	numCs = *(numEachElement+3);

	x = (double)numFA/(numFA+numCs);
	y = (double)numI/(numI+numCl);
	
	return 3.04 - 0.2*x - 1.3*y;
}

void mch_trajectory(Trajectory *traj)
{
    /*
     * Provides a description of hopping barriers for a mixed halide crystal
     * The barriers can be direction or element specific.  
     * */
   
	LatticeDynamics *LD = traj->LD = malloc(sizeof(LatticeDynamics));
	bg_trajectory(traj);
   
	LD->numHopPairs = 4;
	LD->hopPairs = crys_elementString(8,VX,"Cl",VX,"I",VA,FA,VA,"Cs");	
	LD->hopPairEnergies = malloc(4*sizeof(double));
	
	*(LD->hopPairEnergies) = clHopEnergy; //Cl
	*(LD->hopPairEnergies+1) = iHopEnergy; //I
	*(LD->hopPairEnergies+2) = faHopEnergy; //FA
	*(LD->hopPairEnergies+3) = csHopEnergy; //Cs
	
    traj->crys = mch_crystal(LD);
    
    
    (traj->nnd)->nndistance=0.75*latticeConstant;
	(traj->nnd)->elementList=crys_elementString(7,"Pb","I","Cl",VA,VX,FA,"Cs");
	(traj->nnd)->numEle=7;
	(traj->nnd+1)->nndistance=1.01*latticeConstant;
	(traj->nnd+1)->elementList=crys_elementString(3,VA,FA,"Cs");
	(traj->nnd+1)->numEle=3;
	traj->numnnds=2;
}


crystal *mch_crystal()
{
    //Setup Crystal
	struct crystal *CsPbCl = perovskite_newCrys("Cs","Pb","Cl","Cl","Cl",latticeConstant);
	struct crystal *crys = crys_multiply(CsPbCl,sizez,sizey,sizez);
	crys_replaceRandomAtoms(crys,FARatio*(sizex*sizey*sizez),"Cs",1,FA);	
	crys_replaceRandomAtoms(crys,3*iodineRatio*(sizex*sizey*sizez),"Cl",1,"I");	
		
	int numHalideVacancies = IvacancyRatio*(crys->totalAtoms)*3.0/5.0;
	if(numHalideVacancies < 2)
	numHalideVacancies = 2;
	printf("%d total atoms. Replacing %d atoms by vacancies actual concentration %f \n",crys->totalAtoms,numHalideVacancies,(double)numHalideVacancies/(crys->totalAtoms*3.0/5.0));
	crys_replaceRandomAtoms(crys,numHalideVacancies,crys_elementString(2,"I","Cl"),2,VX);

	int numCationVacancies = FAvacancyRatio*(crys->totalAtoms)*1.0/5.0;
	if(numCationVacancies < 2)
	numCationVacancies = 2;
	printf("%d total atoms. Replacing %d atoms by vacancies actual concentration %f \n",crys->totalAtoms,numCationVacancies,(double)numCationVacancies/(crys->totalAtoms*3.0/5.0));
	crys_replaceRandomAtoms(crys,numCationVacancies,crys_elementString(2,FA,"Cs"),2,VA);

    //The 0.75 in the cn_fillNetwork function is not satisfactory. 
    //cn_allocateCrystalNetwork(crys);
    //cn_fillNetwork(crys,0.75*latticeConstant,crys_elementString(7,"Pb","I","Cl",VA,VX,FA,"Cs"),7,crys_elementString(7,"Pb","I","Cl",VA,VX,FA,"Cs"),7);
    //cn_fillNetwork(crys,1.01*latticeConstant,crys_elementString(3,VA,FA,"Cs"),3,crys_elementString(3,VA,FA,"Cs"),3);
	
    crys_printElements(crys);
    
    return crys;
}

void mch_registerSettings()
{
	bg_registerSettings();
	registerInt(&coordinationNumber,"coordinationNumber",3);
	registerString(dir,"dir",".");
	registerDouble(&FAvacancyRatio,"FAvacancyRatio",0.00);
	registerDouble(&IvacancyRatio,"IvacancyRatio",0.00);
	registerDouble(&iodineRatio,"iodineRatio",0);
	registerDouble(&(IBrRepulsiveEnergy),"IBrRepulsiveEnergy",0.0);
	registerDouble(&(FACsRepulsiveEnergy),"FACsRepulsiveEnergy",0.0);
	registerDouble(&FARatio,"FARatio",0);
	registerDouble(&faHopEnergy,"faHopEnergy",0.25);
	registerDouble(&csHopEnergy,"csHopEnergy",0.25);
	registerDouble(&clHopEnergy,"clHopEnergy",0.25);
	registerDouble(&iHopEnergy,"iHopEnergy",0.25);
	registerInt(&sizex,"sizex",3);
	registerInt(&sizey,"sizey",3);
	registerInt(&sizez,"sizez",3);
	
}
/*
void mch_loadSettings(char *filename)
{
    printf("Loading settings from file %s\n",filename);
	
	FILE *settingsfile = fopen(filename,"r");
    loadpmcSettings(settingsfile);
    bg_loadSettings(settingsfile);
    
	dir = readString(settingsfile,"dir");
	readInt(settingsfile,"coordinationNumber",&(coordinationNumber));
    readDouble(settingsfile,"vacancyRatio",&(vacancyRatio));
	readDouble(settingsfile,"iodineRatio",&(iodineRatio));
	readDouble(settingsfile,"FARatio",&(FARatio));
	readDouble(settingsfile,"faHopEnergy",&(faHopEnergy));
	readDouble(settingsfile,"csHopEnergy",&(csHopEnergy));
	readDouble(settingsfile,"clHopEnergy",&(clHopEnergy));
	readDouble(settingsfile,"iHopEnergy",&(iHopEnergy));
    readInt(settingsfile,"sizex",&(sizex));
	readInt(settingsfile,"sizey",&(sizey));
	readInt(settingsfile,"sizez",&(sizez));	


    fclose(settingsfile);
}

void mch_saveSettings()
{
    printf("Dir is %s\n",dir);
    mkdir2(dir);
    char *filename = strcat2(dir,"settings");
    printf("Saving settings to %s\n",filename);
    FILE *settingsfile = fopen(filename,"w");

    savepmcSettings(settingsfile);
    bg_saveSettings(settingsfile);
    
	fprintf(settingsfile,"coordinationNumber = %d\n",coordinationNumber);
	fprintf(settingsfile,"clHopEnergy = %g\n",clHopEnergy);
	fprintf(settingsfile,"iHopEnergy = %g\n",iHopEnergy);
	fprintf(settingsfile,"faHopEnergy = %g\n",faHopEnergy);
	fprintf(settingsfile,"csHopEnergy = %g\n",csHopEnergy);
	fprintf(settingsfile,"iodineRatio = %g\n",iodineRatio);
	fprintf(settingsfile,"FARatio = %g\n",FARatio);
	fprintf(settingsfile,"vacancyRatio = %g\n",vacancyRatio);
	fprintf(settingsfile,"sizex = %d\n",sizex);
	fprintf(settingsfile,"sizey = %d\n",sizey);
	fprintf(settingsfile,"sizez = %d\n",sizez);
    
    fclose(settingsfile);
	
	char command[1000];
	sprintf(command,"cat %s\n",strcat2(dir,"settings"));
	system(command);
}*/

//A lot of this belongs in bandgap.c
void mch_setup()
{
	char *coordElement = "Pb";
    char *bandgapAlteringElements = crys_elementString(4,"Cl","I",FA,"Cs");	
	int numBandgapAlteringElements = 4;
	int numStates = sizex*sizey*sizez;
	int repulsiveShellSize = 8;
	int maxRepulsiveCoordination = 10;
	char *repulsiveElements = crys_elementString(4,"Br","I",FA,"Cs");
	double *repulsiveEnergies = calloc(2,sizeof(double *));
	*(repulsiveEnergies) = IBrRepulsiveEnergy;
	*(repulsiveEnergies+1) = FACsRepulsiveEnergy;
	int numRepulsiveElements = 4;
	
	bg_setup(mch_gapEnergy,coordElement,bandgapAlteringElements,numBandgapAlteringElements,shells,numStates,repulsiveElements,numRepulsiveElements,repulsiveEnergies,repulsiveShellSize,maxRepulsiveCoordination,mch_trajectory);
	//pmcSetup(mixedCationHalideLatticeDynamics,OS_allocate,OS_free,saveOptoelectronicState,sizeof(OptoelectronicState));
	
}

int main(int argc, char **argv)
{
	allocateSettings();
	mch_registerSettings();
	
	FILE *infile = fopen(argv[1],"r");
	loadSettings(infile);
	fclose(infile);
	mkdir2(dir);
	 
	char outfileName[1000];
	sprintf(outfileName,"%s/settings",dir);
	
	FILE *outfile = fopen(outfileName,"w");
	saveSettings(outfile);
	fclose(outfile);
		    
    mch_setup();   
    simulateTrajectories();
    
	return 0;
}

