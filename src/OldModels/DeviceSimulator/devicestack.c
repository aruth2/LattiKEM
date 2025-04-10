/*
 * chargedcationhalidevacancies.c
 * 
 * Copyright 2019 Anthony Ruth <aruth@heesolar.com>
 * Simulates the movement of FA and I vacancies in a perovskite lattice
 * considering ionic drift and diffusion. 
 * 
 */

#include "perovskite.h"
#include "stacd.h"
//For now the other layers besides the perovskite will not be shown. I will work on drawing them in the future. 
char dir[1000];
double FAVacancyRatio,IVacancyRatio;
double IHopEnergy, FAHopEnergy;
double FAPbI3latticeConstant = 6.3620;//DOI:10.1021/acs.jpclett.5b01432
int	sizex,sizey,sizez;

void ds_trajectory(Trajectory *traj);
crystal *ds_crystal();
void ds_registerSettings();
void setup();

void ds_trajectory(Trajectory *traj)
{
     traj->crys = ds_crystal(traj);
     
     stacd_Trajectory(traj);
     
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
	(traj->nnd)->nndistance=0.75*FAPbI3latticeConstant;
	(traj->nnd)->elementList=crys_elementString(5,"Pb","I",VA,VX,FA);
	(traj->nnd)->numEle=5;
	(traj->nnd+1)->nndistance=1.01*FAPbI3latticeConstant;
	(traj->nnd+1)->elementList=crys_elementString(2,VA,FA);
	(traj->nnd+1)->numEle=2;
	traj->numnnds=2;
}

crystal *ds_crystal()
{
	struct crystal *perovcrys,*crys;

    //Setup Crystal
	struct crystal *FAPbI = perovskite_newCrys(FA,"Pb","I","I","I",FAPbI3latticeConstant);
	crys = crys_multiply(FAPbI,sizex,sizey,sizez);
	
	int numHalideVacancies = IVacancyRatio*(crys->totalAtoms)*3.0/5.0;
	int numCationVacancies = FAVacancyRatio*(crys->totalAtoms)*1.0/5.0;
	
	//printf("%d total atoms. Replacing %d atoms by vacancies actual concentration %f \n",crys->totalAtoms,numHalideVacancies,(double)numHalideVacancies/(crys->totalAtoms*3.0/5.0));
	crys_replaceRandomAtoms(crys,numHalideVacancies,crys_elementString(1,"I"),1,VX);

	//printf("%d total atoms. Replacing %d atoms by vacancies actual concentration %f \n",crys->totalAtoms,numCationVacancies,(double)numCationVacancies/(crys->totalAtoms*1.0/5.0));
	crys_replaceRandomAtoms(crys,numCationVacancies,crys_elementString(1,FA),1,VA);

	//This eliminates periodic boundary conditions in the z direction.
	*(crys->latticeVectors+8) *= 2;
	
	int numNiOLayers = 2 * 400 / FAPbI3latticeConstant;//Assumes 40 nm
	int numC60Layers = 2 * 400 / FAPbI3latticeConstant;//Assumes 40 nm
	
	return crys;
}
crystal *ds_crystalPerovskiteOnly()
{
	struct crystal *perovcrys,*crys;

    //Setup Crystal
	struct crystal *FAPbI = perovskite_newCrys(FA,"Pb","I","I","I",FAPbI3latticeConstant);
	crys = crys_multiply(FAPbI,sizex,sizey,sizez);
	
	int numHalideVacancies = IVacancyRatio*(crys->totalAtoms)*3.0/5.0;
	int numCationVacancies = FAVacancyRatio*(crys->totalAtoms)*1.0/5.0;
	//I think we will let the charge balance on its own for now. It will help confirm the validity of the methods use
	
	printf("%d total atoms. Replacing %d atoms by vacancies actual concentration %f \n",crys->totalAtoms,numHalideVacancies,(double)numHalideVacancies/(crys->totalAtoms*3.0/5.0));
	crys_replaceRandomAtoms(crys,numHalideVacancies,crys_elementString(1,"I"),1,VX);

	printf("%d total atoms. Replacing %d atoms by vacancies actual concentration %f \n",crys->totalAtoms,numCationVacancies,(double)numCationVacancies/(crys->totalAtoms*1.0/5.0));
	crys_replaceRandomAtoms(crys,numCationVacancies,crys_elementString(1,FA),1,VA);

	//This eliminates periodic boundary conditions in the z direction.
	*(crys->latticeVectors+8) *= 2;

    crys_printElements(crys);
    
    return crys;
}

void ds_registerSettings()
{
	registerString(dir,"dir",".");
	registerDouble(&FAVacancyRatio,"FAVacancyRatio",0.0);
	registerDouble(&IVacancyRatio,"IVacancyRatio",0.0);
	registerDouble(&FAHopEnergy,"FAHopEnergy",0.25);
	registerDouble(&IHopEnergy,"IHopEnergy",0.25);
	registerInt(&sizex,"sizex",4);
	registerInt(&sizey,"sizey",4);
	registerInt(&sizez,"sizez",4);
	stacd_registerSettings();
}
void setup()
{
	
    char *ionNames; 	
	int numIonTypes,numIons;
	int numMaterials=5;
	numIonTypes = 2;
	int totalDefectLevels = 2;
	
	double *ionCharge = malloc(numIonTypes*sizeof(double));
	int *materialTransitions = calloc(numMaterials,sizeof(int));
	char **materialNames = calloc(numMaterials,sizeof(char*));
	double *fermiLevels = calloc(numMaterials,sizeof(double));
	double *bandOffsets = calloc(2*numMaterials,sizeof(double));
	int *numIonDefectLevels = calloc(numIonTypes,sizeof(double));
	double *defectLevels = calloc(totalDefectLevels,sizeof(double));
	double *defectLevelCharge = calloc(totalDefectLevels,sizeof(double));
	
	double *dielectricConstants = calloc(numMaterials,sizeof(double));
	int *absorbers=calloc(numMaterials,sizeof(int));
	
	DeviceStack *architecture = malloc(sizeof(DeviceStack));
	
	//ionNames = crys_elementString(4,VX,VA,"D","A");
	ionNames = crys_elementString(2,VX,VA);
	*ionCharge = 1; *(ionCharge+1) = -1;
	//*ionCharge = 0; *(ionCharge+1) = 0;
	//*(ionCharge+2) = 1; *(ionCharge+3) = -1;
	//*(ionCharge+2) = 0; *(ionCharge+3) = 0;

	crystal *crys = ds_crystalPerovskiteOnly();

	//This allocation is not very precise
	double *zValues = malloc(crys->totalAtoms*sizeof(double));
	double *formulaUnits = calloc(crys->totalAtoms,sizeof(double));
	double *thickness = calloc(crys->totalAtoms,sizeof(double));
	int numzValues;
	crys_zValues(crys,zValues,&numzValues);	
	
	*(materialNames+0)="ITO";
	*(materialNames+1)="NiO";
	*(materialNames+2)="FAPbI3";
	*(materialNames+3)="C60";
	*(materialNames+4)="Cu";
	
	int numTiO2Layers = 6;
	int numNiOLayers = 2 * 400 / FAPbI3latticeConstant;//Assumes 40 nm
	int numPerovskiteLayers = 2 * sizez;
	int numC60Layers = 2 * 400 / FAPbI3latticeConstant;//Assumes 40 nm
	int numCuLayers = 6;
	int iLayer;
	for(iLayer=0;iLayer<numTiO2Layers;iLayer++)//2 layers are added already by padding in crys_zValues
	padListLeft(zValues,&numzValues);
	for(iLayer=0;iLayer<numNiOLayers;iLayer++)
	padListLeft(zValues,&numzValues);
	for(iLayer=0;iLayer<numC60Layers;iLayer++)
	padListRight(zValues,&numzValues);
	for(iLayer=0;iLayer<numCuLayers;iLayer++)//2 layers are added already by padding in crys_zValues
	padListRight(zValues,&numzValues);
	
	//This should not be necessary needs cleanup.
	crys = ds_crystal();
	
	
	int TiO2start = 0;
	int NiOstart = TiO2start+numTiO2Layers;
	int Perovskitestart = NiOstart+numNiOLayers;
	int C60start = Perovskitestart+numPerovskiteLayers;
	int Custart = C60start + numC60Layers;
	
	*(materialTransitions+0)=TiO2start;
	*(materialTransitions+1)=NiOstart;
	*(materialTransitions+2)=Perovskitestart;
	*(materialTransitions+3)=C60start;
	*(materialTransitions+4)=Custart;
	
	double *defectDensities = calloc(numMaterials,sizeof(double));
	double *bandDefectLevels = calloc(numMaterials,sizeof(double));
	double *bandDefectCharge = calloc(numMaterials,sizeof(double));
	
	//This should be corrected later to include differences in DOS
	for(iLayer=0;iLayer<numzValues;iLayer++)
	{
		*(formulaUnits+iLayer) = sizex*sizey/2.0;
		*(thickness + iLayer) = FAPbI3latticeConstant/2;
	}
	
	//*(bandOffsets+0*2) = -4.1; *(bandOffsets+0*2+1) = -4.1;
	//*(bandOffsets+1*2) = -1.8; *(bandOffsets+1*2+1) = -5.4;
	//*(bandOffsets+1*2) = -1.8; *(bandOffsets+1*2+1) = -5.4;
	//*(bandOffsets+2*2) = -3.9; *(bandOffsets+2*2+1) = -5.5;
	//*(bandOffsets+3*2) = -4.0; *(bandOffsets+3*2+1) = -7.6; //C60
	//*(bandOffsets+4*2) = -5.3; *(bandOffsets+4*2+1) = -5.3;//For Cu	

	*(bandOffsets+0*2) = -5.0; *(bandOffsets+0*2+1) = -5.0;
	//*(bandOffsets+1*2) = -1.8; *(bandOffsets+1*2+1) = -5.4;
	*(bandOffsets+1*2) = -1.55; *(bandOffsets+1*2+1) = -5.15;
	*(bandOffsets+2*2) = -4.07; *(bandOffsets+2*2+1) = -5.6;
	*(bandOffsets+3*2) = -4.1; *(bandOffsets+3*2+1) = -6.4; //C60
	*(bandOffsets+4*2) = -4.2; *(bandOffsets+4*2+1) = -4.2;//For Cu	
	
	/*
	*(bandOffsets+0*2) = -5.6; *(bandOffsets+0*2+1) = -5.6;
	*(bandOffsets+1*2) = -2.6; *(bandOffsets+1*2+1) = -5.6;
	*(bandOffsets+2*2) = -3.9; *(bandOffsets+2*2+1) = -5.6;
	*(bandOffsets+3*2) = -3.9; *(bandOffsets+3*2+1) = -6.9; //C60
	*(bandOffsets+4*2) = -3.9; *(bandOffsets+4*2+1) = -3.9;//For Cu	
	*/
	/*
	*(bandOffsets+0*2) = -4.8; *(bandOffsets+0*2+1) = -4.8;
	*(bandOffsets+1*2) = -1.8; *(bandOffsets+1*2+1) = -5.4;
	*(bandOffsets+2*2) = -4; *(bandOffsets+2*2+1) = -5.6;
	*(bandOffsets+3*2) = -4.2; *(bandOffsets+3*2+1) = -7.8; //BCP
	*(bandOffsets+4*2) = -4.8; *(bandOffsets+4*2+1) = -4.8;//For Cu	
	*/

	*(numIonDefectLevels+0) = 1;//VI Electron Trap
	*(numIonDefectLevels+1) = 1;//VFa Hole Trap
	//*(numIonDefectLevels+2) = 1;//ETL Electron Trap
	//*(numIonDefectLevels+3) = 1;//HTL Hole Trap
	
	*(defectLevels+0) = -0.2; *(defectLevelCharge+0) = -1;//VI Electron Trap
	*(defectLevels+1) = 0.2; *(defectLevelCharge+1) = 1;//VFa Hole Trap
	//*(defectLevels+2) = -0.4; *(defectLevelCharge+2) = -1;//ETL Electron Trap //Is this still used?
	//*(defectLevels+3) = 0.4; *(defectLevelCharge+3) = 1;//HTL Hole Trap
	
	*(defectDensities+0) = 0; *(bandDefectCharge+0) = 0.0; *(bandDefectLevels+0) = 0;//ITO
	//*(defectDensities+1) = 1.39e-5; *(bandDefectCharge+1) = -1; *(bandDefectLevels+1) = 0.3;//NiO
	*(defectDensities+1) = 1e-6; *(bandDefectCharge+1) = -1; *(bandDefectLevels+1) = 0.3;//NiO
	//*(defectDensities+1) = 0; *(bandDefectCharge+1) = -1; *(bandDefectLevels+1) = 0.4;//NiO
	//*(defectDensities+1) = 6e-6; *(bandDefectCharge+1) = -1; *(bandDefectLevels+1) = 0.4;//NiO
	*(defectDensities+2) = 0; *(bandDefectCharge+2) = 0.0; *(bandDefectLevels+2) = 0;//Perovskite
	*(defectDensities+3) = 1e-8; *(bandDefectCharge+3) = 1; *(bandDefectLevels+3) = -0.3;//BCP
	//*(defectDensities+3) = 0; *(bandDefectCharge+3) = 1; *(bandDefectLevels+3) = -0.4;//BCP
	//*(defectDensities+3) = 6e-6; *(bandDefectCharge+3) = 1; *(bandDefectLevels+3) = -0.4;//BCP
	*(defectDensities+4) = 0; *(bandDefectCharge+4) = 0.0; *(bandDefectLevels+4) = 0;//Cu
	
	architecture->defectDensities = defectDensities;
	architecture->bandDefectCharge = bandDefectCharge;
	architecture->bandDefectLevels = bandDefectLevels;
	
	architecture->numDefectLevels = 0;
	architecture->numDefectStates = 0;
	
	int iEle;
	for(iEle=0;iEle<numIonTypes;iEle++)
	{
	architecture->numDefectLevels += *(numIonDefectLevels+iEle);
	architecture->numDefectStates += *(numIonDefectLevels+iEle) * crys_elementCount(crys,ionNames+namelength*iEle);
	}
	
	architecture->numStates = 2*numzValues + architecture->numDefectStates;
	
	printf("The number of states is %d\n",architecture->numStates);
	
	*(dielectricConstants+0) = 10000;
	//*(dielectricConstants+1) = 9.1;
	*(dielectricConstants+1) = 10.7;
	*(dielectricConstants+2) = 10;
	//*(dielectricConstants+2) = 18;
	//*(dielectricConstants+3) = 4;
	*(dielectricConstants+3) = 4.4;
	//*(dielectricConstants+3) = 9.1;
	*(dielectricConstants+4) = 10000;
	
	*(absorbers+0) = 0;
	*(absorbers+1) = 0;
	*(absorbers+2) = 1;
	*(absorbers+3) = 0;
	*(absorbers+4) = 0;
	//int numAbsorberLayers= 2*sizez; //This could be figured out from the absorbers
	
	for(iEle=0,numIons=0;iEle<numIonTypes;iEle++)
	numIons += crys_elementCount(crys,ionNames+iEle*namelength);
	
	architecture->ionAtoms = calloc(numIons,sizeof(int));
	architecture->ionEles = calloc(numIons,sizeof(int));
	
	architecture->numIonTypes=numIonTypes;
	architecture->numIons=numIons;
	
	architecture->ionCharge=ionCharge;
	architecture->ionNames=ionNames;
	
	architecture->materialTransitions=materialTransitions;
	architecture->materialNames=materialNames;
	architecture->bandOffsets=bandOffsets;
	architecture->defectLevels=defectLevels;
	architecture->dielectricConstants=dielectricConstants;
	architecture->formulaUnits=formulaUnits;
	architecture->thickness=thickness;
	architecture->zValues=zValues;
	architecture->numzValues=numzValues;
	
	architecture->defectLevelCharge = defectLevelCharge;
	architecture->numIonDefectLevels = numIonDefectLevels;
	architecture->numMaterials=numMaterials;
	architecture->absorbers=absorbers;
	
	architecture->volumePerUnitCell = pow(FAPbI3latticeConstant,3);//DOI:10.1021/acs.jpclett.5b01432
	
	architecture->initializationCrystal = crys;
	fillIonAtoms(architecture);
	compressZProfile(architecture);
	printDeviceStack(architecture);
	
	stacd_setup(architecture,ds_trajectory);
}

int main(int argc, char **argv)
{
	allocateSettings();
	ds_registerSettings();
	
	
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



