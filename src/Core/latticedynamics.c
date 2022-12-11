#include "latticedynamics.h"
double (*energyFunction)(Configuration *);
void (*savingFunction)(Configuration *);
void (*breakpoint_allocator)(void *);
void (*breakpoint_freer)(void *);
//This should be removed in favor of autocreating the weights and calling weightedAveragingFunction
void (*averagingFunction)(Configuration *,int, Configuration *);
void (*weightedAveragingFunction)(Configuration *,int, Configuration *, double *);
size_t dataSize;
void (*traj_generator)(Trajectory *);
//Trajectory *settingsTrajectory;
int breakpointmethod;
int startAtStep=-1;

int numSteps;
int hopMode;
double temperature;

void LD_listMoves(crystal *crys, LatticeDynamics *LD, MarkovChainMonteCarlo *mcmc)
{
	//printf("Enumerating moves\n");
	mcmc->numMoves = 0;
	
	int moveType;
	int species1start,species1end;
	int iatom1,ineighbor,iatom2;
	
	for(moveType = 0;moveType<LD->numHopPairs;moveType++)
	{
		species1start = crys_elementOffset(crys,LD->hopPairs+2*moveType*namelength);
		species1end = species1start+crys_elementCount(crys,LD->hopPairs+2*moveType*namelength);
		//printf("species %s goes from %d to %d\n",LD->hopPairs+2*moveType*namelength,species1start,species1end);
		for(iatom1=species1start;iatom1<species1end;iatom1++)
		{
			for(ineighbor=0;ineighbor<*(crys->network->numAdjacent+iatom1);ineighbor++)
			{
				iatom2=*(crys->network->adjacencyList+iatom1*maxconnections+ineighbor);
				
				if(!strcmp(crys->species+iatom2*namelength,LD->hopPairs+(2*moveType+1)*namelength))
				{
					//printf("Atom %d %s and atom %d %s is a valid move\n",iatom1,crys->species+iatom1*namelength,iatom2,crys->species+iatom2*namelength);
					*(mcmc->moves+2*mcmc->numMoves) = iatom1;
					*(mcmc->moves+2*mcmc->numMoves+1) = iatom2;
					*(mcmc->moveBarriers+mcmc->numMoves)=*(LD->hopPairEnergies+moveType);
					(mcmc->numMoves)++;
				}
			}
		}
		
	}
	
}

void mcmc_printMoves(crystal *crys, MarkovChainMonteCarlo *mcmc)
{
	int iMove;
	for(iMove = 0;iMove<mcmc->numMoves;iMove++)
	{
		printf("Move %d is swapping %d %s and %d %s with a barrier of %g\n",iMove,*(mcmc->moves+2*iMove),crys->species+*(mcmc->moves+2*iMove)*namelength,*(mcmc->moves+2*iMove+1),crys->species+*(mcmc->moves+2*iMove+1)*namelength,*(mcmc->moveBarriers+iMove));
	}
}

void mcmc_allocate(MarkovChainMonteCarlo *mcmc, int maxMoves)
{
	mcmc->moves= malloc(2*maxMoves*sizeof(int));
	mcmc->moveBarriers= malloc(maxMoves*sizeof(double));
    mcmc->moveEnthalpies = malloc(maxMoves*sizeof(double));
    mcmc->moveRates = malloc(maxMoves*sizeof(double));
    mcmc->moveProbabilityRanges = malloc(maxMoves*sizeof(double));
    mcmc->initialenergy=1e308;
}
/*
void traj_logarithmicBreakpoints(Trajectory *traj)
{
    int step, stepsize, enthalpyBreakpoint;
    for(enthalpyBreakpoint=0;enthalpyBreakpoint<traj->numEnthalpyBreakpoints;enthalpyBreakpoint++)
    {
		for(step=1+*(traj->enthalpyStateBreakpoints+enthalpyBreakpoint),stepsize=1;step<*(traj->enthalpyStateBreakpoints+enthalpyBreakpoint+1);step+=stepsize,stepsize*=2)
		{
			*(traj->breakpoints+traj->numbreakpoints) = step;
			(traj->numbreakpoints)++;
		}
	}
}*/

void traj_allBreakpoints(Trajectory *traj)
{
    int step;
    int numSteps = getNumSteps();
    for(step=0;step<numSteps;step++)
    {
        *(traj->breakpoints+step) = step;
    }
    (traj->numbreakpoints) = numSteps;
}
/*
void traj_biasBreakpoints(Trajectory *traj)
{
    int step;
    traj->numbreakpoints=0;
    int stategetter;
    for(step=0;step<traj->numSteps;step++)
    {
		stategetter = binarysearch(traj->enthalpyStateBreakpoints,step,traj->numEnthalpyBreakpoints);
		if(stategetter <0)
		continue;
				 
		if(*(traj->enthalpyStates+stategetter-1))
		{
			
        *(traj->breakpoints+traj->numbreakpoints) = step;
        (traj->numbreakpoints)++;
		}
    }
}*/

//Update the crystal stored on the thread until it reaches the described step
void traj_crysAtStep(crystal *crys, int *usedMoves, int initialStep, int targetStep, int hasNetwork)
{
	//We could add an automatic performance boost here. 
	//Either every swap is performed until the structure is updated, O(maxconnections^2*nSteps) or
	//Perform every swap without updating the network and then rebuild the network at end O(nAtoms^2)
	int iStep;
	//Move forward until reaching the correct step
	//printf("Moving crystal from step %d to step %d\n",initialStep,targetStep);
    for(iStep=initialStep;iStep<targetStep;iStep++)
    {
		cn_swap(crys,*(usedMoves+2 * iStep),*(usedMoves+2 * iStep+1),hasNetwork);
	}
    
    //Move backwards until reaching the correct step        
    for(iStep=initialStep;iStep>targetStep;iStep--)
    {
		cn_swap(crys,*(usedMoves+2 * (iStep-1)),*(usedMoves+2 * (iStep-1)+1),hasNetwork);
	}    
            
}

void traj_registerSettings()
{
	//settingsTrajectory = malloc(sizeof(Trajectory));
	//registerInt(&(settingsTrajectory->numSteps),"numSteps",1000);
	registerInt(&numSteps,"numSteps",1000);
	registerInt(&startAtStep,"startAtStep",-1);	
		
	registerDouble(&temperature,"temperature",0.025);
	
	registerEnum(2,&hopMode,"hopMode",HOP_MODE_KINETIC,"metropolis","kinetic");
	registerEnum(3,&breakpointmethod,"breakpoints",ALL_BREAKPOINTS,"all","log","bias");
}

int getNumSteps()
{
	return numSteps;
}
int getHopMode()
{
	return hopMode;
}
double getTemperature()
{
	return temperature;
}

void LD_setup(double (*newEnergyFunction)(Configuration *), void (*newSavingFunction)(Configuration *), void (*newbreakpoint_allocator)(void *), void (*newbreakpoint_freer)(void *), void (*newAveragingFunction)(Configuration *,int, Configuration *), void (*newWeightedAveragingFunction)(Configuration *,int, Configuration *, double *), void (*newtraj_generator)(Trajectory *), int newDataSize)
{
	energyFunction=newEnergyFunction;
	savingFunction=newSavingFunction;
	dataSize = newDataSize;
	printf("Size of a single configuration is %d\n",dataSize);
	breakpoint_allocator=newbreakpoint_allocator;
	breakpoint_freer = newbreakpoint_freer;
	averagingFunction=newAveragingFunction;
	weightedAveragingFunction = newWeightedAveragingFunction;
	
	traj_generator=newtraj_generator;

	//traj_generator(settingsTrajectory);
	
}

void config_allocate(Configuration *config, Trajectory *traj)
{
	//This assumes the configuration itself has been allocated by malloc(sizeof(Configuration)).
	//The configuration is typically tied to a trajectory.
	config->data=malloc(dataSize);
	//This is causing a memory leak. This pointer is dropped in kineticThread and replaced with a pointer to the trajectory. 
	//However, during config_weightedAverage, it is expected that the pointer is there for the target config. 
	config->seriesData=malloc(traj->numSeriesData*sizeof(double));
	breakpoint_allocator(config->data);
	
}

void config_free(Configuration *config)
{
	//This assumes the configuration itself will also be freed
	breakpoint_freer(config->data);
	free(config->data);
	free(config->seriesData);
}

void config_save(Configuration *config, Trajectory *traj, int step, int maxSteps)
{
	//if(!config_toBeSaved(config, traj, step, maxSteps))
	//return;
	
    char *stepString = formatStep(step, maxSteps);
    sprintf(config->dataFileName,"%s/traj/state%s",traj->dir,stepString);
	//printf("The save flag is %d\n",config->savedata);
	//This should no longer be necessary. If we are calling this function, we are saving. Ditto for the saveData flag in input files
	//if(config->savedata)
		savingFunction(config);
	//if(config->savecrys)
	//	crys_makexyz(config->crys,name);
	free(stepString);
}

void config_savecrys(Configuration *config, Trajectory *traj, int step, int maxSteps)
{
	char *stepString = formatStep(step, maxSteps);
	sprintf(config->crysFileName,"%s/traj/%s.xyz",traj->dir,stepString); 
	crys_makexyz(config->crys,config->crysFileName);
	free(stepString);
}
double config_energy(Configuration *config)
{
	//This function acts as a wrapper so that all files which include latticedynamics.h can call
	//the same energy function which is only set once.
	return energyFunction(config);
}

void config_average(Configuration * configs,int numConfigs, Configuration *outConfig)
{
	//This function acts as a wrapper so that all files which include latticedynamics.h can call
	//the same averaging function which is only set once.
	double *unitArray = unit(numConfigs);
	
	//STACD does not have an averaging function only a weighted averaging function. We are going to remove the general averaging function
	//averagingFunction(configs,numConfigs,outConfig);
	weightedAveragingFunction(configs,numConfigs,outConfig,unitArray);
	free(unitArray);
}

void config_weightedAverage(Configuration * configs,int numConfigs, Configuration *outConfig, double *weights)
{
	//This function acts as a wrapper so that all files which include latticedynamics.h can call
	//the same averaging function which is only set once.
	weightedAveragingFunction(configs,numConfigs,outConfig, weights);
	
	//Average the energy. I could not find any other place that did this, and I am surprised. It's not clear the energy needs to be averaged for config. It looks like the energy gets averaged for trajectories.
	/*iConfig;
	double sumWeight = 0;
	outConfig->energy = 0;
	for(iConfig = 0; iConfig < numConfigs; iConfig++)
		{
		outConfig->energy += *(weights + iConfig) * (configs+iConfig)->energy;	
		sumWeight += *(weights + iConfig);
	}
	outConfig->energy /= sumWeight;
	*/	
}

//This is no longer used since all saving is done by post-processing. Saving is achieved using the data swap which should already have the save flag set
/*
int config_toBeSaved(Configuration *config, Trajectory *traj, int step, int maxSteps)
{
	//Returns whether or not the config will get saved. This is useful to reduced unnecessary calculations for information that will not get saved.
	//This will also set the saveData and savecrys flags
	if(binarysearch(traj->breakpoints,step,traj->numbreakpoints)<0)
	{
	//config->savecrys = traj->savecrys;
	config->savedata = traj->savedata;
	return 1;
	}
	else
	{
	//config->savecrys = 0;
	config->savedata = 0;
	return 0;
	}
}*/

void traj_average(Trajectory *trajectories, int numRuns, Trajectory *outTraj)
{
	//Combine series data, energy, and time series
	int iStep,iData,iTraj;
	int numSteps = getNumSteps();
	for(iStep=0;iStep<numSteps;iStep++)
	{
		for(iData=0;iData<outTraj->numSeriesData;iData++)
		{
			*(outTraj->seriesData+iStep*(outTraj->numSeriesData)+iData)=0;
				for(iTraj=0;iTraj<numRuns;iTraj++)
				*(outTraj->seriesData+iStep*outTraj->numSeriesData+iData)+=*((trajectories+iTraj)->seriesData+iStep*((trajectories+iTraj)->numSeriesData)+iData);
				
				*(outTraj->seriesData+iStep*outTraj->numSeriesData+iData)/=numRuns;

		}
		*(outTraj->energies+iStep)=0;
		*(outTraj->timeseries+iStep)=0;
		*(outTraj->timeSteps+iStep)=0;
		for(iTraj=0;iTraj<numRuns;iTraj++)
		{
			*(outTraj->energies+iStep)+=*((trajectories+iTraj)->energies+iStep);
			//We could align the times series many different ways. 
			*(outTraj->timeseries+iStep)+=*((trajectories+iTraj)->timeseries+iStep);
			*(outTraj->timeSteps+iStep)+=*((trajectories+iTraj)->timeSteps+iStep);
		}
		*(outTraj->energies+iStep)/=numRuns;
		*(outTraj->timeseries+iStep)/=numRuns;
		*(outTraj->timeSteps+iStep)/=numRuns;
	}
	
}

void traj_weightedAverage(Trajectory *trajectories, int numRuns, Trajectory *outTraj, int longestRun )
{
	//Combine series data, energy, and time series
	int longestRunStep,thisRunStep,iData,iTraj;
	//int numSteps = getNumSteps();
	double weights[2*numRuns];
	double longestRunTime;
	printf("There are %d seriesData to average\n",outTraj->numSeriesData);
		
	for(longestRunStep=0;longestRunStep<numSteps;longestRunStep++)
	{
		longestRunTime = *((trajectories+longestRun)->timeseries+longestRunStep);
		//for(iTraj=0;iTraj<numRuns;iTraj++)
		//	linearInterpolationWeights((trajectories+iTraj)->timeseries,longestRunTime,numSteps,&thisRunStep,weights+2*iTraj,weights+2*iTraj+1);
		
		//if(thisRunStep == numSteps-1)
		//thisRunStep--;		
		for(iData=0;iData<outTraj->numSeriesData;iData++)
		{
			
			*(outTraj->seriesData+longestRunStep*(outTraj->numSeriesData)+iData)=0;
			for(iTraj=0;iTraj<numRuns;iTraj++)
				{
				linearInterpolationWeights((trajectories+iTraj)->timeseries,longestRunTime,numSteps,&thisRunStep,weights+2*iTraj,weights+2*iTraj+1);
		
				if(thisRunStep == numSteps-1)
					thisRunStep--;
				
				//printf("Using series data %g and %g from steps %d and %d of run %d data %d\n",*((trajectories+iTraj)->seriesData+thisRunStep*((trajectories+iTraj)->numSeriesData)+iData),*((trajectories+iTraj)->seriesData+(thisRunStep+1)*((trajectories+iTraj)->numSeriesData)+iData),thisRunStep,thisRunStep+1,iTraj,iData);
				*(outTraj->seriesData+longestRunStep*outTraj->numSeriesData+iData)+=*(weights+2*iTraj)**((trajectories+iTraj)->seriesData+thisRunStep*((trajectories+iTraj)->numSeriesData)+iData);
				*(outTraj->seriesData+longestRunStep*outTraj->numSeriesData+iData)+=*(weights+2*iTraj+1)**((trajectories+iTraj)->seriesData+(thisRunStep+1)*((trajectories+iTraj)->numSeriesData)+iData);
				}
				
			*(outTraj->seriesData+longestRunStep*outTraj->numSeriesData+iData)/=numRuns;

		}
		*(outTraj->energies+longestRunStep)=0;
		*(outTraj->timeseries+longestRunStep)=0;
		*(outTraj->timeSteps+longestRunStep)=0;
		for(iTraj=0;iTraj<numRuns;iTraj++)
		{
			linearInterpolationWeights((trajectories+iTraj)->timeseries,longestRunTime,numSteps,&thisRunStep,weights+2*iTraj,weights+2*iTraj+1);
			
			if(thisRunStep == numSteps-1)
			thisRunStep--;	
			
			*(outTraj->energies+longestRunStep)+=*(weights+2*iTraj)**((trajectories+iTraj)->energies+thisRunStep);
			*(outTraj->energies+longestRunStep)+=*(weights+2*iTraj+1)**((trajectories+iTraj)->energies+thisRunStep+1);
			//We could align the times series many different ways. 
			//*(outTraj->timeseries+longestRunStep)+=*(weights+2*iTraj)**((trajectories+iTraj)->timeseries+thisRunStep);
			//*(outTraj->timeseries+longestRunStep)+=*(weights+2*iTraj+1)**((trajectories+iTraj)->timeseries+thisRunStep+1);
		}
		*(outTraj->energies+longestRunStep)/=numRuns;
		//The entire timeseries from the longestRun is what gets used. This could be changed, but it would need to be performed earlier
		*(outTraj->timeseries+longestRunStep)=*((trajectories+longestRun)->timeseries+longestRunStep);
		*(outTraj->timeSteps+longestRunStep)=*((trajectories+longestRun)->timeSteps+longestRunStep);
	}
	
}


void traj_loadSeries(Trajectory *traj)
{
	traj_loadEnergy(traj);
	traj_loadMoves(traj);
	traj_loadSeriesData(traj);
}

void traj_saveSeries(Trajectory *traj)
{
	traj_saveEnergy(traj);
	traj_saveMoves(traj);
	traj_saveSeriesData(traj);
}

//The number of steps cannot exceed this
#define maxfiles 10000000
int traj_load(Trajectory *traj)
{
	traj_loadSeries(traj);
	
	if(traj->step==0)
	{
		printf("No previous record found for run %d\n",traj->iTraj);
		return traj->state = TRAJECTORY_STATE_RUNNING;  
	}
	
	//Starting from the final step appeared more efficient. However this creates a problem if the number of steps is reduced, and the final structure does not actually correspond to the highest step
	//Using only the initial should not be too costly.
	//traj->crys=crys_readxyz(strcat2(traj->dir,"final.xyz"));
	
	//We could load saved .xyz files instead of only allowing first/last to be read.
	//if(traj->crys == NULL)
	//{
		traj->crys = crys_readxyz(strcat2(traj->dir,"initial.xyz"));
		traj_crysAtStep(traj->crys,traj->selectedMoves,0,traj->step,CN_NO_NETWORK);
		printf("Crystal %s loaded\n",strcat2(traj->dir,"initial.xyz"));
	//}
	//else
		//printf("Crystal %s loaded\n",strcat2(traj->dir,"final.xyz"));
	
	if(traj->crys == NULL)
		printf("No crystal found, unsure what to do.\n");
	
	if(startAtStep != -1)
	{
		printf("It has been specified to start at step %d\n",startAtStep);
		traj_crysAtStep(traj->crys,traj->selectedMoves,traj->step,startAtStep,CN_NO_NETWORK);
		traj->step=startAtStep;
	}
	
	if(traj->step<numSteps)
	{
		printf("Trajectory is on step %d of %d. It will be continued to finish the job\n",traj->step,numSteps);
		return traj->state = TRAJECTORY_STATE_RUNNING;
	}
	else
	{
		printf("Trajectory %d has finished\n",traj->iTraj);
		return traj->state = TRAJECTORY_STATE_FINISHED;
	}
}

void traj_loadSeriesData(Trajectory *traj)
{
	int iStep,iData;
	int hopMode = getHopMode();
	for(iData=0;iData<traj->numSeriesData;iData++)
	{	
		FILE *seriesDataFile = fopen(strcat2(traj->dir,*(traj->seriesDataNames+iData)),"r");
		if(seriesDataFile == NULL)
		{
			printf("Could not load series Data file %s\n",strcat2(traj->dir,*(traj->seriesDataNames+iData)));
			continue;
		}
		int numLines = countLines(seriesDataFile);
		printf("seriesDataFile file %s has %d lines\n",strcat2(traj->dir,*(traj->seriesDataNames+iData)),numLines);
		
		for(iStep=0;iStep<fmin(numLines,numSteps);iStep++)//This should probably be replaced by an integer min in case the number of steps exceeds the point where a float can exactly represent an int.
			{
			if(hopMode == HOP_MODE_METROPOLIS)
				fscanf(seriesDataFile,"%*d %lf\n",(traj->seriesData+iStep*traj->numSeriesData+iData));
			else if(hopMode == HOP_MODE_KINETIC)//Consider storing timeSteps here in the future.
				fscanf(seriesDataFile,"%*d %lf %lf\n",(traj->timeseries+iStep),(traj->seriesData+iStep*traj->numSeriesData+iData));
			}
		fclose(seriesDataFile);
	}
}


void traj_saveSeriesData(Trajectory *traj)
{
	int iStep,iData;
	int hopMode = getHopMode();
	//printf("Saving %d series data files\n",traj->numSeriesData);
	for(iData=0;iData<traj->numSeriesData;iData++)	
	{
		//printf("Saving series data file %s\n",*(traj->seriesDataNames+iData));
		FILE *seriesDataFile = fopen(strcat2(traj->dir,*(traj->seriesDataNames+iData)),"w");
		for(iStep=0;iStep<traj->step;iStep++)
		{
		if(hopMode == HOP_MODE_METROPOLIS)
			fprintf(seriesDataFile,"%d %.15g\n",iStep,*(traj->seriesData+iStep*traj->numSeriesData+iData));
		else if(hopMode == HOP_MODE_KINETIC)
			fprintf(seriesDataFile,"%d %.15g %.15g\n",iStep,*(traj->timeseries+iStep),*(traj->seriesData+iStep*traj->numSeriesData+iData));
		}
		fclose(seriesDataFile);
	}
	
}

void traj_loadMoves(Trajectory *traj)
{
	int i;
	FILE *movesFile = fopen(strcat2(traj->dir,"moves"),"r");
	if(movesFile == NULL)
	{
		printf("Could not load moves file %s\n",strcat2(traj->dir,"moves"));
		return;
	}
	int numLines = countLines(movesFile);
	printf("Moves file has %d lines\n",numLines);
	for(i=0;i<fmin(numLines,numSteps);i++)//This should probably be replaced by an integer min in case the number of steps exceeds the point where a float can exactly represent an int.
            fscanf(movesFile,"%*d %d %d\n",(traj->selectedMoves+2*i),(traj->selectedMoves+2*i+1));
	
	traj->step=fmin(numLines,numSteps);
	
	fclose(movesFile);
}

void traj_saveMoves(Trajectory *traj)
{
	int i;
	FILE *movesfile = fopen(strcat2(traj->dir,"moves"),"w");
	for(i=0;i<traj->step;i++)
            fprintf(movesfile,"%d %d %d\n",i,*(traj->selectedMoves+2*i),*(traj->selectedMoves+2*i+1));
	fclose(movesfile);
}

void traj_loadEnergy(Trajectory *traj)
{
	int i;
	int hopMode = getHopMode();
	
	FILE *energyfile = fopen(strcat2(traj->dir,"energy"),"r");
	if(energyfile == NULL)
	{
		printf("Could not load energy file %s\n",strcat2(traj->dir,"energy"));
		return;
	}
	int numLines = countLines(energyfile);
	int numColumns = countColumns(energyfile);
	printf("Energy file has %d lines and %d columns\n",numLines,numColumns);
	for(i=0;i<fmin(numLines,numSteps);i++)//This should probably be replaced by an integer min in case the number of steps exceeds the point where a float can exactly represent an int.
	{
        if(hopMode == HOP_MODE_METROPOLIS)
            fscanf(energyfile,"%*d %lf\n",(traj->energies+i));
        else if(hopMode == HOP_MODE_KINETIC)//3 or 4 columns are allowed for Kinetic since the timesteps are saved.
        {
			//The numColumns == 3option should eventually be removed. At present I have some data that was stored with the old version of traj_saveEnergy.
			if(numColumns == 3)
			{
				fscanf(energyfile,"%*d %lf %lf\n",(traj->timeseries+i),(traj->energies+i));
				if(i != 0)
					*(traj->timeSteps+i) = *(traj->timeseries+i)-*(traj->timeseries+i-1);
			}
			if(numColumns == 4)
				fscanf(energyfile,"%*d %lf %lf %lf\n",(traj->timeseries+i),(traj->timeSteps+i),(traj->energies+i));
		}
		//printf("read energy %g\n",*(traj->energies+i));
     }       
   traj->step=fmin(numLines,numSteps);//This is repeated between loadEnergy and loadMoves but not loadSeriesData. 
   //This also prevents using the "final.xyz" file because "final.xyz" always corresponds to step numLines, but if the numSteps is less than the numLine, this does not work.
            
   fclose(energyfile);         
}

void traj_saveEnergy(Trajectory *traj)
{
	int i;
	int hopMode = getHopMode();
	FILE *energyfile = fopen(strcat2(traj->dir,"energy"),"w");
	for(i=0;i<traj->step;i++)
        if(hopMode == HOP_MODE_METROPOLIS)
            fprintf(energyfile,"%d %.15g\n",i,*(traj->energies+i));
        else if(hopMode == HOP_MODE_KINETIC)
            //fprintf(energyfile,"%d %.15g %.15g\n",i,*(traj->timeseries+i),*(traj->energies+i));
            fprintf(energyfile,"%d %.15g %.15g %.15g\n",i,*(traj->timeseries+i),*(traj->timeSteps+i),*(traj->energies+i));
	fclose(energyfile);
}

//This does not have all it needs to free trajectory
void traj_free(Trajectory *traj)
{
	//printf("Freeing run arrays\n");
	free(traj->energies);
	free(traj->timeseries);
	free(traj->timeSteps);
	free(traj->breakpoints);
	free(traj->seriesData);
}
/*
void traj_copySettings(Trajectory *dest)
{
	if(dest==settingsTrajectory)
	return;
	//dest->numSteps=settingsTrajectory->numSteps;
	//dest->hopMode=settingsTrajectory->hopMode;
	//dest->savecrys=settingsTrajectory->savecrys;
	//dest->temperature=settingsTrajectory->temperature;
	dest->numSeriesData=settingsTrajectory->numSeriesData;
	//dest->savedata=settingsTrajectory->savedata;
	
}*/

void traj_identifySelf(Trajectory *traj, int iTraj, char *dir, int willBeLoaded)
{
	traj->dir = calloc(1000,sizeof(char));
    sprintf(traj->dir,"%s/%03d/",dir,iTraj);
    mkdir2(traj->dir);
	mkdir2(strcat2(traj->dir,"/traj"));
	traj->iTraj = iTraj;
	traj->willBeLoaded = willBeLoaded;

}

 void traj_initialize(Trajectory *traj)
{
	//traj_copySettings(traj);
	traj_generator(traj);
		
	traj->state=TRAJECTORY_STATE_RUNNING;

	printf("numsteps %d\n",numSteps);
	
	LatticeDynamics *LD = traj->LD;
	
	int iHop,maxMoves=0;
	for(iHop=0,maxMoves=0;iHop<traj->LD->numHopPairs;iHop++)
	maxMoves += crys_elementCount(traj->crys,LD->hopPairs+2*iHop*namelength);
	
	maxMoves *= maxconnections;
	MarkovChainMonteCarlo *mcmc = traj->mcmc = malloc(sizeof(MarkovChainMonteCarlo));
	mcmc_allocate(mcmc,maxMoves);
	
	//cn_printAdjacencyList(traj->crys);
	
	traj->energies = calloc(numSteps,sizeof(double));
	traj->timeseries = calloc(numSteps,sizeof(double));
	traj->timeSteps = calloc(numSteps,sizeof(double));
	traj->selectedMoves = calloc(2*numSteps,sizeof(double));
	traj->numrejectedhops=0;
    
    *(traj->timeseries) = 0;
    *(traj->timeSteps) = 0;
    traj->step = 0;
    traj->breakpoints = malloc(numSteps*sizeof(int));//lil overkill but it works
	traj->numbreakpoints = 0;
	
    switch(breakpointmethod)
	{
		case ALL_BREAKPOINTS:
		traj_allBreakpoints(traj);
		break;
		case LOGARITHMIC_BREAKPOINTS:
//This function has been disabled. Maybe later it will be restored.			
		//traj_logarithmicBreakpoints(traj);
		break;
		case BIAS_BREAKPOINTS:
//This function has been disabled. Maybe later it will be restored.		
		//traj_biasBreakpoints(traj);
		break;
	}
	
}

//Trajectory * traj_getSettings()
//{
//	return settingsTrajectory;
//}
