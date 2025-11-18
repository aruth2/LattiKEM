#include "lattikem.h"

char dir[1000];
int numRuns;
int hardcodedMaxThreads;//If this is set it will override the use of the number of logical cores to determine number of threads.
int maxThreads;//The maximum number of threads to use in the threadpool. This does not count control threads. Typically the number of logical cores. 	
//int saveCrys,saveData, saveLightCube;
int saveData;
int jobType;
int replaceMolecules;
int savePeriod;

//For Post-Processing
int averagingSteps;
int equivocationStep;
int averagingMethod;
int equivocationMethod;

int saveCrys;
int saveCrysMode;
int saveCrysStride;
int *saveCrysList;
int numSaveSteps;
//int lightCubeFrameRepeat;
//int lightCubeFrameStride;

void lattikem_startupMessage()
{
	char *message = 
    "  _                                                   _     _               ____        ____                                     \n"
    " | |                        _           _            | |   / |  _________   |   \\      /   |                                    \n"
    " | |                       | |         | |           | |  / /   | _______|  | |\\ \\    / /| |                                             \n"
    " | |                       | |         | |           | | / /    | |         | | \\ \\  / / | |                                       \n"
    " | |                    ___| |___   ___| |___        | |/ /     | |         | |  \\ \\/ /  | |                                       \n"
    " | |                   |____ ____| |____ ____|  O    |   /      | |______   | |   \\  /   | |                                    \n"
    " | |           ___   _     | |         | |      _    |   \\      | _______|  | |    \\/    | |                                            \n"
    " | |          / __ \\| |    | |         | |     | |   | |\\ \\     | |         | |          | |                                     \n"
    " | |         | |  | | |    | |         | |     | |   | | \\ \\    | |         | |          | |                                     \n"
    " | |_______  | |__| | |    | |         | |     | |   | |  \\ \\   | |______   | |          | |                                          \n"
    " |_________|  \\____/|_|    |_|         |_|     |_|   |_|   \\_|  |________|  |_|          |_|                                               \n"
	"\n"
	"A software package for the simulation of lattice chemistry by kinetic monte carlo methods\n"
	"\n"
	"Copyright 2016-2025\n"
	"Anthony Ruth\n"
	"\n";
	printf(message);
}

char *getDir()
{
	return dir;
}

void lattikem_registerSettings()
{
	lattikem_startupMessage();
	pkmc_registerSettings();
	
    registerString(dir,"dir",".");
	registerInt(&numRuns,"numRuns",1);
	registerInt(&hardcodedMaxThreads,"hardcodedMaxThreads",0);
	registerInt(&replaceMolecules,"replaceMolecules",0);
	registerInt(&averagingSteps,"averagingSteps",1);
	registerInt(&saveCrys,"saveCrys",0);
	//registerInt(&saveLightCube,"saveLightCube",0);
	registerInt(&saveData,"saveData",1);
	registerInt(&savePeriod,"savePeriod",1000);
	registerInt(&equivocationStep,"equivocationStep",0);
	registerInt(&saveCrysStride,"saveCrysStride",1);
	saveCrysList = registerIntList("saveCrysList",&numSaveSteps);
	//registerInt(&lightCubeFrameRepeat,"lightCubeFrameRepeat",1);
	//registerInt(&lightCubeFrameStride,"lightCubeFrameStride",1);
		
	registerEnum(6,&jobType,"job",JOB_AUTO,"auto","restart","continue","postProcess","main", "setup");
	registerEnum(3,&averagingMethod,"averagingMethod", TIME_AVERAGE, "step","time","none");
	registerEnum(2,&equivocationMethod,"equivocationMethod", STEP_SUM, "stepSum","timeElapsed");
	registerEnum(2,&saveCrysMode,"saveCrysMode", CRYS_SAVE_STRIDE, "stride","list");
		
	traj_registerSettings();
	
	
	//srand(time(0)) should preferably be used for pseudorandom numbers. However, for debugging using a specific value can be helpful.
	srand(time(0));	
	//srand(5);
}

void performTrajectory(Trajectory *traj)
{
	int hopMode = getHopMode();
	//LatticeDynamics *LD = traj->LD;
	int numSteps = getNumSteps();
    printf("\nStarting trajectory %d from step %d Performing trajectory to %d steps\n\n",traj->iTraj,traj->step,numSteps);
    crystal *crys = traj->crys;
	//int stateSetter;

	//Evaluate the energy of initial structure
	energyjob(&(traj->mcmc->initialenergy),traj->iTraj,traj->step,0,0,NULL);
	
	threadwait(traj->iTraj);
	
	for(;(traj->step)<(numSteps);(traj->step)++)
	{
        //pause the main thread until worker threads have finished their jobs and freed memory
		if((traj->step % savePeriod) == 0 && traj->step != 0)
		{
			//printf("waiting to clear job queue\n");
			threadwait(traj->iTraj);
			traj_saveSeries(traj);
		}
        
        //Advance the system one step forward in time
		if(hopMode == HOP_MODE_KINETIC)
            traj_parallelKineticHopping(traj);
        
        *(traj->selectedMoves+2*(traj->step)) = *(traj->mcmc->moves+2*traj->mcmc->chosenMove); 
        *(traj->selectedMoves+2*(traj->step)+1) = *(traj->mcmc->moves+2*traj->mcmc->chosenMove+1); 
        //printf("Move is %d %d\n",*(traj->selectedMoves+2*(traj->step)),*(traj->selectedMoves+2*(traj->step)+1));
        
        if(traj->step != 0)
        {
			*(traj->timeSteps+traj->step) = (traj->mcmc->timestep);
			*(traj->timeseries+traj->step) = *(traj->timeseries+traj->step-1) + (traj->mcmc->timestep);
		}
        traj->mcmc->initialenergy = *(traj->energies+traj->step);
	}
    
	if(hopMode == HOP_MODE_KINETIC)
        threadwait(traj->iTraj);
	if(hopMode == HOP_MODE_METROPOLIS)
        printf("Trajectory finished with %d rejected hops\n",traj->numrejectedhops);
	//system("date\n");
	
    traj_saveSeries(traj);
	crys_makexyz(crys,strcat2(traj->dir,"final.xyz"));

	traj->state=TRAJECTORY_STATE_FINISHED;
	printf("\nTrajectory %d finished \n",traj->iTraj);
}

//This should be renamed
void simulateTrajectories()
{
	printf("\n\nSimulating Trajectories\n\n");
	system("date\n");		
	mkdir2(dir);
    mkdir2(strcat2(dir,"traj"));
	
	//Disabling this will save a lot of time by not saving the FA in full
	//Also, if the molecules are saved there is currently no way to read the the structure back in.
	if(replaceMolecules)
	perovskite_initializeMolecules();
	
	Trajectory *trajectories = (Trajectory *)malloc(numRuns*sizeof(Trajectory));
	pthread_t *threads = (pthread_t *)malloc(numRuns*sizeof(pthread_t));
	
	
	if(hardcodedMaxThreads)
		maxThreads = hardcodedMaxThreads;
	else
		maxThreads = sysconf(_SC_NPROCESSORS_ONLN);
	printf("Num threads %d\n",maxThreads);
	
	traj_parallelInitialization(trajectories, numRuns, maxThreads, dir, (jobType == JOB_AUTO) || (jobType == JOB_CONTINUE) || (jobType == JOB_POSTPROCESS) || (jobType == JOB_MAIN) || (jobType == JOB_SETUP));
	startWorkerThreadPool(maxThreads, trajectories, numRuns);
	
	int irun,numactiveruns=0;
	int iTraj=0;
	if((jobType == JOB_AUTO) || (jobType == JOB_CONTINUE) || (jobType == JOB_RESTART) || (jobType == JOB_MAIN))//Only case which does not perform trajectories is JOB_POSTPROCESS
	while(numactiveruns > 0 || iTraj < numRuns)
	{
		numactiveruns = 0;
		for(irun=0;irun<iTraj;irun++)
		{
			if((trajectories+irun)->state == TRAJECTORY_STATE_RUNNING)
			numactiveruns++;
		}
	//	printf("Num active threads %d\n",numactivesims);
		while(numactiveruns < maxThreads && iTraj<numRuns)
		{
			if((trajectories+iTraj)->state == TRAJECTORY_STATE_RUNNING)
			{
			//printf("Starting run %d\n",iTraj);
			//system("date\n");
            pthread_create(threads + iTraj,NULL,(void *(*)(void *))performTrajectory,trajectories+iTraj);
			numactiveruns++;
			}
			iTraj++;
		}
		sleep(1);
	}

	printf("All runs finished\n");
			system("date\n");
	if((jobType == JOB_AUTO) || (jobType == JOB_CONTINUE) || (jobType == JOB_RESTART) || (jobType == JOB_POSTPROCESS))
		postProcess(trajectories);		
			
}

void postProcess(Trajectory *trajectories)
{
	//This only saves the structure of the first run. In the future, this may be updated to combine the crystals.
	if(saveCrys)
	{
		printf("Saving Trajectory Crystals\n");
		saveTrajectoryCrystals(trajectories);
	}
	//if(saveLightCube)
	//saveTrajectoryLightCube(trajectories);
	if(saveData)
	{
		printf("Saving Trajectory Data\n");
	if(averagingMethod == STEP_AVERAGE || averagingMethod == NOT_COMBINED)
		stepAverage(trajectories);
	if(averagingMethod == TIME_AVERAGE)
		timeAverage(trajectories);
	}
}

void saveTrajectoryCrystals(Trajectory *traj)
{
	Configuration *config = (Configuration *)malloc(sizeof(Configuration));
	config->crys=traj->crys;
	int iStep;
	int numSteps = getNumSteps();
	for(iStep=0;iStep<(numSteps);iStep++)
	{
		if(iStep==0)
		//traj_crysAtStep(traj->crys,traj->selectedMoves,traj->step-1,iStep,CN_NO_NETWORK);
		traj_crysAtStep(traj->crys,traj->selectedMoves,traj->step,iStep,CN_NO_NETWORK);
		else
		traj_crysAtStep(traj->crys,traj->selectedMoves,iStep-1,iStep,CN_NO_NETWORK);
		switch(saveCrysMode)
		{
			case(CRYS_SAVE_STRIDE):
			if(iStep % saveCrysStride == 0)
				config_savecrys(config, traj, iStep, numSteps);
			break;
			
			case(CRYS_SAVE_LIST):
			if(intcontains(saveCrysList,iStep,numSaveSteps))
				config_savecrys(config, traj, iStep, numSteps);
			break;	
		}
	}
}

//This code is for a cube of LEDs that shows a LattChem trajectory. The code was only used once to program a cube, and will not be supported going forward.
/*
void saveTrajectoryLightCube(Trajectory *traj)
{	
	CrystalConverter *CC = calloc(1,sizeof(CrystalConverter));
	CC->elements = crys_elementString(5,FA,"Pb","I",VA,VX);
	CC->colors = calloc(5*3,sizeof(int));

	CC->offset[0] = 0;
	CC->offset[1] = 4;
	CC->offset[2] = 0;
	*(CC->colors+0*3+0) = 63;
	*(CC->colors+0*3+1) = 63;
	*(CC->colors+1*3+0) = 63;
	*(CC->colors+1*3+2) = 63;
	*(CC->colors+2*3+1) = 63;
	*(CC->colors+2*3+2) = 63;
	
	*(CC->colors+3*3+2) = 255;
	*(CC->colors+4*3+0) = 255;

	CC->numElements = 5;
	CC->cellLength = 3.181; //FAPbI3 lattice constant divided by 2

	Frame frame;
	char filename[1000];
	sprintf(filename,"%s/traj/perovskite.CUBE12",traj->dir);
	FILE *outfile = fopen(filename,"w");

	Configuration *config = malloc(sizeof(Configuration));
	
	config->crys=traj->crys;
	int iStep,iRepeat;
	for(iStep=0;iStep<(traj->numSteps);iStep++)
	{
		if(iStep==0)
		traj_crysAtStep(traj->crys,traj->selectedMoves,traj->step,iStep,CN_NO_NETWORK);
		else
		traj_crysAtStep(traj->crys,traj->selectedMoves,iStep-1,iStep,CN_NO_NETWORK);
		crystalToFrame(traj->crys,&frame,CC);
		
		if((iStep % lightCubeFrameStride) != 0)
		continue;
		for(iRepeat=0;iRepeat<lightCubeFrameRepeat;iRepeat++)
		printFrame(outfile,&frame);
	}
	fclose(outfile);
}*/



void stepAverage(Trajectory *trajectories)
{
	/* This averages trajectories so that trajectories on the same step number are combined into one.
	 * First multiple steps for the same trajectory are averaged, the the combined steps of all trajectories are averaged
	 * */
	int numSteps = getNumSteps();
	Trajectory *combinedTrajectory = (Trajectory *)malloc(sizeof(Trajectory));
	int iMod,iStep,iStepMod,iTraj;

	int modulationSteps = (numSteps)/averagingSteps;
	printf("Post Processing\n");
	printf("There will be %d modulation steps of %d each\n",modulationSteps,averagingSteps);

	printf("Allocating space for post-processing\n");
	Configuration *singleStepConfigs = (Configuration *)malloc(averagingSteps*sizeof(Configuration));
	Configuration *stepAveragedConfigs = (Configuration *)malloc(numRuns*sizeof(Configuration));
	Configuration *runAveragedConfig = (Configuration *)malloc(sizeof(Configuration));
	
	if(averagingMethod == STEP_AVERAGE)
	{
		traj_identifySelf(combinedTrajectory,0,dir,0);
		traj_initialize(combinedTrajectory);
		combinedTrajectory->dir=dir;
		combinedTrajectory->step=trajectories->step;
	}

	for(iStepMod=0;iStepMod<averagingSteps;iStepMod++)
	{
		(singleStepConfigs+iStepMod)->threadnumber=-1;//Do not bind this thread to a specific core/GPU
		config_allocate(singleStepConfigs+iStepMod,combinedTrajectory);
		(singleStepConfigs +iStepMod)->savedata = 1;//Set flag so configs know to record all data that is to be saved/plotted.
	}			
	for(iTraj=0;iTraj<numRuns;iTraj++)
	{
		(stepAveragedConfigs+iTraj)->threadnumber=-1;//Do not bind this thread to a specific core/GPU
		config_allocate(stepAveragedConfigs+iTraj,trajectories+iTraj);
	}
	runAveragedConfig->threadnumber=-1;//Do not bind this thread to a specific core/GPU	
	config_allocate(runAveragedConfig,combinedTrajectory);
		
	
	printf("Combining runs\n");
	for(iMod=0;iMod<modulationSteps;iMod++)
	{
		for(iTraj=0;iTraj<numRuns;iTraj++)
		{
			for(iStepMod=0;iStepMod<averagingSteps;iStepMod++)
			{
				iStep = iStepMod+iMod*averagingSteps;
				//Not sure if this is necessary or kinetic thread will set properly.
				//(singleStepConfigs+iStepMod)->step = iStep;
				energyjob((trajectories+iTraj)->energies+iStep,iTraj,iStep,0,0,singleStepConfigs+iStepMod);
			}
			threadwait(iTraj);//This will be poorly-parallelized when averagingSteps<numthreads
			config_average(singleStepConfigs,averagingSteps,stepAveragedConfigs+iTraj);
		}
		
		if(averagingMethod == STEP_AVERAGE)
		{
			config_average(stepAveragedConfigs,numRuns,runAveragedConfig);
            config_save(runAveragedConfig,combinedTrajectory,iMod,modulationSteps);
		}
		else
			for(iTraj=0;iTraj<numRuns;iTraj++)
				config_save(stepAveragedConfigs+iTraj,trajectories+iTraj,iMod,modulationSteps);
		
	}
		
	if(averagingMethod == STEP_AVERAGE)
	{
		traj_average(trajectories,numRuns,combinedTrajectory);
		//There are no moves here and thus traj_saveSeries should not be called
		traj_saveEnergy(combinedTrajectory);
		traj_saveSeriesData(combinedTrajectory);
	}
	
	for(iTraj=0;iTraj<numRuns;iTraj++)
	{
		traj_saveEnergy(trajectories+iTraj);
		traj_saveSeriesData(trajectories+iTraj);
	}
}

void timeAverage(Trajectory *trajectories)
{
	/* This averages trajectories based on time
	 * First the times are equilibrated at a specific step (time zero)
	 * Then the longest trajectory is identified. It's time points are used
	 * as the time for the combined trajectory.
	 * This follows a different averaging order than the step averaging.
	 * First multiple runs are averaged to the same time, then multiple steps from the
	 * combined trajectory are averaged 
	 * */
	 int numSteps = getNumSteps();
	Trajectory *combinedTrajectory=(Trajectory *)calloc(1,sizeof(Trajectory));
	int modulationSteps = numSteps/averagingSteps;
	
	int iStep, iMod, iStepMod, iTraj, iConfig, longestRunStep, thisRunStep;
	int longestRun=-1;
	double longestRunTime=0;
	double equivocationTime;
	
	printf("Post Processing by time averaging\n");
	printf("There will be %d modulation steps of %d each\n",modulationSteps,averagingSteps);

	printf("Allocating space for post-processing\n");
	printf("Allocating %d configs\n",2*numRuns*averagingSteps);	

	Configuration *singleStepConfigs = (Configuration *)calloc(2*numRuns*averagingSteps,sizeof(Configuration));
	Configuration *timeAveragedConfigs = (Configuration *)calloc(averagingSteps,sizeof(Configuration));
	Configuration *runAveragedConfig = (Configuration *)calloc(1,sizeof(Configuration));
	
	double weights[2*numRuns];
	for(iTraj=0;iTraj<numRuns;iTraj++)
	{
		equivocationTime = *((trajectories+iTraj)->timeseries+equivocationStep);
		
		*((trajectories+iTraj)->timeseries+equivocationStep) = 0;
		
		if(equivocationMethod == STEP_SUM)
		{
			//Set all steps before the equivocation step to negative times. The time series cannot have overlapping times, so this is necessary.
		for(iStep=0;iStep<equivocationStep;iStep++)
			*((trajectories+iTraj)->timeseries+iStep)-= equivocationTime;
			
		for(iStep=equivocationStep+1;iStep<numSteps;iStep++)
		{
			*((trajectories+iTraj)->timeseries+iStep) = *((trajectories+iTraj)->timeseries+iStep-1)  + *((trajectories+iTraj)->timeSteps+iStep);
		}
		}
		if(equivocationMethod == TIME_ELAPSED)
		for(iStep=0;iStep<numSteps;iStep++)
		{
			//When an equivocationTime is used, the timeseries should be contructed by summing the timesteps.
			//This helps to alleviate a numerical precision error from a sudden change in the external conditions. 
			//Moreover, this fix is not an excellent one as the post-processing would need to be re-ran with every sudden change. 
			*((trajectories+iTraj)->timeseries+iStep) -= equivocationTime;
		}
		printf("Run time is %g Equivocation time is %g\n",*((trajectories+iTraj)->timeseries+numSteps-1),equivocationTime);

		if(*((trajectories+iTraj)->timeseries+numSteps-1)>longestRunTime)
		{
			longestRun=iTraj;
			longestRunTime=*((trajectories+iTraj)->timeseries+numSteps-1);
		}
	}
	

	//This has been producing Seg Faults when the above mentioned numerical precision error is occuring. 
	//Time averaging is VERY poorly defined when the numerical precision error occurs. 
	//The only real option to do with that data is run step averaging instead.
	if(longestRun == -1)
	{
		printf("All runs appear to have zero time elapsed.\n Numerical error in timeSeries expected. Try adjusting equivocationStep or rerun with step averaging\n");
		exit(0);
	}
	
	printf("Longest run is %d at time %g Time of first step is %g\n", longestRun,longestRunTime,*((trajectories+longestRun)->timeseries));
		

	traj_identifySelf(combinedTrajectory,0,dir,0);
	traj_initialize(combinedTrajectory);
	combinedTrajectory->dir=dir;
	combinedTrajectory->step=trajectories->step;

	for(iStepMod=0;iStepMod<averagingSteps;iStepMod++)
	{
		(timeAveragedConfigs+iStepMod)->threadnumber=-1;//Do not bind this thread to a specific core/GPU	
		config_allocate(timeAveragedConfigs+iStepMod,combinedTrajectory);
	}			
	for(iTraj=0;iTraj<numRuns;iTraj++)
		{
			for(iStepMod=0;iStepMod<averagingSteps;iStepMod++)
			{
			iConfig = 2*iTraj+2*numRuns*iStepMod;
			(singleStepConfigs+iConfig)->threadnumber=-1;//Do not bind this thread to a specific core/GPU
			(singleStepConfigs+iConfig+1)->threadnumber=-1;//Do not bind this thread to a specific core/GPU
			config_allocate(singleStepConfigs+iConfig,trajectories+iTraj);
			config_allocate(singleStepConfigs+iConfig+1,trajectories+iTraj);
			(singleStepConfigs+iConfig)->savedata = 1;//Set flag so configs know to record all data that is to be saved/plotted.
			(singleStepConfigs+iConfig+1)->savedata = 1;
			}
		}
	(runAveragedConfig)->threadnumber=-1;//Do not bind this thread to a specific core/GPU
	config_allocate(runAveragedConfig,combinedTrajectory);
	
	//The series data appears to get messed up by the energy jobs below. This call here actually uses the data that was loaded from the files to generate the average series data.
	//Unfortunately this prevents regenerating the series data if it is lost, but it should not be for a good run.
	//It would be preferable if the energy jobs below filled in the series data, but not all steps of all runs are evaluated.
	printf("Averaging the trajectories\n");
	traj_weightedAverage(trajectories,numRuns,combinedTrajectory,longestRun);
	//There are no moves here and thus traj_saveSeries should not be called
	traj_saveEnergy(combinedTrajectory);
	traj_saveSeriesData(combinedTrajectory);
	
	
	printf("Combining runs\n");
	for(iMod=0;iMod<modulationSteps;iMod++)
	{
		for(iStepMod=0;iStepMod<averagingSteps;iStepMod++)
		{
			longestRunStep = iStepMod+iMod*averagingSteps;
			longestRunTime = *((trajectories+longestRun)->timeseries+longestRunStep);
			for(iTraj=0;iTraj<numRuns;iTraj++)
			{
				linearInterpolationWeights((trajectories+iTraj)->timeseries, longestRunTime, numSteps, &thisRunStep, weights+2*iTraj, weights+2*iTraj+1);
				//A significant cost savings could be implemented here where the configurations only need to be calculated once if they are at the very beggining or 
				//very end of the series, Since these points get reused.
				if(thisRunStep == numSteps-1)
				{
				thisRunStep--;
				}
				//printf("Performing energy job for run %d steps %d and %d with weights %g and %g based on longestRunStep %d and thisRunStep %d longestRunTime %g and thisRunTime %g\n",iTraj, thisRunStep,thisRunStep+1,*(weights+2*iTraj),*(weights+2*iTraj+1),longestRunStep,thisRunStep,longestRunTime,*((trajectories+iTraj)->timeseries+thisRunStep));

				iConfig = 2*iTraj+2*numRuns*iStepMod;
				energyjob((trajectories+iTraj)->energies+thisRunStep,iTraj,thisRunStep,0,0,singleStepConfigs+iConfig);
				energyjob((trajectories+iTraj)->energies+thisRunStep+1,iTraj,thisRunStep+1,0,0,singleStepConfigs+iConfig+1);
			}
			//threadwait();//This will be poorly-parallelized when numRuns<numthreads
			
		}
		threadwait_all();//This will be much better parallelized for when numRuns<numThreads, however it will use a signficantly higher amount of memory for a large number of runs. Scaling is numRuns*averagingSteps*configSize.
		for(iStepMod=0;iStepMod<averagingSteps;iStepMod++)
		{
			config_weightedAverage(singleStepConfigs+2*numRuns*iStepMod,2*numRuns,timeAveragedConfigs+iStepMod,weights);
		}
		
		//printf("Combining %d steps into a single step\n",averagingSteps);
		config_average(timeAveragedConfigs,averagingSteps,runAveragedConfig);
		//printf("Saving the averaged configuration\n");
        config_save(runAveragedConfig,combinedTrajectory,iMod,modulationSteps);
	}

	
	printf("Time Averaging has finished\n");
}

