/*
 * pkmc.c - a library for parallel kinetic monte carlo
 * in each step all possible moves are enumerated and the enthalpy change of the moves
 * are evaluated in parallel by a pool of worker threads.
 * 
 * */

#include "pkmc.h"

pthread_t *kineticthreads;
int numkineticthreads;
int nextNewKineticJob = 0;
int nextAcquiredKineticJob = 0;
pthread_mutex_t kineticmutex;
pthread_cond_t kineticcond;
pthread_cond_t threadwaitcond;
pthread_mutex_t threadwaitmutex;
int *thread_status;
kineticjob *kineticjobs;
Trajectory *trajectories;
int numTrajectories;


void traj_parallelKineticHopping(Trajectory *traj)
{
	/* Finds and performs a hop, swapping two species. Only hops which involve the allowed hop elements can be performed.
	 * The energy for all possible hops must be calculated.
	 * energy jobs will be popped out to calculate the energy of a particular configuration.
	 * */
	LatticeDynamics *LD = traj->LD;	  	 
    MarkovChainMonteCarlo *mcmc = traj->mcmc;
    crystal *crys = traj->crys;
    double temperature = getTemperature(); 
    
    //double *timestep=&(mcmc->timestep);  
	
    //Step 1 - enumerate all hops; the total energy of the system after the hop is submitted as an energyjob
	
	LD_listMoves(crys,LD,mcmc);
	double *newenergy=(traj->energies+traj->step);
	double initialenergy = mcmc->initialenergy; 
    int *moves = mcmc->moves;
    double *moveRates = mcmc->moveRates;
	double *moveEnthalpies=mcmc->moveEnthalpies;
	double *moveBarriers=mcmc->moveBarriers;
	double *moveProbabilityRanges=mcmc->moveProbabilityRanges;
	
	//mcmc_printMoves(crys, mcmc);
	int numMoves = mcmc->numMoves;
	//printf("There are %d possible moves\n",numMoves);
	int iMove;
	for(iMove=0;iMove<numMoves;iMove++)
	{
		energyjob(moveEnthalpies+iMove,traj->iTraj,traj->step,*(moves+2*iMove),*(moves+2*iMove+1),NULL);
	}
    //Wait until the energy of all post-hop structures has been evaluated
	threadwait();
			
	//For each hop, sum the enthalpic change with the energy barrier.
	mcmc->rate =0;
	int mostFavorableMove=0;
	for (iMove=0;iMove<numMoves;iMove++)
	{
		*(moveRates+iMove) = exp(-(*(moveEnthalpies+iMove)-initialenergy+*(moveBarriers+iMove))/temperature);
		//printf("move %d has an enthalpy of %g, a barrier of %g, a previous energy of %g, and results in a rate of %g\n",iMove,*(moveEnthalpies+iMove),*(moveBarriers+iMove),initialenergy,*(moveRates+iMove));
		if(*(moveEnthalpies+iMove)+*(moveBarriers+iMove)<*(moveEnthalpies+mostFavorableMove)+*(moveBarriers+mostFavorableMove))
			mostFavorableMove = iMove; 

		mcmc->rate += *(moveRates+iMove);
	}
	
	if(!isnormal(mcmc->rate)) // If there is an overflow or underflow in the hop rate set the rate to be the largest value allowed
		{
		//printf("Rates broken, highest rate %g\n",*(moveRates+iMove));
		
		for (iMove=0;iMove<numMoves;iMove++)
		*(moveRates+iMove) = 0;
		
		*(moveRates+mostFavorableMove)=1;
		mcmc->rate=1;
		}
       
	weightToProbabilityRange(moveRates, moveProbabilityRanges, numMoves);
	int chosenMove = mcmc->chosenMove = chooseItem(moveProbabilityRanges,numMoves);
	
	//printf("Chosen move is %d involving swap of %d %s and %d %s \n",chosenMove,*(moves+2*chosenMove),crys->species+*(moves+2*chosenMove)*namelength,*(moves+2*chosenMove+1),crys->species+*(moves+2*chosenMove+1)*namelength);

    //perform the hop
	cn_swap(crys,*(moves+2*(mcmc->chosenMove)),*(moves+2*(mcmc->chosenMove)+1),CN_HAS_NETWORK);
    
	//The energy of the new configuration was already calculated, but we need to calculate the properties of the chosen structure.
    //This can be done asynchronously by adding an energyjob    
	*newenergy = *(moveEnthalpies+(mcmc->chosenMove));
	
	energyjob(newenergy,traj->iTraj,traj->step,0,0,NULL);
	
	//printf("In pkmc energy is %g saving it to %d\n",*newenergy,(int)newenergy);
	mcmc->timestep = 1/mcmc->rate;

}

void threadwait()
{
	//This function will not return until all jobs have been acquired and all kineticthreads have status 0(waiting for work)
	//This may be able to be improved by having a separate function which waits for all worker threads to be finished with a specific trajectory.
	int i;
	int done=0;
	//printf("Waiting for threads to finish\n");
	while(!done)
	{        
	pthread_mutex_lock(&kineticmutex);
	pthread_mutex_lock(&threadwaitmutex);
	for(done=1,i=0;i<numkineticthreads;i++)
	{
		if(*(thread_status+i) == THREAD_WORKING)
		done = 0;
	}
    if(nextNewKineticJob!=nextAcquiredKineticJob)
		done = 0;	 
	pthread_mutex_unlock(&kineticmutex);
	if(!done)
	{
		pthread_cond_wait(&threadwaitcond,&threadwaitmutex);
	}
	pthread_mutex_unlock(&threadwaitmutex);
	}
	
}

void energyjob(double *energy, int iTraj, int step,int additionalMoveAtom1, int additionalMoveAtom2, Configuration *dataSwap)
{
	/* This function adds a new job to the queue for the kinetic threads to handle
	 * */
     
	pthread_mutex_lock(&kineticmutex);
	kineticjob *job = kineticjobs+nextNewKineticJob;
		
	 job->energy = energy;
     job->iTraj=iTraj;
     job->step=step;
     job->additionalMoveAtom1=additionalMoveAtom1;
     job->additionalMoveAtom2=additionalMoveAtom2;
     job->dataSwap=dataSwap;
     
	nextNewKineticJob++;
	if(nextNewKineticJob == maxJobs)
	nextNewKineticJob=0;
	
	//Wake up any waiting threads.
	pthread_cond_signal(&kineticcond);
	
	//printf("Added job %d to the queue\n",numkineticjobs);
	pthread_mutex_unlock(&kineticmutex);
}

void kineticthread(workerData *wd)
{
	/* This thread operates in a threadpool. It acquires a mutex, grabs a job, releases the mutex, and does work
	 * 
	 * */
     int threadnumber = wd->threadnumber;
     Configuration *config;
	 Trajectory *traj;
	 printf("\nStarted worker thread %d\n",threadnumber);
	 kineticjob *job;
	 void *dataHolder;
	 
	 while(1)
	 {
		 pthread_mutex_lock(&kineticmutex);
		*(thread_status+threadnumber) = THREAD_WORKING;
		if(nextNewKineticJob!=nextAcquiredKineticJob)
		 {
			job = (kineticjobs+nextAcquiredKineticJob);
			nextAcquiredKineticJob++;
			if(nextAcquiredKineticJob==maxJobs)
				nextAcquiredKineticJob=0;
			pthread_mutex_unlock(&kineticmutex);
			
			 //printf("Performing Job %d which is step %d on run %d finding energy and pushing it to %d\n",acquiredkineticjobs,job->step,job->run,job->energy);
			
			traj = trajectories+job->iTraj;
			config = wd->workspace+job->iTraj;

			//This sets the save flags so the thread knows whether to calculate the configuration in full, or whether to only calculate enough to get the energy
            if(job->dataSwap != NULL)
			{
				//printf("Swapping data to save\n");
				dataHolder=config->data;
				config->data=job->dataSwap->data;
				config->savedata = job->dataSwap->savedata;//Set the save flag
			}

			traj_crysAtStep(config->crys,traj->selectedMoves,*(wd->trajectoriesteps+job->iTraj),job->step,CN_HAS_NETWORK);
			*(wd->trajectoriesteps+job->iTraj) = (job->step);
            
			config->externalConditions=(traj->externalConditions+traj->numExternalConditions*job->step);


			//This section is used only for moves which are being considered but have not yet been selected.
            if((job->additionalMoveAtom1 != 0) || (job->additionalMoveAtom2 != 0))
            {
            cn_swap(config->crys,job->additionalMoveAtom1,job->additionalMoveAtom2,CN_HAS_NETWORK);
			
			//Never save when it is not the cannonical move
			config->seriesData=NULL;
			}
			else
			config->seriesData=traj->seriesData+traj->numSeriesData*job->step;
			
            *(job->energy) = config_energy(config);
            
            if((job->additionalMoveAtom1 != 0) || (job->additionalMoveAtom2 != 0))
            cn_swap(config->crys,job->additionalMoveAtom1,job->additionalMoveAtom2,CN_HAS_NETWORK);
            
            //Set thread status to zero here so that the main thread does not need to wait for file writing.
            //All known instances of thread_wait only need the energy, not the files, and therefore can continue.
            *(thread_status+threadnumber) = THREAD_IDLE;
            
            if(job->dataSwap != NULL)
				config->data=dataHolder;
			
		 }
		 else
		 {  
			 //printf("Thread %d Sleeping and waiting for work\n",threadnumber);
			*(thread_status+threadnumber) = THREAD_IDLE;
			pthread_mutex_lock(&threadwaitmutex);
			pthread_cond_broadcast(&threadwaitcond);
			pthread_mutex_unlock(&threadwaitmutex);
			pthread_cond_wait(&kineticcond,&kineticmutex);
			pthread_mutex_unlock(&kineticmutex);
		}	
	 }
	
}

void startWorkerThreadPool(int numThreads, Trajectory *newtraj, int newNumTrajectories)
{
	
	trajectories = newtraj;
	numTrajectories = newNumTrajectories;
	
	
	printf("Starting thread pool with %d threads\n",numThreads);
	kineticjobs = malloc(maxJobs*sizeof(kineticjob));
	kineticthreads = malloc(numThreads*sizeof(pthread_t));
	thread_status = malloc(numThreads*sizeof(int));
	numkineticthreads = numThreads;
	pthread_mutex_init(&kineticmutex,NULL);
	pthread_cond_init(&kineticcond,NULL);
	
	pthread_mutex_init(&threadwaitmutex,NULL);
	pthread_cond_init(&threadwaitcond,NULL);
	//printf("Allocated workspace to %d\n",(int)workspace);
	
	int ithread,itraj;
	workerData *wd;
    for(ithread=0;ithread<numThreads;ithread++)
		{
		wd = malloc(sizeof(workerData));
		wd->threadnumber = ithread;
		wd->workspace = malloc(numTrajectories*sizeof(Configuration));
		wd->trajectoriesteps=malloc(numTrajectories*sizeof(int));
   
		for(itraj=0;itraj<numTrajectories;itraj++)
			{
			(wd->workspace+itraj)->crys=crys_duplicate((trajectories+itraj)->crys);
			config_allocate(wd->workspace+itraj,(trajectories+itraj));
			*(wd->trajectoriesteps+itraj)=(trajectories+itraj)->step;
			}
    
		pthread_create(kineticthreads+ithread,NULL,(void *(*)(void *))kineticthread,(void *)wd);
		}
	
}

void singleTrajectoryInitialization(void *data)
{
	Trajectory *traj = (Trajectory *)data;
	int innd;
	
	printf("\n\nInitializing Administrator Thread %d\n\n",traj->iTraj);
	traj_initialize(traj);
	if(traj->willBeLoaded)//Only case which does not load trajectory is JOB_RESTART
		traj_load(traj);
		
	cn_allocateCrystalNetwork(traj->crys);
	for(innd=0;innd<traj->numnnds;innd++)
		cn_fillFromnnd(traj->crys,traj->nnd+innd);
	
	//cn_printAdjacencyList(traj->crys);
		
	printf("Run %d initialized\n",traj->iTraj);	
}

void traj_parallelInitialization(Trajectory *trajectories, int numRuns, int maxThreads, char *dir, int willBeLoaded )
{
	/*
	 * This function initilizes all trajectories in a parallel manner. 
	 * Initialization of trajectories is assigned to threads and each thread performs one
	 * initialization. 
	 * */
	pthread_t *initializationThreads = malloc(numRuns*sizeof(pthread_t));
	int iTraj;
	int numGroups = ceil((double)numRuns/maxThreads);
	printf("\n\nInitializing %d trajectories in %d groups of %d Administrator threads each\n\n",numRuns,numGroups,maxThreads);
	int iGroup;
	for(iGroup = 0;iGroup<numGroups;iGroup++)
	{
		for(iTraj=iGroup*maxThreads;iTraj<numRuns && (iTraj - iGroup*maxThreads) <maxThreads;iTraj++)
		{
			traj_identifySelf((trajectories+iTraj),iTraj,dir,willBeLoaded);
			pthread_create(initializationThreads+iTraj,NULL,(void *(*)(void *))singleTrajectoryInitialization,(void *)(trajectories+iTraj));
		}
		for(iTraj=iGroup*maxThreads;iTraj<numRuns && (iTraj - iGroup*maxThreads) <maxThreads;iTraj++)
			pthread_join(*(initializationThreads+iTraj),NULL);
	}
	
}
