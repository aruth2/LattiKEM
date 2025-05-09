/*
 * pkmc.c - a library for parallel kinetic monte carlo
 * in each step all possible moves are enumerated and the enthalpy change of the moves
 * are evaluated in parallel by a pool of worker threads.
 * 
 * */

#include "pkmc.h"

pthread_t *kineticthreads;
int numkineticthreads;
int *nextNewKineticJob;
int *nextAcquiredKineticJob;
kineticjob **kineticjobs;
pthread_mutex_t *trajmutex;
pthread_mutex_t kineticmutex;
pthread_cond_t kineticcond;
pthread_cond_t *threadwaitcond;


int *thread_status;
int *thread_traj;
Trajectory *trajectories;
int numTrajectories;
int threadWaitTimeout;

int bindThreads;
int hyperthreads;
int sockets;
int nodesPerSocket;

int coreOffset;

int core_number(int threadnumber)
{
	
	int numCores = sysconf(_SC_NPROCESSORS_ONLN);
	int coresPerSocket = numCores/sockets;
	int coresPerNode = coresPerSocket/nodesPerSocket;
	//int socket = (threadnumber+coreOffset)/coresPerSocket;
	int hyperthread = ((threadnumber+coreOffset) % coresPerNode)/(coresPerNode/hyperthreads);
	
	int node = (threadnumber+coreOffset) / (coresPerNode); //the numa node
	int offset = (threadnumber+coreOffset) % (coresPerNode/hyperthreads);
	int core_id = offset + node*(coresPerNode/hyperthreads) + hyperthread * (coresPerSocket/hyperthreads);
	//int core_id = offset + node*(coresPerNode/hyperthreads) +  socket*coresPerSocket + hyperthread * (coresPerSocket/hyperthreads);
	
	return core_id;
}

int bind_thread_to_core(int threadnumber) {
   
   int num_cores = sysconf(_SC_NPROCESSORS_ONLN);
   int core_id = core_number(threadnumber);
   if (core_id < 0 || core_id >= num_cores)
   {
      printf("Thread %d Core id error %d\n",threadnumber,core_id);
      return 1;
	}
   cpu_set_t cpuset;
   CPU_ZERO(&cpuset);
   CPU_SET(core_id, &cpuset);

   pthread_t current_thread = pthread_self();
   printf("Bound thread %d to core %d\n",threadnumber,core_id);    
   return pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset);
}

double thread_affinity(int threadnumber){
	//Returns thread_affinity between 0 and 1 for the relative enumeration of hardware threads.
	int num_HDW_threads = sysconf(_SC_NPROCESSORS_ONLN);
	return (double)(threadnumber+coreOffset)/(num_HDW_threads);
	//int cores = num_HDW_threads/hyperthreads;
	//int core_id = core_number(threadnumber);
	//return ((double)(threadnumber % cores))/(cores);
}

int get_NUMA_Node(int threadnumber){
	return sockets * nodesPerSocket * thread_affinity(threadnumber);
}

void pkmc_registerSettings()
{
	registerInt(&bindThreads,"bindThreads",1);
	registerInt(&hyperthreads,"hyperthreads",2);
	registerInt(&sockets,"sockets",1);
	registerInt(&nodesPerSocket,"nodesPerSocket",4);
	registerInt(&threadWaitTimeout,"threadWaitTimeout",30);
	registerInt(&coreOffset,"coreOffset",0);
}

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
    
    //Step 1 - enumerate all hops; the total energy of the system after the hop is submitted as an energyjob
	
	LD_listMoves(crys,LD,mcmc);
	double *newenergy=(traj->energies+traj->step);
	double initialenergy = mcmc->initialenergy; 
    int *moves = mcmc->moves;
    double *moveRates = mcmc->moveRates;
	double *moveEnthalpies=mcmc->moveEnthalpies;
	double *moveBarriers=mcmc->moveBarriers;
	double *moveProbabilityRanges=mcmc->moveProbabilityRanges;
	
	int numMoves = mcmc->numMoves;
	//printf("There are %d possible moves\n",numMoves);
	int iMove;
	for(iMove=0;iMove<numMoves;iMove++)
	{
		energyjob(moveEnthalpies+iMove,traj->iTraj,traj->step,*(moves+2*iMove),*(moves+2*iMove+1),NULL);
	}
    //Wait until the energy of all post-hop structures has been evaluated
	threadwait(traj->iTraj);
			
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

void threadwait(int iTraj)
{
	//This function will not return until all jobs have been acquired and all kineticthreads have status 0(waiting for work)
	//This may be able to be improved by having a separate function which waits for all worker threads to be finished with a specific trajectory.
	int i;
	int done=0;
	//printf("Waiting for threads to finish\n");
	while(!done)
	{        
	pthread_mutex_lock(trajmutex+iTraj);
	
	for(done=1,i=0;i<numkineticthreads;i++)
	{
		if((*(thread_status+i) == THREAD_WORKING) && (*(thread_traj+i) == iTraj))
		{
		done = 0;
		//printf("Thread %d working on traj %d, cannot return yet\n",i,iTraj);
		}
	}
    if(*(nextNewKineticJob+iTraj)!=*(nextAcquiredKineticJob+iTraj))
		done = 0;	 
	
	if(!done)
	{
		//printf("waiting for traj %d to wake\n",iTraj);
		
		//printf("waiting for traj %d to wake\n",iTraj);
		
		pthread_mutex_lock(&kineticmutex);
		pthread_cond_broadcast(&kineticcond);
		pthread_mutex_unlock(&kineticmutex);
		
		int retVal;
		struct timespec ts;
		timespec_get(&ts,TIME_UTC);
		ts.tv_sec = ts.tv_sec + threadWaitTimeout; //specify the timeout
		retVal = pthread_cond_timedwait(threadwaitcond+iTraj,trajmutex+iTraj,&ts);
		if(retVal)
		{
			printf("Traj %d errored with retVal %d\n",iTraj,retVal);
		}
		//pthread_cond_wait(threadwaitcond+iTraj,trajmutex+iTraj);
		
		//printf("Traj %d awake\n",iTraj);
	}
	pthread_mutex_unlock(trajmutex+iTraj);
	
	}
	
}

void threadwait_all()
{
	int iTraj;
	for(iTraj=0;iTraj<numTrajectories;iTraj++)
		threadwait(iTraj);
}

void energyjob(double *energy, int iTraj, int step,int additionalMoveAtom1, int additionalMoveAtom2, Configuration *dataSwap)
{
	/* This function adds a new job to the queue for the kinetic threads to handle
	 * */
    //printf("Energy job 0\n"); 
	
    //printf("Energy job for traj %d with %d unacquired jobs\n",iTraj,*(nextNewKineticJob+iTraj)-*(nextAcquiredKineticJob+iTraj)); 
	pthread_mutex_lock(trajmutex+iTraj);
	
	kineticjob *job = *(kineticjobs+iTraj)+*(nextNewKineticJob+iTraj);
		
	job->energy = energy;
    job->iTraj=iTraj;
    job->step=step;
    job->additionalMoveAtom1=additionalMoveAtom1;
    job->additionalMoveAtom2=additionalMoveAtom2;
    job->dataSwap=dataSwap;
     
     
	(*(nextNewKineticJob+iTraj))++;
	if(*(nextNewKineticJob+iTraj) == maxJobs)
		(*(nextNewKineticJob+iTraj))=0;
	
	pthread_mutex_unlock(trajmutex+iTraj);
	//printf("Energy job 2\n"); 
	//Wake up any waiting threads.
	//pthread_cond_signal(&kineticcond);
	
	pthread_mutex_lock(&kineticmutex);
	//This is only enough work for a single thread so only one thread is woken up.
	pthread_cond_signal(&kineticcond);
	//pthread_cond_broadcast(&kineticcond);
	pthread_mutex_unlock(&kineticmutex);

	//printf("Added job %d to the queue\n",numkineticjobs);
	
}

void kineticthread(workerData *wd)
{
	/* This thread operates in a threadpool. It acquires a mutex, grabs a job, releases the mutex, and does work
	 * 
	 * */
     int threadnumber = wd->threadnumber;
     
     if(bindThreads)
		bind_thread_to_core(threadnumber);
     
     int iTraj;
     for(iTraj=0;iTraj<numTrajectories;iTraj++)
	 {
		(wd->workspace+iTraj)->threadnumber = threadnumber;
		config_allocate(wd->workspace+iTraj,(trajectories+iTraj));
     }
     
     Configuration *config;
	 Trajectory *traj;
	 printf("\nStarted worker thread %d\n",threadnumber);
	 kineticjob *job;
	 void *dataHolder;
	 int stepHolder;
	 
	 int jobFound;
	 while(1)
	 {
		jobFound = 0;
		for(iTraj = 0;iTraj<numTrajectories;iTraj++)
		{
		traj = trajectories+iTraj;
		config = wd->workspace+iTraj;
		//This is helpful during kinetic steps, but not during post processing.
			
		if(config->step != traj->step)
		traj_configAtStep(config,traj->selectedMoves,config->step,traj->step,CN_HAS_NETWORK);
		
		pthread_mutex_lock(trajmutex+iTraj);
		//pthread_mutex_lock(trajmutex+iTraj);
		if(*(nextNewKineticJob+iTraj)!=*(nextAcquiredKineticJob+iTraj))
		{
			*(thread_status+threadnumber) = THREAD_WORKING;
			jobFound = 1;
			*(thread_traj+threadnumber) = iTraj;
			job = (*(kineticjobs+iTraj)+*(nextAcquiredKineticJob+iTraj));
			(*(nextAcquiredKineticJob+iTraj))++;
			
			if(*(nextAcquiredKineticJob+iTraj)==maxJobs)
				(*(nextAcquiredKineticJob+iTraj))=0;
			pthread_mutex_unlock(trajmutex+iTraj);
			//printf("Acquired job on traj %d with %d unacquired jobs\n",iTraj,*(nextNewKineticJob+iTraj)-*(nextAcquiredKineticJob+iTraj));
			 //printf("Performing Job %d which is step %d on run %d finding energy and pushing it to %d\n",acquiredkineticjobs,job->step,job->run,job->energy);
			
			//traj = trajectories+job->iTraj;
			//config = wd->workspace+job->iTraj;

			//This sets the save flags so the thread knows whether to calculate the configuration in full, or whether to only calculate enough to get the energy
            if(job->dataSwap != NULL)
			{
				//printf("Swapping data to save\n");
				dataHolder=config->data;
				stepHolder=config->step;
				config->step=job->dataSwap->step;
				config->data=job->dataSwap->data;
				config->savedata = job->dataSwap->savedata;//Set the save flag
			}

			traj_crysAtStep(config->crys,traj->selectedMoves,*(wd->trajectoriesteps+job->iTraj),job->step,CN_HAS_NETWORK);
			traj_configAtStep(config,traj->selectedMoves,config->step,job->step,CN_HAS_NETWORK);
			*(wd->trajectoriesteps+job->iTraj) = (job->step);
            
			config->externalConditions=(traj->externalConditions+traj->numExternalConditions*job->step);


			//This section is used only for moves which are being considered but have not yet been selected.
            if((job->additionalMoveAtom1 != 0) || (job->additionalMoveAtom2 != 0))
            {
            cn_swap(config->crys,job->additionalMoveAtom1,job->additionalMoveAtom2,CN_HAS_NETWORK);
			config_performSwap(config,job->additionalMoveAtom1,job->additionalMoveAtom2,1);
			//Never save when it is not the cannonical move
			config->seriesData=NULL;
			}
			else
			config->seriesData=traj->seriesData+traj->numSeriesData*job->step;
			
            *(job->energy) = config_energy(config);
            
            if((job->additionalMoveAtom1 != 0) || (job->additionalMoveAtom2 != 0))
            {
				cn_swap(config->crys,job->additionalMoveAtom1,job->additionalMoveAtom2,CN_HAS_NETWORK);
				config_performSwap(config,job->additionalMoveAtom1,job->additionalMoveAtom2,-1);

            }

            if(job->dataSwap != NULL)
			{
				config->data=dataHolder;
				job->dataSwap->step = config->step;
				config->step = stepHolder;
			}
			//Set thread status to idle here so that the main thread does not need to wait for file writing.
            //All known instances of thread_wait only need the energy, not the files, and therefore can continue.
            *(thread_status+threadnumber) = THREAD_IDLE;
            
			//printf("Thread %d done with job on traj %d, wakeTraj %d\n",threadnumber,iTraj,wakeTraj);
			//printf("Waking up traj %d\n",iTraj);
			pthread_mutex_lock(trajmutex+iTraj);				
			pthread_cond_broadcast(threadwaitcond+iTraj);
			pthread_mutex_unlock(trajmutex+iTraj);
			
		 }
		 else
		pthread_mutex_unlock(trajmutex+iTraj);
		
		}
		if(!jobFound)
		{  
			
			pthread_mutex_lock(&kineticmutex);
			//printf("waiting for thread %d to wake\n",threadnumber);
			pthread_cond_wait(&kineticcond,&kineticmutex);
			pthread_mutex_unlock(&kineticmutex);
			//printf("Thread %d awake\n",threadnumber);
		}
	}	
}

void startWorkerThreadPool(int numThreads, Trajectory *newtraj, int newNumTrajectories)
{
	
	trajectories = newtraj;
	numTrajectories = newNumTrajectories;
	
	
	printf("Starting thread pool with %d threads\n",numThreads);
	//kineticjobs = malloc(maxJobs*sizeof(kineticjob));
	kineticthreads = (pthread_t *)malloc(numThreads*sizeof(pthread_t));
	thread_status = (int *)malloc(numThreads*sizeof(int));
	thread_traj = (int *)malloc(numThreads*sizeof(int));
	numkineticthreads = numThreads;
	
	nextNewKineticJob = (int *)calloc(numTrajectories,sizeof(int));
	nextAcquiredKineticJob = (int *)calloc(numTrajectories,sizeof(int));
	kineticjobs = (kineticjob **)calloc(numTrajectories,sizeof(kineticjob *));
	trajmutex = (pthread_mutex_t *)malloc(numTrajectories*sizeof(pthread_mutex_t));
	threadwaitcond = (pthread_cond_t *)malloc(numTrajectories*sizeof(pthread_cond_t));
	pthread_mutex_init(&kineticmutex,NULL);
	
	for(int iTraj = 0;iTraj<numTrajectories;iTraj++)
	{
		pthread_mutex_init(trajmutex+iTraj,NULL);
		pthread_cond_init(threadwaitcond+iTraj,NULL);
		*(kineticjobs+iTraj) = (kineticjob *)malloc(maxJobs*sizeof(kineticjob));
		
	}
	pthread_cond_init(&kineticcond,NULL);
	
	//pthread_mutex_init(&threadwaitmutex,NULL);
	//pthread_cond_init(&threadwaitcond,NULL);
	//printf("Allocated workspace to %d\n",(int)workspace);
	
	int ithread,itraj;
	workerData *wd;
    for(ithread=0;ithread<numThreads;ithread++)
		{
		wd = (workerData *)malloc(sizeof(workerData));
		wd->threadnumber = ithread;
		wd->workspace = (Configuration *)malloc(numTrajectories*sizeof(Configuration));
		wd->trajectoriesteps=(int *)malloc(numTrajectories*sizeof(int));
   
		for(itraj=0;itraj<numTrajectories;itraj++)
			{
			(wd->workspace+itraj)->crys=crys_duplicate((trajectories+itraj)->crys);
			(wd->workspace+itraj)->step=(trajectories+itraj)->step;
			//config_allocate(wd->workspace+itraj,(trajectories+itraj));
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
	{
		cn_printnnd(traj->nnd+innd);
		//crys_printElements(traj->crys);
		//crys_printAllAtoms(traj->crys);
		cn_fillFromnnd(traj->crys,traj->nnd+innd);
	}
	//cn_printAdjacencyList(traj->crys,traj->crys->network);
	
	if(traj->step == 0)//Only save the crystal as initial if it is on the first step.
		crys_makexyz(traj->crys,strcat2(traj->dir,"initial.xyz"));
	
	traj_saveExternalConditions(traj);	
	printf("Run %d initialized\n",traj->iTraj);	
}

void traj_parallelInitialization(Trajectory *trajectories, int numRuns, int maxThreads, char *dir, int willBeLoaded )
{
	/*
	 * This function initilizes all trajectories in a parallel manner. 
	 * Initialization of trajectories is assigned to threads and each thread performs one
	 * initialization. 
	 * */
	pthread_t *initializationThreads = (pthread_t *)malloc(numRuns*sizeof(pthread_t));
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
