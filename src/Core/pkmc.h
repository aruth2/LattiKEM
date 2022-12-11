#ifndef _PKMC_H_
#define _PKMC_H_

#include "latticedynamics.h"

enum worker_status {THREAD_IDLE, THREAD_WORKING};

//The workerdata structure contains information unique to a single thread.
//This includes the worker's id so it can declare its state to the main thread.
//The worker is also given workspace for when only an energy, no data is required from it.
typedef struct workerData
{
    int threadnumber;
    Configuration *workspace;
	int *trajectoriesteps;
} workerData;

//A kinetic job is a single task which any worker thread can take on and complete. Following completion of the task,
//The energy variable is filled in, and if requested the data is also filled in. The worker thread may be tasked
//with saving the data and also possibly the crystal following its completion of the job.
typedef struct kineticjob{
	int iTraj;
	int step;
	double *energy;
	Configuration *dataSwap;
	int additionalMoveAtom1, additionalMoveAtom2;
} kineticjob;

#define maxJobs 100000

void traj_parallelKineticHopping(Trajectory *traj);
void threadwait();
void energyjob(double *energy, int iTraj, int step,int additionalMoveAtom1, int additionalMoveAtom2, Configuration *dataswap);
void kineticthread(workerData *wd);
void startWorkerThreadPool(int numthreads, Trajectory *newtraj, int newNumTrajectories);
void traj_parallelInitialization(Trajectory *trajectories, int numRuns, int maxThreads, char *dir, int jobType);
#endif
