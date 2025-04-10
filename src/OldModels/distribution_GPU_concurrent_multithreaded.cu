//#include "distribution_GPU.h"
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <cuda_runtime.h>
#include <cuda.h>
#include <nvtx3/nvToolsExt.h>
#include <numa.h>
enum MEMCPY_BITMASK {ELEVELS_IN, FENERGY_IN, WEIGHTS_OUT, FENERGY_OUT, ENERGY_OUT};

#define optics_llimit -1.0
#define optics_ulimit 1.0
__device__ void fermiIteration(double *energyLevels, double fermiGuess, double temperature, int numEnergyLevels, double *weights, double *energyWeights, double *weightSum)
{
	double threadEnergy=0;//Sum for each thread
	double threadWeight=0;
	int iLevel; 
	for(iLevel = threadIdx.x;iLevel < numEnergyLevels;iLevel+=blockDim.x)
	{
		weights[iLevel] = 1.0/(1.0+exp((*(energyLevels+iLevel)-fermiGuess)/temperature));
		threadEnergy +=weights[iLevel] * energyLevels[iLevel];
		threadWeight += weights[iLevel];
	}

	energyWeights[threadIdx.x] = threadEnergy; //Store sum from each thread at base value
	weightSum[threadIdx.x] = threadWeight; //Store sum from each thread at base value
	for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) 
	{
		__syncthreads();

		if(threadIdx.x < stride)
		{
			weightSum[threadIdx.x ] += weightSum[threadIdx.x+stride];
			energyWeights[threadIdx.x] += energyWeights[threadIdx.x+stride];
		}
	}
	__syncthreads();
}

__global__ void fermiKernel(double *energyLevels, double *fermiEnergy, double currentExcitations, double fermiConvergence, double temperature, int numEnergyLevels,  double *weights, double *energyWeights, double *d_energy, double *weightSum)
{
	//double fermiGuess = 0.5;
	double fermiGuess = *fermiEnergy;
	double fermiMin = fermiGuess - 7.5*fermiConvergence;
	double fermiMax = fermiGuess + 8.49*fermiConvergence;//This is slightly offset so that the fermi level could be almost unchanged. 
	//double fermiEnergy;
	double totalCarriersMin=0;
	double totalCarriersMax=0; 
	double totalCarriers;
	//Try to bound the fermi level within a small box around the previous result which requires not more than 6 iterations to converge.
	fermiIteration(energyLevels, fermiMin, temperature, numEnergyLevels, weights, energyWeights, weightSum);
	//fermiIterationMultiblock(energyLevels, fermiMin, temperature, numEnergyLevels, weights, energyWeights, weightSum);
	totalCarriersMin = *weightSum;
	fermiIteration(energyLevels, fermiMax, temperature, numEnergyLevels, weights, energyWeights, weightSum);
	//fermiIterationMultiblock(energyLevels, fermiMax, temperature, numEnergyLevels, weights, energyWeights, weightSum);
	totalCarriersMax = *weightSum;
	
	if((totalCarriersMin > currentExcitations) || (totalCarriersMax < currentExcitations))//If the solution is not within the bounds expand the bounds.  	
	{
		//printf("Fermi level outside bounds, Min %.10g max %.10g, min carriers %.10g, max carriers %.10g, target carriers %g\n",fermiMin,fermiMax,totalCarriersMin,totalCarriersMax,currentExcitations);
		fermiMin = optics_llimit;
		fermiMax = optics_ulimit;
	}
	
	while(fermiMax - fermiMin > fermiConvergence)
	{
 
		fermiGuess = (fermiMax+fermiMin)/2.0;	
		
		__syncthreads();
		fermiIteration(energyLevels, fermiGuess, temperature, numEnergyLevels, weights, energyWeights, weightSum);
		//fermiIterationMultiblock(energyLevels, fermiGuess, temperature, numEnergyLevels, weights, energyWeights, weightSum);
		totalCarriers = *weightSum;
	
		if(totalCarriers > currentExcitations) //Lower search range
			fermiMax = fermiGuess;
		else
			fermiMin = fermiGuess;
		//printf("Fermi range %g %g Carriers %g %g\n",fermiMin,fermiMax,totalCarriers,currentExcitations);__syncthreads();
	}
	
	if(threadIdx.x == 0)
	{
		*fermiEnergy = fermiGuess;
		*d_energy = *energyWeights;
	}
}


void fermiSum( double *energyLevels, double fermiEnergy, double temperature, int numEnergyLevels, double *weights, double *weightSum, double *energySum)
{
 	/*
	Evaluates the fermi Sum, N=sum(fi) and E=sum(fi * Ei)
	where fi is the fermi function fi = 1/(1+exp((Ei-Ef)/kT))
	
	Inputs are:
	energyLevels (Ei)
	fermiEnergy (Ef)
	temperature (kT)
	numEnergyLevels (n)
	
	outputs are:
	weights (fi)
	weightSum (N)
	energySum (E)
	*/
	*weightSum=0;
	*energySum=0;
	for(int i = 0;i<numEnergyLevels;i++)
	{
		*(weights+i) = 1.0/(1.0+exp((*(energyLevels+i)-fermiEnergy)/temperature));
		*weightSum += *(weights+i);
		*energySum += *(energyLevels+i) * *(weights+i);
	}
}


double solveFermiEnergy(double currentExcitations, double fermiGuess, double fermiConvergence, double *energyLevels,double *weights, int numEnergyLevels, double temperature, double *energy)
{
	double fermiMin = fermiGuess - 7.5*fermiConvergence;
	double fermiMax = fermiGuess + 8.49*fermiConvergence;//This is slightly offset so that the fermi level could be almost unchanged. 
	//double fermiEnergy;
	double totalCarriersMin=0;
	double totalCarriersMax=0; 
	//Try to bound the fermi level within a small box around the previous result which requires not more than 6 iterations to converge.
	fermiSum( energyLevels, fermiMin, temperature, numEnergyLevels, weights, &totalCarriersMin, energy);
	fermiSum( energyLevels, fermiMax, temperature, numEnergyLevels, weights, &totalCarriersMax, energy);
	
	
	if((totalCarriersMin > currentExcitations) || (totalCarriersMax < currentExcitations))//If the solution is not within the bounds expand the bounds.  	
	{
		//printf("Fermi level outside bounds, Min %.10g max %.10g, min carriers %.10g, max carriers %.10g, target carriers %g\n",fermiMin,fermiMax,totalCarriersMin,totalCarriersMax,currentExcitations);
		fermiMin = optics_llimit;
		fermiMax = optics_ulimit;
	}
		
	double totalCarriers=0;
	while(fermiMax-fermiMin > fermiConvergence)
	{
		fermiGuess = (fermiMax+fermiMin)/2.0;	
		fermiSum( energyLevels, fermiGuess, temperature, numEnergyLevels, weights, &totalCarriers, energy);
		if(totalCarriers > currentExcitations) //Lower search range
			fermiMax = fermiGuess;
		else
			fermiMin = fermiGuess;
		//printf("Fermi range %g %g Carriers %g %g\n",fermiMin,fermiMax,totalCarriers,currentExcitations);
	}
	return fermiGuess;
}

double RandDouble(double low, double high, unsigned int *state) {
	double t = (double)rand_r(state) / (double)RAND_MAX;
	return (1.0f - t) * low + t * high;
}

//https://stackoverflow.com/questions/1407786/how-to-set-cpu-affinity-of-a-particular-pthread
// core_id = 0, 1, ... n-1, where n is the system's number of cores

int stick_this_thread_to_core(int core_id) {
   int num_cores = sysconf(_SC_NPROCESSORS_ONLN);
   if (core_id < 0 || core_id >= num_cores)
      return EINVAL;

   cpu_set_t cpuset;
   CPU_ZERO(&cpuset);
   CPU_SET(core_id, &cpuset);

   pthread_t current_thread = pthread_self();    
   return pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset);
}

// Total number of input vector pairs; 
//const int num_pthreads = 256;
int num_pthreads = 256;
// Number of elements per vector; arbitrary,
// but strongly preferred to be a multiple of warp size
// to meet memory coalescing constraints
int ELEMENT_N = 64000;
//const int ELEMENT_N = 64000;
// Total number of data elements
int sDATA_SZ;
//const int sDATA_SZ = ELEMENT_N * sizeof(double);
int num_blocks = 1;
int num_threads = 128;
int run_mask = 3;
int repeats = 1;
int memcpy_mask = 31;
double fermiConvergence = 1e-7;

#define hyperthreads 2
#define sockets 1
#define nodesPerSocket 4	
void * single_thread(void *data)
{
	int thread_number = *((int *)data);
	stick_this_thread_to_core(thread_number);
	int num_gpu;
	cudaGetDeviceCount(&num_gpu);
	//int iNUMA = (double)(thread_number % 64)/16;
	int cores = num_pthreads/hyperthreads;
	int iNUMA = sockets * nodesPerSocket * ((double)(thread_number % cores))/(cores);
	int iGPU = num_gpu * ((double)(thread_number % cores))/(cores);
	//printf("Thread %d cores %d nodesPerSocket %d NUMA %d\n",thread_number,cores,nodesPerSocket,iNUMA);
	//cudaSetDevice((double)thread_number/num_threads*num_gpu);
	//cudaSetDevice(3-iNUMA);
	cudaSetDevice(iGPU);
	double  *h_weights, *weights_CPU, fermiEnergy_CPU,  energy_CPU;
	double *h_energy, *h_fermiEnergy;
	double *h_energyLevels;
	double *d_energyLevels, *d_weights, *d_energyWeights, *d_weightSum, *d_fermiEnergy, *d_energy;
	int i,iRepeat,ielevel;
	int num_elevel_update = 16;
	int *elevel_update = (int *)malloc(num_elevel_update*sizeof(int));	
	int saccum_size = num_threads*sizeof(double);
	if(thread_number==0)
		printf("Initializing data...\n");
	if(thread_number==0)
		printf("...allocating CPU memory.\n");
	
	//cudaMallocManaged((void**)&h_energyLevels, sDATA_SZ);
	cudaMallocHost((void**)&h_energyLevels, sDATA_SZ);
	//h_energyLevels = (double *)malloc(sDATA_SZ);
	//h_energyLevels = (double *)numa_alloc_onnode(sDATA_SZ,iNUMA);
	//cudaMallocManaged((void**)&h_weights, sDATA_SZ);
	cudaMallocHost((void**)&h_weights, sDATA_SZ);
	//h_weights = (double *)malloc(sDATA_SZ);
	//h_weights = (double *)numa_alloc_onnode(sDATA_SZ,iNUMA);
	//cudaMallocManaged((void**)&h_energy, sizeof(double));
	cudaMallocHost((void**)&h_energy, sizeof(double));
	//h_energy = (double *)malloc(sizeof(double));
	//h_energy = (double *)numa_alloc_onnode(sizeof(double),iNUMA);
	//cudaMallocManaged((void**)&h_fermiEnergy, sizeof(double));
	cudaMallocHost((void**)&h_fermiEnergy, sizeof(double));
	//h_fermiEnergy = (double *)malloc(sizeof(double));
	//h_fermiEnergy = (double *)numa_alloc_onnode(sizeof(double),iNUMA);
	weights_CPU = (double *)malloc(sDATA_SZ);

	if(thread_number==0)
		printf("...allocating GPU memory.\n");
	
	nvtxRangePushA(":CUDAMalloc");
	cudaMalloc((void **)&d_energyLevels, sDATA_SZ);
	cudaMalloc((void **)&d_weights, sDATA_SZ);
	cudaMalloc((void **)&d_energyWeights, saccum_size);
	cudaMalloc((void **)&d_energy, sizeof(double));
	cudaMalloc((void **)&d_weightSum, saccum_size);
	cudaMalloc((void **)&d_fermiEnergy, sizeof(double));
	
  	cudaStreamSynchronize(0);//Finish Mallocs before starting Memcpy
  	nvtxRangePop();
  	
  	if(thread_number==0)
		printf("...initializing data.\n");
	
	unsigned int state;	
	for (i = 0; i < ELEMENT_N; i++) 
	{
		h_energyLevels[i] = RandDouble(0.0f, 1.0f,&state);
	}
	*h_fermiEnergy = 0.5;
	
  
	if((run_mask / 2) % 2 == 1)
	{
		if(thread_number==0)
			printf("...Beginning CPU calculation.\n");
		for(iRepeat=0;iRepeat<repeats;iRepeat++)
		{
			fermiEnergy_CPU =  solveFermiEnergy( ELEMENT_N/2, fermiEnergy_CPU, fermiConvergence, h_energyLevels,weights_CPU,ELEMENT_N,0.025,&energy_CPU);
		
		}
		if(thread_number==0)	
			printf("Fermi Level from CPU is %g energy %g\n",fermiEnergy_CPU,energy_CPU); 
	}

	if ((run_mask % 2 == 1))
	{
		if(thread_number==0)
		printf("...copying input data to GPU mem.\n");
		
		// Copy options data to GPU memory for further processing
		cudaMemcpyAsync(d_energyLevels, h_energyLevels, sDATA_SZ, cudaMemcpyHostToDevice,0);
		//cudaMemcpyAsync(d_fermiEnergy, h_fermiEnergy, sizeof(double), cudaMemcpyHostToDevice,0);
  		cudaStreamSynchronize(0);
  		
  		if(thread_number==0)
		printf("...Running GPU Kernel.\n");
			
		for(iRepeat=0;iRepeat<repeats;iRepeat++)
		{
			if((memcpy_mask >> ELEVELS_IN) %2)
			{
				//Generate random list of updates outside NVTX timing
				for(ielevel=0;ielevel<num_elevel_update;ielevel++)
					elevel_update[ielevel]=RandDouble(0,1,&state);
				nvtxRangePushA(":ELEVELS_IN");
				//for(ielevel=0;ielevel<num_elevel_update;ielevel++) 
					//cudaMemcpyAsync(d_energyLevels+elevel_update[ielevel], h_energyLevels+elevel_update[ielevel], sizeof(double), cudaMemcpyHostToDevice,0);
				cudaMemcpyAsync(d_energyLevels, h_energyLevels, sDATA_SZ, cudaMemcpyHostToDevice,0);
				//cudaMemcpyAsync(d_energyLevels, h_energyLevels, 16*sizeof(double), cudaMemcpyHostToDevice,0);
				cudaStreamSynchronize(0);
				nvtxRangePop();
			}
			
			if((memcpy_mask >> FENERGY_IN) %2)
			{
				nvtxRangePushA(":FENERGY_IN"); 
				cudaMemcpyAsync(d_fermiEnergy, h_fermiEnergy, sizeof(double), cudaMemcpyHostToDevice,0);
				cudaStreamSynchronize(0);
				nvtxRangePop();
			}
			nvtxRangePushA(":KERNEL");	
			fermiKernel<<<1, num_threads,0,0>>>(d_energyLevels, d_fermiEnergy,  ELEMENT_N/2, fermiConvergence , 0.025, ELEMENT_N, d_weights,d_energyWeights,d_energy,d_weightSum);
 			cudaStreamSynchronize(0);		
			nvtxRangePop();
				
			
			if((memcpy_mask >> WEIGHTS_OUT) %2)
			{
				nvtxRangePushA(":WEIGHTS_OUT");
				cudaMemcpyAsync(h_weights, d_weights, sDATA_SZ, cudaMemcpyDeviceToHost,0);
				cudaStreamSynchronize(0);
				nvtxRangePop();
			}
			
			if((memcpy_mask >> FENERGY_OUT) %2)
			{
				nvtxRangePushA(":FENERGY_OUT");
				cudaMemcpyAsync(h_fermiEnergy, d_fermiEnergy,  sizeof(double), cudaMemcpyDeviceToHost,0);
				cudaStreamSynchronize(0); 
				nvtxRangePop();
			}
			
			if((memcpy_mask >> ENERGY_OUT) %2)
			{ 
				nvtxRangePushA( ":ENERGY_OUT");
				cudaMemcpyAsync(h_energy, d_energy, sizeof(double), cudaMemcpyDeviceToHost,0);
				cudaStreamSynchronize(0);
				nvtxRangePop(); 
			}
  		}

  		cudaMemcpyAsync(h_weights, d_weights, sDATA_SZ, cudaMemcpyDeviceToHost,0);
		cudaMemcpyAsync(h_fermiEnergy, d_fermiEnergy,  sizeof(double), cudaMemcpyDeviceToHost,0);  
		cudaMemcpyAsync(h_energy, d_energy, sizeof(double), cudaMemcpyDeviceToHost,0); 	
  		
  		cudaStreamSynchronize(0);	
 
  		//cudaDeviceSynchronize();	
  		if(thread_number==0)
  		printf("Fermi Level GPU is %g Energy GPU is %g\n",*h_fermiEnergy,*h_energy);
  		//printf("Fermi Level GPU is %g Energy GPU is %g\n",*d_fermiEnergy,*d_energy);
  		cudaStreamSynchronize(0);//Finish Memcpys before Cuda_free
  	}
  	
	if(thread_number==0)
  	printf("...freeing GPU memory.\n");
  	nvtxRangePushA(":CUDAFree");
	cudaFree(d_energyLevels);
	cudaFree(d_weights);
	cudaFree(d_energyWeights);
	cudaFree(d_weightSum);
	cudaFree(d_fermiEnergy);
	cudaFree(d_energy);
	cudaStreamSynchronize(0);
	nvtxRangePop(); 
	return 0;
}

int main(int argc, char **argv)
{

	if (argc == 8)
	{
		num_pthreads = atoi(argv[1]);
		num_blocks = atoi(argv[2]);
		num_threads = atoi(argv[3]);
		run_mask = atoi(argv[4]);
		repeats = atoi(argv[5]);
		ELEMENT_N = atoi(argv[6]);
		memcpy_mask = atoi(argv[7]);
	}
	else
		printf("Usage: distribution_GPU  num_pthreads num_blocks num_threads run_mask repeats elements memcpy_mask[0-31]\n");
	
	
	printf("Running with %d pthreads %d blocks %d threads, %d run_mask and %d repeats %d elements\n",num_pthreads,num_blocks,num_threads,run_mask,repeats,ELEMENT_N);
	if((memcpy_mask >> ELEVELS_IN) %2)
		printf("Copying energy levels in\n");
	if((memcpy_mask >> FENERGY_IN) %2)
		printf("Copying fermi energy in\n");
	if((memcpy_mask >> WEIGHTS_OUT) %2)
		printf("Copying weights out\n");
	if((memcpy_mask >> FENERGY_OUT) %2)
		printf("Copying fermi energy out\n");
	if((memcpy_mask >> ENERGY_OUT) %2)
		printf("Copying energy out\n");
	printf("%s Starting...\n\n", argv[0]);
	sDATA_SZ = ELEMENT_N * sizeof(double);

  	int iThread = 0;
  	int *threadNumber = (int *)malloc(num_pthreads*sizeof(int));
  	void *retval;
  	pthread_t *threads = (pthread_t *)malloc(num_pthreads*sizeof(pthread_t));
  	for(iThread=0;iThread<num_pthreads;iThread++)
  	{
  		*(threadNumber+iThread) = iThread;
  		pthread_create(threads+iThread,NULL,(void *(*)(void *))single_thread,(void *)(threadNumber+iThread));
  		//single_thread((void *)(&iThread));
  	}
  	for(iThread=0;iThread<num_pthreads;iThread++)
  	{
  		pthread_join(*(threads+iThread),&retval);
  	}
}


