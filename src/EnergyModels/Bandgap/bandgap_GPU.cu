#include "bandgap_GPU.h"
int GPU_threads;
int sDATA_SZ;

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
	//printf("Running kernel on GPU\n");
	double fermiGuess = *fermiEnergy;
	double fermiMin = fermiGuess - 7.5*fermiConvergence;
	double fermiMax = fermiGuess + 8.49*fermiConvergence;//This is slightly offset so that the fermi level could be almost unchanged. 

	double totalCarriersMin=0;
	double totalCarriersMax=0; 
	double totalCarriers;
	
	//Try to bound the fermi level within a small box around the previous result which requires not more than 6 iterations to converge.
	fermiIteration(energyLevels, fermiMin, temperature, numEnergyLevels, weights, energyWeights, weightSum);
	totalCarriersMin = *weightSum;
	fermiIteration(energyLevels, fermiMax, temperature, numEnergyLevels, weights, energyWeights, weightSum);
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

__host__ void bgg_runKernel(OptoelectronicState *OS, int copyWeights, int num_states, double currentExcitations)
{
	double temperature = getTemperature();
	
	nvtxRangePushA(":ELEVELS_IN");
	cudaMemcpyAsync(OS->d_energyLevels, OS->energystates, sDATA_SZ, cudaMemcpyHostToDevice,0);
	//cudaStreamSynchronize(0);
	nvtxRangePop();
	
	nvtxRangePushA(":KERNEL");	
	fermiKernel<<<1, GPU_threads,0,0>>>(OS->d_energyLevels, OS->d_fermiEnergy,  currentExcitations, getFermiConvergence() , temperature, num_states, OS->d_weights,OS->d_energyWeights,OS->d_energy,OS->d_weightSum);
 	//cudaStreamSynchronize(0);		
	nvtxRangePop();
	
	if(copyWeights)
	{
	nvtxRangePushA(":WEIGHTS_OUT");
	cudaMemcpyAsync(OS->weights, OS->d_weights, sDATA_SZ, cudaMemcpyDeviceToHost,0);
	//cudaStreamSynchronize(0);
	nvtxRangePop();
	}
	
	nvtxRangePushA(":FENERGY_OUT");
	cudaMemcpyAsync(OS->fermiEnergy, OS->d_fermiEnergy,  sizeof(double), cudaMemcpyDeviceToHost,0);
	//cudaStreamSynchronize(0); 
	nvtxRangePop();
	
	nvtxRangePushA( ":ENERGY_OUT");
	cudaMemcpyAsync(OS->photocarrierEnergy, OS->d_energy, sizeof(double), cudaMemcpyDeviceToHost,0);
	cudaStreamSynchronize(0);
	nvtxRangePop(); 
}

__host__ void bindGPU(int threadnumber)
{
	if(threadnumber == -1)
	return;
	
	int num_gpu;
	cudaGetDeviceCount(&num_gpu);
	int iGPU = num_gpu * thread_affinity(threadnumber);
	cudaSetDevice(iGPU);
}

__host__ void bgg_allocate(Configuration *config, int num_states)
{
	//printf("Allocating on GPU\n");
	OptoelectronicState *OS = (OptoelectronicState *)(config->data);
	bindGPU(config->threadnumber);
	
	int saccum_size = GPU_threads*sizeof(double);
	sDATA_SZ = num_states * sizeof(double);
	nvtxRangePushA("Allocate");
	cudaMalloc((void **)&(OS->d_energyLevels), sDATA_SZ);
	cudaMalloc((void **)&(OS->d_weights), sDATA_SZ);
	cudaMalloc((void **)&(OS->d_energyWeights), saccum_size);
	cudaMalloc((void **)&(OS->d_energy), sizeof(double));
	cudaMalloc((void **)&(OS->d_weightSum), saccum_size);
	cudaMalloc((void **)&(OS->d_fermiEnergy), sizeof(double));
	
	cudaMallocHost((void**)&(OS->energystates), sDATA_SZ);
	cudaMallocHost((void**)&(OS->weights), sDATA_SZ);
	cudaMallocHost((void**)&(OS->photocarrierEnergy), sizeof(double));
	cudaMallocHost((void**)&(OS->fermiEnergy), sizeof(double));
	
	cudaStreamSynchronize(0);//Finish Mallocs before starting Memcpy
  	nvtxRangePop();
}

__host__ void bgg_free(OptoelectronicState *OS)
{
	cudaFree(OS->d_energyLevels);
	cudaFree(OS->d_weights);
	cudaFree(OS->d_energyWeights);
	cudaFree(OS->d_energy);
	cudaFree(OS->d_weightSum);
	cudaFree(OS->d_fermiEnergy);
	
	cudaFree(OS->energystates);
	cudaFree(OS->weights);
	cudaFree(OS->photocarrierEnergy);
	cudaFree(OS->fermiEnergy);
	
}

__host__ void bgg_registerSettings()
{
	registerInt(&GPU_threads,"GPU_threads",1024);
	//registerDouble(&fermiConvergence,"fermiConvergence",1e-5);
}
