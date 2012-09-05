//fir filtering via fft with cuda

//includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>


// includes, project
#include <cufft.h>
#include <cutil.h>

// Complex data type
typedef float Complex[2]; 
typedef float Real;
static __device__ __host__ inline Complex ComplexAdd(Complex, Complex);
static __device__ __host__ inline Complex ComplexScale(Complex, float);
static __device__ __host__ inline Complex ComplexMul(Complex, Complex);
static __global__ 			void  ComplexPointwiseMulAndScale(Complex*, const Complex*, int, float);

// Filtering functions
void fftFilter(float* hsignal, float* kernel, int n)
{// filteredSignal will hold the convolution of signal and kernel, all of which have n elements.
	
    int memSizeReal    = sizeof(Real) * n;
	int memSizeComplex = sizeof(Complex) * (n/2+1);
 
    // Allocate device memory for signal
    Real* dSignal;
    CUDA_SAFE_CALL(cudaMalloc((void**)&dSignal, memSizeReal));
    // Copy host memory to device
    CUDA_SAFE_CALL(cudaMemcpy(dSignal, hSignal,memSizeReal,
                              cudaMemcpyHostToDevice));

    // Allocate device memory for filter kernel signal, and transforms
    Real* dKernel;
	CUDA_SAFE_CALL(cudaMalloc((void**)&dKernel, memSizeReal));

    // Copy host memory to device
    CUDA_SAFE_CALL(cudaMemcpy(dKernel, hKernel, memSizeReal,
                              cudaMemcpyHostToDevice));

	//allocate device memory for transforms
	Complex* dKernelTransform;
	Complex* dSignalTransform;
	CUDA_SAFE_CALL(cudaMalloc((void**)&dKernelTransform, memSizeComplex));
	CUDA_SAFE_CALL(cudaMalloc((void**)&dSignalTransform, memSizeComplex));
		
    // CUFFT plan
    cufftHandle fPlan;//forward plan
	cufftHandle rPlan;//reverse plan
	
    CUFFT_SAFE_CALL(cufftPlan1d(&fPlan, n, CUFFT_R2C, 1));
	CUFFT_SAFE_CALL(cufftPlan1d(&rPlan, n, CUFFT_C2R, 1));

    // Transform signal and kernel
	
    CUFFT_SAFE_CALL(cufftExecR2C(fPlan, (cufftReal *)dSignal, (cufftComplex *)dSignalTransform));
    CUFFT_SAFE_CALL(cufftExecR2C(fPlan, (cufftReal *)dKernel, (cufftComplex *)dKernelTransform));

    // Multiply the coefficients together and normalize the result
    ComplexPointwiseMulAndScale<<<32, 256>>>(dSignalTransform, dKernelTransform, n, 1.0f / n);

    // Check if kernel execution generated and error
    CUT_CHECK_ERROR("Kernel execution failed [ ComplexPointwiseMulAndScale ]");

    // Transform signal back
    CUFFT_SAFE_CALL(cufftExecC2R(rPlan, (cufftComplex *)dSignalTransform, (cufftReal *)dSignal));

    CUDA_SAFE_CALL(cudaMemcpy(hSignal, dSignal, memSizeReal,
                              cudaMemcpyDeviceToHost));

	
    //Destroy CUFFT context
    CUFFT_SAFE_CALL(cufftDestroy(fPlan));
	CUFFT_SAFE_CALL(cufftDestroy(rPlan));

    // cleanup memory
    free(hSignal);
    free(hKernel);
    CUDA_SAFE_CALL(cudaFree(dSignal));
    CUDA_SAFE_CALL(cudaFree(dKernel));
	CUDA_SAFE_CALL(cudaFree(dKernelTransform));
	CUDA_SAFE_CALL(cudaFree(dSignalTransform));
}

int main(int argc,char* argv)
{
	CUT_DEVICE_INIT(argc,argv);
	const int n = 256;
	float signal[n]; 	
	float kernel[n];
		
	//printf("the signal: \n");
	
	//initialize signal and filter
	for (unsigned int i = 0; i<n;i++)
	{
		if(i%2==0)
			signal[i] = 1;
		else
			signal[i] = -1;
			
			//printf("%1.1f ", signal[i]);
			
		if(i<(1+n/2))
			kernel[i]=1;
		else
			kernel[i]=0;
			
	}
	
	/*
	printf("\n \n initialize filteredSignal's memory to constant: \n");
	 
	for (unsigned int i = 0; i<n;i++)
		printf("%1.1f ", filteredSignal[i]);
		
	printf("\n \n the kernel: \n");
	
	for (unsigned int i = 0; i<n;i++)
		printf("%1.1f ", kernel[i]);
		
	printf("\n \n the filtered signal: \n");
	*/
	
	//define timing variables
	time_t startTime;
	time_t endTime;
	double runTime;
	double timePerCall;
	int nCalls = 100000;
	
	
	time(&startTime);
	//call the filtering function
	for(unsigned int i = 0; i<nCalls;i++)
	fftFilter( &(signal[0]), &(kernel[0]), n);
	
	time(&endTime);
	printf("start: %ld ", startTime);
	printf("end: %ld ", endTime);
	runTime = difftime(endTime,startTime);
	timePerCall = runTime*(double)1000/(double)nCalls;
	
	
	/*	
	for (unsigned int i = 0; i<n;i++)
		printf("%1.1f ", filteredSignal[i]);
	*/
		
	printf("\n \n total runtime for %i calls was %f seconds", nCalls, runTime);
	printf("\n \n time per call: %f ms" ,timePerCall);
	getchar();
	return 0;
}

////////////////////////////////////////////////////////////////////////////////
// Complex operations
////////////////////////////////////////////////////////////////////////////////

// Complex addition
static __device__ __host__ inline Complex ComplexAdd(Complex a, Complex b)
{
    Complex c;
    c[0] = a[0] + b[0];
    c[1] = a[1] + b[1];
    return c;
}

// Complex scale
static __device__ __host__ inline Complex ComplexScale(Complex a, float s)
{
    Complex c;
    c[0] = s * a[0];
    c[1] = s * a[1];
    return c;
}

// Complex multiplication
static __device__ __host__ inline Complex ComplexMul(Complex a, Complex b)
{
    Complex c;
    c[0] = a[0] * b[0] - a[1] * b[1];
    c[1] = a[0] * b[1] + a[1] * b[0];
    return c;
}

// Complex pointwise multiplication
static __global__ void ComplexPointwiseMulAndScale(Complex* a, const Complex* b, int size, float scale)
{
    const int numThreads = blockDim.x * gridDim.x;
    const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
    for (int i = threadID; i < size; i += numThreads)
        a[i] = ComplexScale(ComplexMul(a[i], b[i]), scale);     
} 