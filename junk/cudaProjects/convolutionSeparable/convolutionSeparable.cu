/*
 * Copyright 1993-2007 NVIDIA Corporation.  All rights reserved.
 *
 * NOTICE TO USER:   
 *
 * This source code is subject to NVIDIA ownership rights under U.S. and 
 * international Copyright laws.  Users and possessors of this source code 
 * are hereby granted a nonexclusive, royalty-free license to use this code 
 * in individual and commercial software.
 *
 * NVIDIA MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THIS SOURCE 
 * CODE FOR ANY PURPOSE.  IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR 
 * IMPLIED WARRANTY OF ANY KIND.  NVIDIA DISCLAIMS ALL WARRANTIES WITH 
 * REGARD TO THIS SOURCE CODE, INCLUDING ALL IMPLIED WARRANTIES OF 
 * MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.
 * IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL, 
 * OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS 
 * OF USE, DATA OR PROFITS,  WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE 
 * OR OTHER TORTIOUS ACTION,  ARISING OUT OF OR IN CONNECTION WITH THE USE 
 * OR PERFORMANCE OF THIS SOURCE CODE.  
 *
 * U.S. Government End Users.   This source code is a "commercial item" as 
 * that term is defined at  48 C.F.R. 2.101 (OCT 1995), consisting  of 
 * "commercial computer  software"  and "commercial computer software 
 * documentation" as such terms are  used in 48 C.F.R. 12.212 (SEPT 1995) 
 * and is provided to the U.S. Government only as a commercial end item.  
 * Consistent with 48 C.F.R.12.212 and 48 C.F.R. 227.7202-1 through 
 * 227.7202-4 (JUNE 1995), all U.S. Government End Users acquire the 
 * source code with only those rights set forth herein. 
 *
 * Any use of this source code in individual and commercial software must 
 * include, in the user documentation and internal comments to the code,
 * the above Disclaimer and U.S. Government End Users Notice.
 */

/*
 * This sample implements a separable convolution filter 
 * of a 2D signal with a gaussian kernel.
 */



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cutil.h>



////////////////////////////////////////////////////////////////////////////////
// Common host and device functions
////////////////////////////////////////////////////////////////////////////////
//Round a / b to nearest higher integer value
int iDivUp(int a, int b){
    return (a % b != 0) ? (a / b + 1) : (a / b);
}

//Round a / b to nearest lower integer value
int iDivDown(int a, int b){
    return a / b;
}

//Align a to nearest higher multiple of b
int iAlignUp(int a, int b){
    return (a % b != 0) ?  (a - a % b + b) : a;
}

//Align a to nearest lower multiple of b
int iAlignDown(int a, int b){
    return a - a % b;
}


////////////////////////////////////////////////////////////////////////////////
// GPU convolution
////////////////////////////////////////////////////////////////////////////////
//Global macro, controlling innermost convolution loop unrolling
#define UNROLL_INNER
#include <convolutionSeparable_kernel.cu>



////////////////////////////////////////////////////////////////////////////////
// Data configuration
//////////////////////////////////------//////////////////////////////////////////////
//Image width should be aligned to maximum coalesced read/write size
//for best global memory performance in both row and column filter.
#ifdef __DEVICE_EMULATION__
//Reduce problem size to have reasonable emulation time
const int      DATA_W = iAlignUp(256, 16);
const int      DATA_H = 256;
#else
const int      DATA_W = iAlignUp(8192, 16);
const int      DATA_H = 1024;//CHANGE NUMBER OF ROWS HERE!
#endif
const int   DATA_SIZE = DATA_W * DATA_H * sizeof(float);
const int KERNEL_SIZE = KERNEL_W * sizeof(float);

//////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// filtering on the gpu wrapped for labview, matlab, etc /////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __declspec(dllexport) void gpuFilter(float* h_Result,  int dataW, int dataH)
//void gpuFilter(float* h_Result,   float* h_Data,    int dataW, int dataH)
	{
	////////////////////////// kernel init /////////////////////////////////////
	 //declare and allocate for kernel on host
	 float *h_Kernel, kernelSum = 0;
	 h_Kernel    = (float *)malloc(KERNEL_SIZE);
	 
	 //build the kernel
        for(unsigned int i = 0; i < KERNEL_W; i++){
            float dist = (float)(i - KERNEL_RADIUS) / (float)KERNEL_RADIUS;
            h_Kernel[i] = expf(- dist * dist / 2);
            kernelSum += h_Kernel[i];
        }
        for(unsigned int i = 0; i < KERNEL_W; i++)
            h_Kernel[i] /= kernelSum;

	//copy host kernel to device
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(d_Kernel, h_Kernel, KERNEL_SIZE) );
	free(h_Kernel);//free the host kernel memory
	////////////////////////// end kernel init /////////////////////////////////////
	
	////////////////////////////////fitlering procedure/////////////////////////////
	//declare and initialize device variables
	float  *d_Data, *d_Result, *d_Temp;
	
	CUDA_SAFE_CALL( cudaMalloc( (void **)&d_Data, DATA_SIZE) );
    CUDA_SAFE_CALL( cudaMalloc( (void **)&d_Result, DATA_SIZE) );
    CUDA_SAFE_CALL( cudaMalloc( (void **)&d_Temp , DATA_SIZE) );
	
	//transfer the data
	CUDA_SAFE_CALL( cudaMemcpy(d_Data, h_Result, DATA_SIZE, cudaMemcpyHostToDevice) );
    
	//setup for the call to the device
	dim3 blockGridRows(iDivUp(DATA_W, ROW_TILE_W), DATA_H);
    dim3 threadBlockRows(KERNEL_RADIUS_ALIGNED + ROW_TILE_W + KERNEL_RADIUS);
	
	//call the device funtion
	convolutionRowGPU<<<blockGridRows, threadBlockRows>>>( d_Result, d_Data,  dataW, dataH);
	
	//obtain the result
	CUDA_SAFE_CALL( cudaMemcpy(h_Result, d_Result, DATA_SIZE, cudaMemcpyDeviceToHost) );
	
	//clean up device memory
	
	CUDA_SAFE_CALL( cudaFree(d_Data) );
    CUDA_SAFE_CALL( cudaFree(d_Result) );
	CUDA_SAFE_CALL( cudaFree(d_Temp) );
	//////////////////////////// end filtering procedure ////////////////////////////////
	}
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

int __stdcall DllMain(void)
{
return 0;
}
