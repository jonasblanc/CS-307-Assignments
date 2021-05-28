/*
============================================================================
Filename    : implementation.cu
Author      : Jonas Blanc, Mélissa Gehring
SCIPER      : 287508, 264265
============================================================================
*/

#include <iostream>
#include <iomanip>
#include <sys/time.h>
#include <cuda_runtime.h>
using namespace std;

#define INPUT(I,J) input[(I)*length+(J)]
#define S_DATA(I,J) sdata[(I)*s_length+(J)]

// CPU Baseline
void array_process(double *input, double *output, int length, int iterations)
{
    double *temp;

    for(int n=0; n<(int) iterations; n++)
    {
        for(int i=1; i<length-1; i++)
        {
            for(int j=1; j<length-1; j++)
            {
                output[(i)*(length)+(j)] = (input[(i-1)*(length)+(j-1)] +
                                            input[(i-1)*(length)+(j)]   +
                                            input[(i-1)*(length)+(j+1)] +
                                            input[(i)*(length)+(j-1)]   +
                                            input[(i)*(length)+(j)]     +
                                            input[(i)*(length)+(j+1)]   +
                                            input[(i+1)*(length)+(j-1)] +
                                            input[(i+1)*(length)+(j)]   +
                                            input[(i+1)*(length)+(j+1)] ) / 9;

            }
        }
        output[(length/2-1)*length+(length/2-1)] = 1000;
        output[(length/2)*length+(length/2-1)]   = 1000;
        output[(length/2-1)*length+(length/2)]   = 1000;
        output[(length/2)*length+(length/2)]     = 1000;

        temp = input;
        input = output;
        output = temp;
    }
}

__global__ void compute_one_iteration(double *input, double *output, int length) {
	int j = (blockIdx.x * blockDim.x) + threadIdx.x ; 
    int i = (blockIdx.y * blockDim.y) + threadIdx.y ;
    int array_index = (i * length) + j;

	if(0 < i && i < length-1 && 0 < j && j < length - 1){
		output[array_index] = (input[(i-1)*(length)+(j-1)] +
                                input[(i-1)*(length)+(j)]   +
                                input[(i-1)*(length)+(j+1)] +
                                input[(i)*(length)+(j-1)]   +
                                input[(i)*(length)+(j)]     +
                                input[(i)*(length)+(j+1)]   +
                                input[(i+1)*(length)+(j-1)] +
                                input[(i+1)*(length)+(j)]   +
                                input[(i+1)*(length)+(j+1)] ) / 9;
	}
	output[(length/2-1)*length+(length/2-1)] = 1000;
    output[(length/2)*length+(length/2-1)]   = 1000;
    output[(length/2-1)*length+(length/2)]   = 1000;
    output[(length/2)*length+(length/2)]     = 1000;
 }

 __global__ void compute_one_iteration_smart_mid(double *input, double *output, int length) {

    int j = (blockIdx.x * blockDim.x) + threadIdx.x ; 
    int i = (blockIdx.y * blockDim.y) + threadIdx.y ;
    int array_index = (i * length) + j;

    int mid1 = (length/2-1)*length+(length/2-1);
    int mid2 = (length/2)*length+(length/2-1);
    int mid3 = (length/2-1)*length+(length/2);
    int mid4 = (length/2)*length+(length/2);

    if(array_index == mid1 || array_index == mid2 || array_index == mid3 || array_index == mid4){
        return;
    }

	if(0 < i && i < length-1 && 0 < j && j < length - 1){
		output[array_index] = (input[(i-1)*(length)+(j-1)] +
                                input[(i-1)*(length)+(j)]   +
                                input[(i-1)*(length)+(j+1)] +
                                input[(i)*(length)+(j-1)]   +
                                input[(i)*(length)+(j)]     +
                                input[(i)*(length)+(j+1)]   +
                                input[(i+1)*(length)+(j-1)] +
                                input[(i+1)*(length)+(j)]   +
                                input[(i+1)*(length)+(j+1)] ) / 9;
    }
 }

 __global__ void compute_one_iteration_shared(double *input, double *output, int length) {
    extern __shared__ double sdata[]; // Used in macro S_DATA
    
    int j = (blockIdx.x * (blockDim.x - 2)) + threadIdx.x ; 
    int i = (blockIdx.y * (blockDim.y - 2)) + threadIdx.y ;
    int array_index = (i * length) + j;

    int s_i = threadIdx.y;
    int s_j = threadIdx.x;
    int s_length = blockDim.x; // Used in macro S_DATA

    // Load shared memory
    if(0 <= i && i <= length-1 && 0 <= j && j <= length - 1){
        S_DATA(s_i, s_j) = INPUT(i,j);
        __syncthreads();
    }

    if(0 < s_i && s_i < s_length-1 && 0 < s_j && s_j < s_length - 1){
        if(0 < i && i < length-1 && 0 < j && j < length - 1){
            output[array_index] = ( S_DATA(s_i - 1, s_j -1 )    +
                                    S_DATA(s_i - 1, s_j)        +
                                    S_DATA(s_i - 1, s_j + 1)    +
                                    S_DATA(s_i, s_j - 1)        +
                                    S_DATA(s_i, s_j)            +
                                    S_DATA(s_i, s_j + 1)        +
                                    S_DATA(s_i + 1, s_j - 1)    +
                                    S_DATA(s_i + 1, s_j)        +
                                    S_DATA(s_i + 1, s_j + 1)    ) / 9;
        }
    }

	output[(length/2-1)*length+(length/2-1)] = 1000;
    output[(length/2)*length+(length/2-1)]   = 1000;
    output[(length/2-1)*length+(length/2)]   = 1000;
    output[(length/2)*length+(length/2)]     = 1000;
 }


// GPU Optimized function
void GPU_array_process(double *input, double *output, int length, int iterations)
{
    //Cuda events for calculating elapsed time
    cudaEvent_t cpy_H2D_start, cpy_H2D_end, comp_start, comp_end, cpy_D2H_start, cpy_D2H_end;
    cudaEventCreate(&cpy_H2D_start);
    cudaEventCreate(&cpy_H2D_end);
    cudaEventCreate(&cpy_D2H_start);
    cudaEventCreate(&cpy_D2H_end);
    cudaEventCreate(&comp_start);
    cudaEventCreate(&comp_end);

    /* Preprocessing */
    size_t array_size = length * length * sizeof(double);
    double* gpu_array_in;
    cudaMalloc( (void**)&gpu_array_in, array_size);
    double* gpu_array_out;
    cudaMalloc( (void**)&gpu_array_out, array_size);

    cudaEventRecord(cpy_H2D_start);
    /* Copying array from host to device */
    cudaMemcpy((void*)gpu_array_in, (void*)input, array_size, cudaMemcpyHostToDevice);
    cudaMemcpy((void*)gpu_array_out, (void*)output, array_size, cudaMemcpyHostToDevice);

    cudaEventRecord(cpy_H2D_end);
    cudaEventSynchronize(cpy_H2D_end);

    cudaEventRecord(comp_start);
    /* GPU calculation */
    // SQUARE THREADS BLOCKS
    size_t threadBlockSide = 8;
    size_t nbBlockSide = length / threadBlockSide;

    // If not a multiple
    if(length % threadBlockSide != 0){
        nbBlockSide += 1;
    }

    dim3 thrsPerBlock(threadBlockSide,threadBlockSide); 
    dim3 nBlks(nbBlockSide,nbBlockSide);   
    
    // SHARED MEMORY
    size_t threadBlockSide_shared = 32;
    size_t nbBlockSide_shared = length / (threadBlockSide_shared -2);

    // If not a multiple
    if(length % (threadBlockSide_shared -2) != 0){
        nbBlockSide_shared += 1;
    }

    size_t smemSize_shared = threadBlockSide_shared * threadBlockSide_shared * sizeof(double);

    dim3 thrsPerBlock_shared(threadBlockSide_shared, threadBlockSide_shared); 
    dim3 nBlks_shared(nbBlockSide_shared, nbBlockSide_shared); 

    // ROW THREADS BLOCKS
    size_t threadBlockSide_row = length;
    size_t nbBlockSide_row = 1;

    if(threadBlockSide_row > 1024){
        threadBlockSide_row = 512;
        nbBlockSide_row = length / threadBlockSide_row;

        // If not a multiple
        if(length % threadBlockSide_row != 0){
            nbBlockSide_row += 1;
        }
    }

    dim3 thrsPerBlock_row(threadBlockSide_row, 1); 
    dim3 nBlks_row(nbBlockSide_row, length); 
    
    double *temp;
    for(int n = 0; n < iterations; n++)
    {
        //compute_one_iteration <<< nBlks, thrsPerBlock >>> (gpu_array_in, gpu_array_out, length);
        compute_one_iteration_smart_mid <<< nBlks, thrsPerBlock >>> (gpu_array_in, gpu_array_out, length);
        //compute_one_iteration <<< nBlks_row, thrsPerBlock_row >>> (gpu_array_in, gpu_array_out, length);
        //compute_one_iteration_shared <<< nBlks_shared, thrsPerBlock_shared, smemSize_shared >>> (gpu_array_in, gpu_array_out, length);

        temp = gpu_array_in;
        gpu_array_in = gpu_array_out;
        gpu_array_out = temp;
    }
    cudaEventRecord(comp_end);
    cudaEventSynchronize(comp_end);

    cudaEventRecord(cpy_D2H_start);
    
    /* Copying array from device to host goes here */
    cudaMemcpy((void*)output, (void*)gpu_array_in, array_size, cudaMemcpyDeviceToHost);

    cudaEventRecord(cpy_D2H_end);
    cudaEventSynchronize(cpy_D2H_end);

    /* Postprocessing goes here */
    cudaFree((void*)gpu_array_in);
    cudaFree((void*)gpu_array_out);
    
    float time;
    cudaEventElapsedTime(&time, cpy_H2D_start, cpy_H2D_end);
    cout<<"Host to Device MemCpy takes "<<setprecision(4)<<time/1000<<"s"<<endl;

    cudaEventElapsedTime(&time, comp_start, comp_end);
    cout<<"Computation takes "<<setprecision(4)<<time/1000<<"s"<<endl;

    cudaEventElapsedTime(&time, cpy_D2H_start, cpy_D2H_end);
    cout<<"Device to Host MemCpy takes "<<setprecision(4)<<time/1000<<"s"<<endl;
}
