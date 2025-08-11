#include <cuda_runtime.h>
#include "memory.hpp"

void MallocDeviceMemory(void **ptr, size_t size)
{
    cudaError_t err = cudaMalloc(ptr, size);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error allocating device memory: %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

void FreeDeviceMemory(void *ptr)
{
    cudaError_t err = cudaFree(ptr);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error freeing device memory: %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

void CopyToDevice(void *dst, const void *src, size_t size)
{
    cudaError_t err = cudaMemcpy(dst, src, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error copying to device memory: %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}