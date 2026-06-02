#include <iostream>
#include <petscvec.h>
#include <cuda_runtime.h>

int main(int argc, char **argv)
{
    PetscInitialize(&argc, &argv, NULL, NULL);

    int ndevices;
    cudaGetDeviceCount(&ndevices);
    std::cout << "Number of CUDA devices: " << ndevices << std::endl;

    int rank, size;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    int device_id = rank % ndevices;
    cudaSetDevice(device_id);

    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, device_id);

    std::cout << "MPI rank " << rank << " assigned to GPU " << device_id << ": " << prop.name << std::endl;

    PetscFinalize();
    return 0;
}