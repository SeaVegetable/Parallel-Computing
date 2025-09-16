#include <petscvec.h>
#include <cuda_runtime.h>

int main(int argc, char **argv)
{
    PetscInitialize(&argc, &argv, NULL, NULL);

    int ndevices;
    cudaGetDeviceCount(&ndevices);
    PetscPrintf(PETSC_COMM_WORLD, "Number of CUDA devices: %d\n", ndevices);

    int rank, size;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    int device_id = rank % ndevices;
    cudaSetDevice(device_id);

    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, device_id);

    PetscPrintf(PETSC_COMM_WORLD,
                "MPI rank %d assigned to GPU %d: %s\n",
                rank, device_id, prop.name);

    PetscFinalize();
    return 0;
}