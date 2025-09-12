#include <petscvec.h>
#include <cuda_runtime.h>

int main(int argc, char **argv)
{
    PetscInitialize(&argc, &argv, NULL, NULL);

    int ndevices;
    cudaGetDeviceCount(&ndevices);
    PetscPrintf(PETSC_COMM_WORLD, "Number of CUDA devices: %d\n", ndevices);

    PetscFinalize();
    return 0;
}