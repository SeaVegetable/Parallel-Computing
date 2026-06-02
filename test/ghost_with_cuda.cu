#include <iostream>
#include <petscvec.h>
#include <cuda_runtime.h>

__global__ void setLocalKernel(PetscScalar *xarr, PetscInt nlocal, PetscInt nghost, PetscInt rank)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < nlocal) {
    xarr[idx+nghost] = rank * nlocal + idx;
  }
}

int get_local_rank()
{
  if (const char* v = getenv("SLURM_LOCALID")) return atoi(v);
  if (const char* v = getenv("OMPI_COMM_WORLD_LOCAL_RANK")) return atoi(v);
  if (const char* v = getenv("MPI_LOCALRANKID")) return atoi(v);
  return 0;
}

int main(int argc, char **argv) {
  int local_rank = get_local_rank();
  int ndevices;
  cudaGetDeviceCount(&ndevices);
  cudaSetDevice(local_rank % ndevices);

  PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));

  PetscInt rank, size;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);

  PetscInt nlocal = 10;
  PetscInt nghost = 1;
  PetscInt ghost[nghost];
  if (rank == 0)
  {
    ghost[0] = 10;
  }
  else
  {
    ghost[0] = 9;
  }
  Vec x;
  VecCreateGhost(PETSC_COMM_WORLD, nlocal, PETSC_DECIDE, nghost, ghost, &x);
  VecSetType(x, VECMPICUDA);
  VecSet(x, 0.0);
  
  Vec xlocal;
  PetscScalar *xarray;
  VecGhostGetLocalForm(x, &xlocal);
  VecCUDAGetArray(xlocal, &xarray);
  setLocalKernel<<<(nlocal + 255) / 256, 256>>>(xarray, nlocal, rank);
  cudaDeviceSynchronize();
  VecCUDARestoreArray(xlocal, &xarray);
  VecGhostRestoreLocalForm(x, &xlocal);
  VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);

  const PetscScalar *harray;
  VecGhostGetLocalForm(x, &xlocal);
  VecGetArrayRead(xlocal, &harray);
  if (rank == 0)
  {
    std::cout << "xlocal (rank 0): ";
    for (PetscInt i = 0; i < nlocal + nghost; ++i) {
      std::cout << harray[i] << " ";
    }
    std::cout << std::endl;
  }
  VecRestoreArrayRead(xlocal, &harray);
  VecGhostRestoreLocalForm(x, &xlocal);

  VecDestroy(&x);

  PetscFinalize();
  return 0;
}
