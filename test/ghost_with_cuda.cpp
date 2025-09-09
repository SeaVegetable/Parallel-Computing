#include <petscvec.h>

int main(int argc, char **argv) {
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
  else (rank == 1)
  {
    ghost[0] = 9;
  }
  Vec x;
  VecCreateGhost(PETSC_COMM_WORLD, nlocal, PETSC_DECIDE, nghost, ghost, &x);
  VecSetType(x, VECCUDA);
  
  PetscScalar *xarray;
  VecCUDAGetArray(x, &xarray);
  for (PetscInt i = 0; i < nlocal; ++i) {
    xarray[i] = rank * nlocal + i;
  }
  VecCUDARestoreArray(x, &xarray);
  
  Vec y;
  VecDuplicate(x, &y);
  VecSet(y, 0.0);

  Vec xlocal, ylocal;
  PetscScalar *xarr, *yarr;
  VecGhostGetLocalVector(x, &xlocal);
  VecGhostGetLocalVector(y, &ylocal);
  VecCUDAGetArray(xlocal, &xarr);
  VecCUDAGetArray(ylocal, &yarr);
  for (PetscInt i = nghost; i < nlocal + nghost; ++i) {
    yarr[i] = 2.0 * xarr[i];
  }
  VecCUDARestoreArray(xlocal, &xarr);
  VecCUDARestoreArray(ylocal, &yarr);
  VecGhostRestoreLocalVector(x, &xlocal);
  VecGhostRestoreLocalVector(y, &ylocal);

  VecView(x, PETSC_VIEWER_STDOUT_WORLD);
  VecView(y, PETSC_VIEWER_STDOUT_WORLD);

  VecDestroy(&x);
  VecDestroy(&y);

  PetscFinalize();
  return 0;
}