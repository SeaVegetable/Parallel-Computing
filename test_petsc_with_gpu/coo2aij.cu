#include <petscmat.h>
#include "FileManager.hpp"

int main (int argc, char *argv[])
{
    PetscInitialize(&argc, &argv, NULL, NULL);

    std::string filename_dr = "dryrun";

    std::vector<int> rows{};
    std::vector<int> cols{};
    int nnz = 0;
    for (int ii = 0; ii < part_num_1d * part_num_1d; ++ii)
    {
        std::vector<int> rows_temp{};
        std::vector<int> cols_temp{};
        std::string filename = fm->GetNonZeroCoordinateFilename(filename_dr, ii);
        fm->ReadNonZeroCoordinate(filename, nnz, rows_temp, cols_temp);
        rows.insert(rows.end(), rows_temp.begin(), rows_temp.end());
        cols.insert(cols.end(), cols_temp.begin(), cols_temp.end());
    }

    Mat K;
    PetscInt * d_rows, * d_cols;
    cudaMalloc((void**)&d_rows, nnz * sizeof(PetscInt));
    cudaMalloc((void**)&d_cols, nnz * sizeof(PetscInt));
    cudaMemcpy(d_rows, rows.data(), nnz * sizeof(PetscInt), cudaMemcpyHostToDevice);
    cudaMemcpy(d_cols, cols.data(), nnz * sizeof(PetscInt), cudaMemcpyHostToDevice);
    MatCreate(PETSC_COMM_WORLD, &K);
    MatSetSizes(K, nlocalfunc, nlocalfunc, PETSC_DETERMINE, PETSC_DETERMINE);
    MatSetType(K, MATAIJCUSPARSE);
    MatSetPreallocationCOO(A, nnz, d_rows, d_cols);
    cudaFree(d_rows);
    cudaFree(d_cols);

    PetscScalar * d_val;
    cudaMalloc((void**)&d_val, nnz * sizeof(PetscScalar));
    cudaMemset(d_val, 1.0, nnz * sizeof(PetscScalar));
    MatSetValuesCOO(K, d_val, INSERT_VALUES);
    cudaFree(d_val);

    MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);

    PetscPrintf(PETSC_COMM_WORLD, "Matrix K assembled with %d non-zero entries.\n", nnz);

    MatView(K, PETSC_VIEWER_STDOUT_WORLD);

    PetscFinalize();
    return 0;
}