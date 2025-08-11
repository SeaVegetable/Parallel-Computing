#include <petscmat.h>
#include "FileManager.hpp"

int main (int argc, char *argv[])
{
    int p, q, nElemX, nElemY, part_num_1d, dim;
    double Lx, Ly;
    std::string base_name;

    std::string file_info = "info.txt";

    FileManager * fm = new FileManager();
    fm->ReadPreprocessInfo(file_info, p, q, Lx, Ly, nElemX, nElemY, part_num_1d, dim, base_name);

    PetscInitialize(&argc, &argv, NULL, NULL);
    PetscMPIInt rank, size;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    if (rank == 0)
    {
        std::cout << "p: " << p << std::endl;
        std::cout << "q: " << q << std::endl;
        std::cout << "Lx: " << Lx << std::endl;
        std::cout << "Ly: " << Ly << std::endl;
        std::cout << "nElemX: " << nElemX << std::endl;
        std::cout << "nElemY: " << nElemY << std::endl;
        std::cout << "part_num_1d: " << part_num_1d << std::endl;
        std::cout << "dim: " << dim << std::endl;
        std::cout << "base_name: " << base_name << std::endl;
    }

    int nlocalfunc;
    std::string base_name_fem = "part_fem";
    std::string filename_fem = fm->GetPartitionFilename(base_name_fem, rank);
    fm->ReadPartition(filename_fem, nlocalfunc);

    std::string filename_dr = "coordinate";

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
    MatSetType(K, MATAIJCUSPARSE);
    MatSetSizes(K, nlocalfunc, nlocalfunc, PETSC_DETERMINE, PETSC_DETERMINE);
    MatSetPreallocationCOO(K, nnz, d_rows, d_cols);
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

    delete fm;

    PetscFinalize();
    return 0;
}