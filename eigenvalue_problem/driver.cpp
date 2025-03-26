#include <petscksp.h>
#include "FileManager.hpp"
#include "GlobalAssembly.hpp"

double GetMaximumEigenvalue(Mat &Kiga, Mat &Kfem, Vec &u0,
    const double &tol)
{   
    double norm_u0;
    VecNorm(u0, NORM_2, &norm_u0);
    
    VecScale(u0, 1.0 / norm_u0);

    Vec u;
    VecDuplicate(u0, &u);

    double lambda_old = 0;
    double lambda_new;

    KSP ksp;
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, Kfem, Kfem);
    KSPSetFromOptions(ksp);
    KSPSetTolerances(ksp, 1e-10, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);

    while (true)
    {
        MatMult(Kiga, u0, u);

        KSPSolve(ksp, u, u);

        KSPSolve(ksp, u, u);

        MatMult(Kiga, u, u0);

        VecNorm(u0, NORM_2, &norm_u0);
        VecScale(u0, 1.0 / norm_u0);

        lambda_new = norm_u0;

        if (fabs(lambda_new - lambda_old) < tol)
        {
            break;
        }

        lambda_old = lambda_new;
    }

    VecDestroy(&u);
    KSPDestroy(&ksp);

    return lambda_new;
}

double GetMaximumEigenvalue(Mat &K, Vec &u0, const double &tol)
{   
    double norm_u0;
    VecNorm(u0, NORM_2, &norm_u0);
    
    VecScale(u0, 1.0 / norm_u0);

    Vec u;
    VecDuplicate(u0, &u);

    double lambda_old = 0;
    double lambda_new;

    while (true)
    {
        MatMult(K, u0, u);

        MatMult(K, u, u0);

        VecNorm(u0, NORM_2, &norm_u0);
        VecScale(u0, 1.0 / norm_u0);

        lambda_new = norm_u0;

        if (fabs(lambda_new - lambda_old) < tol)
        {
            break;
        }

        lambda_old = lambda_new;
    }

    VecDestroy(&u);

    return lambda_new;
}

double GetMinimumEigenvalue(Mat &K, Vec &u0, const double &tol)
{   
    double norm_u0;
    VecNorm(u0, NORM_2, &norm_u0);
    
    VecScale(u0, 1.0 / norm_u0);

    Vec u;
    VecDuplicate(u0, &u);

    double lambda_old = 0;
    double lambda_new;

    KSP ksp;
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, K, K);
    KSPSetFromOptions(ksp);
    KSPSetTolerances(ksp, 1e-10, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);

    while (true)
    {
        KSPSolve(ksp, u0, u);

        KSPSolve(ksp, u, u0);

        VecNorm(u0, NORM_2, &norm_u0);
        VecScale(u0, 1.0 / norm_u0);

        lambda_new = norm_u0;

        if (fabs(lambda_new - lambda_old) < tol)
        {
            break;
        }

        lambda_old = lambda_new;
    }

    VecDestroy(&u);
    KSPDestroy(&ksp);

    return lambda_new;
}

int main(int argc, char *argv[])
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

    std::vector<double> CP;
    std::vector<int> ID;
    std::vector<int> IEN;
    std::vector<double> elem_size1;
    std::vector<double> elem_size2;
    std::vector<double> NURBSExtraction1;
    std::vector<double> NURBSExtraction2;
    int nlocalfunc;
    int nlocalelemx;
    int nlocalelemy;

    std::string filename = fm->GetPartitionFilename(base_name, rank);
    fm->ReadPartition(filename, nlocalfunc,
        nlocalelemx, nlocalelemy,
        elem_size1, elem_size2,
        CP, ID, IEN, NURBSExtraction1, NURBSExtraction2);
    
    Element * elem = new Element(p, q);
    const int nLocBas = elem->GetNumLocalBasis();
    LocalAssembly * locassem = new LocalAssembly(p, q);
    GlobalAssembly * globalassem = new GlobalAssembly(IEN, ID, locassem,
        nLocBas, nlocalfunc, nlocalelemx, nlocalelemy);
    
    MatSetOption(globalassem->K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

    globalassem->AssemStiffnessLoad(locassem, IEN, ID, CP,
        NURBSExtraction1, NURBSExtraction2,
        elem_size1, elem_size2, elem);

    MPI_Barrier(PETSC_COMM_WORLD);

    std::vector<double> CP_fem;
    std::vector<int> ID_fem;
    std::vector<int> IEN_fem;
    int nlocalfunc_fem;
    int nlocalelemx_fem;
    int nlocalelemy_fem;

    std::string base_name_fem = "part_fem";
    std::string filename_fem = fm->GetPartitionFilename(base_name_fem, rank);
    fm->ReadPartition(filename_fem, nlocalfunc_fem,
        nlocalelemx_fem, nlocalelemy_fem,
        CP_fem, ID_fem, IEN_fem);

    ElementFEM * elem_fem = new ElementFEM(1, 1);
    LocalAssembly * locassem_fem = new LocalAssembly(1, 1);
    GlobalAssembly * globalassem_fem = new GlobalAssembly(IEN_fem, ID_fem, locassem_fem,
        4, nlocalfunc_fem, nlocalelemx_fem, nlocalelemy_fem);
    
    MatSetOption(globalassem_fem->K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    
    globalassem_fem->AssemStiffnessLoad(locassem_fem, IEN_fem, ID_fem, CP_fem, elem_fem);

    MPI_Barrier(PETSC_COMM_WORLD);

    Vec x;
    VecDuplicate(globalassem->F, &x);
    VecSet(x, 1.0);

    const double tol = 1.0e-8;

    const double max_eigen = sqrt(
        GetMaximumEigenvalue(globalassem->K, globalassem_fem->K, x, tol));
    const double min_eigen = 1.0/sqrt(
        GetMaximumEigenvalue(globalassem_fem->K, globalassem->K, x, tol));
    
    const double cond = max_eigen / min_eigen;

    const double max_eigen_iga = sqrt(
        GetMaximumEigenvalue(globalassem->K, x, tol));
    const double min_eigen_iga = 1.0/sqrt(
        GetMinimumEigenvalue(globalassem->K, x, tol));
    
    const double cond_iga = max_eigen_iga / min_eigen_iga;

    if (rank == 0)
    {
        std::cout << std::setprecision(15);
        std::cout << "cond_iga: " << cond_iga << std::endl;
        std::cout << "cond: " << cond << std::endl;
    }

    VecDestroy(&x);
    
    delete fm; fm = nullptr;
    delete elem; elem = nullptr;
    delete locassem; locassem = nullptr;
    delete globalassem; globalassem = nullptr;
    delete elem_fem; elem_fem = nullptr;
    delete locassem_fem; locassem_fem = nullptr;
    delete globalassem_fem; globalassem_fem = nullptr;
    return 0;
}