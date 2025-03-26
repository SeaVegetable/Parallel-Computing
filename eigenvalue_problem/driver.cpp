#include <petscksp.h>
#include <slepcsvd.h>
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

typedef struct
{
    Mat K;
    KSP ksp;
} UserCtx;

void MyMatMult(Mat K, Vec x, Vec y)
{
    UserCtx *ctx;
    Vec t;

    MatShellGetContext(K, &ctx);
    MatMult(ctx->K, x, t);
    KSPSolve(ctx->ksp, t, y);
    VecDestroy(&t);
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

    // Compute the condition number
    KSP ksp;
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, globalassem_fem->K, globalassem_fem->K);
    KSPSetFromOptions(ksp);

    UserCtx ctx;
    Mat M;
    PetscInt mm, nn;
    MatGetSize(globalassem->K, &mm, &nn);
    MatCreateShell(PETSC_COMM_WORLD, mm, nn, mm, nn, &ctx, &M);
    MatShellSetOperation(M, MATOP_MULT, (void(*)(void))MyMatMult);

    ctx.ksp = ksp;
    ctx.K = globalassem->K;

    SVD svd;
    SVDCreate(PETSC_COMM_WORLD, &svd);
    SVDSetOperators(svd, M, NULL);
    SVDSetProblemType(svd, SVD_STANDARD);

    SVDSetDimensions(svd, 2, PETSC_DEFAULT, PETSC_DEFAULT);
    SVDSetWhichSingularTriplets(svd, SVD_LARGEST);

    PetscInt nconv;
    SVDSolve(svd);
    SVDGetConverged(svd, &nconv);

    PetscReal smax;
    if (nconv > 0)
        SVDGetSingularTriplet(svd, 0, &smax, NULL, NULL);
    else
        PetscPrintf(PETSC_COMM_WORLD, "Maximum singular value: Diverge\n");

    SVDSetWhichSingularTriplets(svd, SVD_SMALLEST);
    SVDSolve(svd);
    SVDGetConverged(svd, &nconv);

    PetscReal smin;
    if (nconv > 0)
        SVDGetSingularTriplet(svd, 0, &smin, NULL, NULL);
    else
        PetscPrintf(PETSC_COMM_WORLD, "Maximum singular value: Diverge\n");

    PetscPrintf(PETSC_COMM_WORLD, "Maximum singular value: %g\n", smax);
    PetscPrintf(PETSC_COMM_WORLD, "Minimum singular value: %g\n", smin);
    PetscPrintf(PETSC_COMM_WORLD, "Condition number: %g\n", smax/smin);
    
    delete fm; fm = nullptr;
    delete elem; elem = nullptr;
    delete locassem; locassem = nullptr;
    delete globalassem; globalassem = nullptr;
    delete elem_fem; elem_fem = nullptr;
    delete locassem_fem; locassem_fem = nullptr;
    delete globalassem_fem; globalassem_fem = nullptr;
    
    PetscFinalize();
    return 0;
}