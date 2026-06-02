#include <petscksp.h>
#include <petscpc.h>
#include "FileManager.hpp"
#include "GlobalAssemblyMF.cuh"
#include "GlobalAssembly.cuh"
#include "ElementMF.hpp"
#include "QuadraturePoint.hpp"
#include "BernsteinBasis.hpp"
#include "Elem2COOGenerator.hpp"
#include "AbscissaeGenerator.hpp"
#include "IENGenerator.hpp"
#include "IDGenerator.hpp"

typedef struct {
    KSP innerksp;
} MyPCCtx;

PetscErrorCode MyPCApply(PC pc, Vec x, Vec y)
{
    MyPCCtx *ctx;
    PCShellGetContext(pc, (void**)&ctx);
    KSPSolve(ctx->innerksp, x, y);
    return 0;
}

PetscErrorCode MyPCDestroy(PC pc)
{
    MyPCCtx *ctx;
    PCShellGetContext(pc, (void**)&ctx);
    KSPDestroy(&ctx->innerksp);
    PetscFree(ctx);
    return 0;
}

typedef struct {
    int nlocalfunc;
    std::vector<double> CP;
    std::vector<int> ID;
    std::vector<int> Dir;
    std::vector<int> IEN;
    std::vector<double> elem_size1;
    std::vector<double> elem_size2;
    std::vector<double> NURBSExtraction1;
    std::vector<double> NURBSExtraction2;
    GlobalAssemblyMF * globalassem;
    QuadraturePoint * quad1;
    QuadraturePoint * quad2;
    BernsteinBasis * bernstein;
    ElementMF * elem;
} MyMeshData;

PetscErrorCode MyMatMult(Mat A, Vec x, Vec y)
{
    MyMeshData *data;
    MatShellGetContext(A, (void**)&data);

    VecSet(y, 0.0);
    
    data->globalassem->MatMulMF(data->quad1, data->quad2,
        data->IEN, data->ID, data->Dir, data->CP,
        data->NURBSExtraction1, data->NURBSExtraction2,
        data->elem_size1, data->elem_size2,
        data->elem, data->bernstein, x, y);

    return 0;
}

PetscErrorCode MyMatCreateVecs(Mat A, Vec *x, Vec *y)
{
    MyMeshData *data;
    MatShellGetContext(A, (void**)&data);

    VecCreate(PETSC_COMM_WORLD, x);
    VecSetSizes(*x, data->nlocalfunc, PETSC_DETERMINE);
    VecSetType(*x, VECCUDA);
    if(y) VecDuplicate(*x, y);
    return 0;
}

int main(int argc, char **argv)
{
    int p, q, nElemX, nElemY, part_num_1d, dim;
    double Lx, Ly;
    std::string base_name;

    std::string file_info = "fem_info.txt";

    FileManager * fm = new FileManager();
    fm->ReadPreprocessInfo(file_info, p, q, Lx, Ly, nElemX, nElemY, part_num_1d, dim, base_name);

    PetscInitialize(&argc, &argv, NULL, NULL);

    PetscInt rank, size;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);


    if(rank == 0)
    {
        std::cout << "Preprocessing information of FEM dry run:" << std::endl;
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
    std::vector<int> localID;
    std::vector<int> ghostID;
    std::vector<int> Dir;
    std::vector<int> IEN;
    std::vector<double> elem_size1;
    std::vector<double> elem_size2;
    std::vector<double> NURBSExtraction1;
    std::vector<double> NURBSExtraction2;
    int nlocalfunc;
    int nlocalelemx;
    int nlocalelemy;

    std::string filename = fm->GetPartitionFilename(base_name, rank);
    fm->ReadPartition(filename, nlocalfunc, nlocalelemx, nlocalelemy,
        elem_size1, elem_size2,
        CP, ID, localID, ghostID, Dir, IEN,
        NURBSExtraction1, NURBSExtraction2);

    std::string filename_dr = "coordinate";
    std::vector<int> rows{};
    std::vector<int> cols{};
    int nnz = 0;
    filename = fm->GetNonZeroCoordinateFilename(filename_dr, rank);

    fm->ReadNonZeroCoordinate(filename, nnz, rows, cols);

    Elem2COOGenerator * elem2coogen = new Elem2COOGenerator(4, nnz,
        nFuncX, nFuncY);

    std::vector<int> elem2coo{};
    std::vector<int> dir2coo{};
    elem2coogen->GenerateElem2COO(IEN_fem, ID_fem, rows, cols, elem2coo);
    elem2coogen->GenerateDir2COO(ID_fem, rows, dir2coo);

    PetscPrintf(PETSC_COMM_WORLD, "Finished generating elem2coo and dir2coo.\n");

    GlobalAssembly * assembly_fem = new GlobalAssembly(4,
        nnz, nlocalfunc, nElemX, nElemY, old_rows, old_cols);
    QuadraturePoint * quad1_fem = new QuadraturePoint(2, 0, 1);
    QuadraturePoint * quad2_fem = new QuadraturePoint(2, 0, 1);
    ElementFEM * elemfem = new ElementFEM(1, 1);

    assembly_fem->AssemStiffness(quad1_fem, quad2_fem,
        IEN_fem, dir2coo, CP_fem, elem2coo, elemfem);
    PetscPrintf(PETSC_COMM_WORLD, "Assembling FEM stiffness matrix...done\n");

    file_info = "info.txt";

    fm->ReadPreprocessInfo(file_info, p, q, Lx, Ly, nElemX, nElemY, part_num_1d, dim, base_name);

    std::cout << "Preprocessing information of matrix-free:" << std::endl;
    std::cout << "p: " << p << std::endl;
    std::cout << "q: " << q << std::endl;
    std::cout << "Lx: " << Lx << std::endl;
    std::cout << "Ly: " << Ly << std::endl;
    std::cout << "nElemX: " << nElemX << std::endl;
    std::cout << "nElemY: " << nElemY << std::endl;
    std::cout << "part_num_1d: " << part_num_1d << std::endl;
    std::cout << "dim: " << dim << std::endl;
    std::cout << "base_name: " << base_name << std::endl;

    int nlocalfunc;
    int nlocalelemx;
    int nlocalelemy;
    MyMeshData * data;
    PetscNew(&data);

    std::string filename = fm->GetPartitionFilename(base_name, 0);
    fm->ReadPartition(filename, nlocalfunc,
        nlocalelemx, nlocalelemy,
        data->elem_size1, data->elem_size2,
        data->CP, data->ID, data->Dir, data->IEN,
        data->NURBSExtraction1, data->NURBSExtraction2);

    data->nlocalfunc = nlocalfunc;

    ElementMF * elemmf = new ElementMF(p, q);
    int nLocBas = elemmf->GetNumLocalBasis();

    data->globalassem = new GlobalAssemblyMF(
        nLocBas, nlocalfunc, nlocalelemx, nlocalelemy);
    data->quad1 = new QuadraturePoint(p+1, 0, 1);
    data->quad2 = new QuadraturePoint(q+1, 0, 1);
    data->bernstein = new BernsteinBasis(p);
    data->elem = elemmf;
    
    data->globalassem->AssemLoad(data->quad1, data->quad2,
        data->IEN, data->ID, data->Dir, data->CP,
        data->NURBSExtraction1, data->NURBSExtraction2,
        data->elem_size1, data->elem_size2, data->elem, data->bernstein);
    
    MPI_Barrier(PETSC_COMM_WORLD);

    Mat K;
    MatCreateShell(PETSC_COMM_WORLD, nlocalfunc, nlocalfunc,
        PETSC_DETERMINE, PETSC_DETERMINE, data, &K);
    MatShellSetContext(K, data);
    MatShellSetOperation(K, MATOP_MULT, (void(*)(void))MyMatMult);
    MatShellSetOperation(K, MATOP_CREATE_VECS, (void(*)(void))MyMatCreateVecs);

    KSP ksp;
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, K, K);
    KSPSetFromOptions(ksp);

    PC pc;
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCSHELL);

    MyPCCtx *ctx;
    PetscNew(&ctx);
    KSPCreate(PETSC_COMM_WORLD, &ctx->innerksp);
    KSPSetOperators(ctx->innerksp, assembly_fem->K, assembly_fem->K);
    KSPSetType(ctx->innerksp, KSPCG);

    PC innerpc;
    KSPGetPC(ctx->innerksp, &innerpc);
    PCSetType(innerpc, PCJACOBI);
    PCSetFromOptions(innerpc);

    PetscReal rtol = 1e-10;
    PetscReal abstol = 1e-10;
    PetscReal divtol = 1e4;
    PetscInt maxits = 10000;
    KSPSetTolerances(ctx->innerksp, rtol, abstol, divtol, maxits);
    KSPSetTolerances(ksp, rtol, abstol, divtol, maxits);

    PCShellSetContext(pc, ctx);
    PCShellSetApply(pc, MyPCApply);
    PCShellSetDestroy(pc, MyPCDestroy);

    Vec u;
    VecDuplicate(data->globalassem->F, &u);
    VecSet(u, 0.0);

    PetscLogDouble tstart, tend;
    PetscTime(&tstart);
    KSPSolve(ksp, data->globalassem->F, u);
    PetscTime(&tend);
    PetscLogDouble time = tend - tstart;

    VecView(u, PETSC_VIEWER_STDOUT_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "Time: %f\n", time);

    delete fm; fm = nullptr;
    delete absgen; absgen = nullptr;
    delete iengen; iengen = nullptr;
    delete idgen; idgen = nullptr;
    delete elem2coogen; elem2coogen = nullptr;
    delete assembly_fem; assembly_fem = nullptr;
    delete elemfem; elemfem = nullptr;
    delete elemmf; elemmf = nullptr;
    delete data->quad1; data->quad1 = nullptr;
    delete data->quad2; data->quad2 = nullptr;
    delete data->bernstein; data->bernstein = nullptr;
    delete data->globalassem; data->globalassem = nullptr;

    VecDestroy(&u);
    MatDestroy(&K);
    KSPDestroy(&ksp);
    PetscFree(data);

    PetscFinalize();
    return 0;
}