#include "ControlPointGenerator.hpp"
#include "IENGenerator.hpp"
#include "IDGenerator.hpp"
#include "NURBSExtractionGenerator.hpp"
#include "FileManager.hpp"

int main(int argc, char *argv[])
{
    int p, q, nElemX, nElemY, part_num_1d, dim;
    double Lx, Ly;
    std::string base_name;

    std::string file_info = "info.txt";

    FileManager * fm = new FileManager();
    fm->ReadPreprocessInfo(file_info, p, q, Lx, Ly, nElemX, nElemY, part_num_1d, dim, base_name);

    double hx = Lx / nElemX;
    double hy = Ly / nElemY;

    int nElem = nElemX * nElemY;
    int nLocBas = (p + 1) * (q + 1);
    
    int nFuncX = p + nElemX;
    int nFuncY = q + nElemY;
    int nFunc = nFuncX * nFuncY;

    std::vector<double> S;
    std::vector<double> T;

    for (int i = 0; i < p + 1; ++i) S.push_back(0.0);
    for (int i = 1; i < nElemX; ++i) S.push_back(i * hx);
    for (int i = 0; i < p + 1; ++i) S.push_back(Lx);

    for (int i = 0; i < q + 1; ++i) T.push_back(0.0);
    for (int i = 1; i < nElemY; ++i) T.push_back(i * hy);
    for (int i = 0; i < q + 1; ++i) T.push_back(Ly);

    std::vector<double> elem_size1(nElemX, hx);
    std::vector<double> elem_size2(nElemY, hy);

    ControlPointGenerator * cpg = new ControlPointGenerator();
    BSplineBasis * basis1 = new BSplineBasis(p, S);
    BSplineBasis * basis2 = new BSplineBasis(q, T);

    double P1min = 0.0;
    double P1max = Lx;
    double P2min = 0.0;
    double P2max = Ly;

    std::cout << "Generating IEN..." << std::endl;
    IENGenerator * igen = new IENGenerator();
    std::vector<int> IEN = igen->GenerateIEN2D(basis1, basis2);

    std::cout << "Generating ID..." << std::endl;
    IDGenerator * idgen = new IDGenerator();
    std::vector<int> ID = idgen->GenerateID2D(basis1, basis2);

    std::vector<int> Dir{};
    for (int ii = 0; ii < ID.size(); ++ii)
    {
        if (ID[ii] == -1)
        Dir.push_back(ii);
    }

    std::cout << "Generating control points..." << std::endl;
    std::vector<double> CP = cpg->GenerateControlPoints2D(basis1, basis2, P1min, P1max, P2min, P2max);

    std::cout << "Generating NURBS extraction..." << std::endl;
    NURBSExtractionGenerator * neg = new NURBSExtractionGenerator();
    std::cout << "Generating extraction 1..." << std::endl;
    std::vector<double> NURBSExtraction1 = neg->GenerateExtraction1D(basis1);
    std::cout << "Generating extraction 2..." << std::endl;
    std::vector<double> NURBSExtraction2 = neg->GenerateExtraction1D(basis2);

    Element * elem = new Element(p, q);
    LocalAssembly * locassem = new LocalAssembly(p, q);
    GlobalAssembly * globalassem = new GlobalAssembly(IEN, ID, Dir, locassem, nLocBas,
            nFunc, nElemX, nElemY);

    MatSetOption(globalassem->K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

    std::vector<double> jacobian_all_qp{};
    std::vector<double> invJac_all_qp{};
    
    globalassem->AssemStiffnessLoad(locassem, IEN, ID, Dir, CP,
        NURBSExtraction1, NURBSExtraction2, elem_size1, elem_size2, elem,
        jacobian_all_qp, invJac_all_qp);

    MPI_Barrier(PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "Assembling stiffness matrix and load vector...done\n");
    PetscPrintf(PETSC_COMM_WORLD, "Stored jacobian count: %d\n",
        static_cast<int>(jacobian_all_qp.size()));
    PetscPrintf(PETSC_COMM_WORLD, "Stored invJac count: %d\n",
        static_cast<int>(invJac_all_qp.size()));
}
