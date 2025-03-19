#include "GlobalAssembly.hpp"

GlobalAssembly::GlobalAssembly(const std::vector<int> &IEN, const std::vector<int> &ID,
    const LocalAssembly * const &locassem, const int &nLocBas,
    const int &nlocalfunc, const int &nlocalelem)
    : nLocBas(nLocBas), nlocalfunc(nlocalfunc), nlocalelem(nlocalelem)
{
    const int dnz = nlocalfunc;
    const int onz = dnz;
    MatCreateAIJ(PETSC_COMM_WORLD, nlocalfunc, nlocalfunc, PETSC_DETERMINE,
        PETSC_DETERMINE, dnz, PETSC_NULL, dnz, PETSC_NULL, &K);
    VecCreate(PETSC_COMM_WORLD, &F);
    VecSetSizes(F, nlocalfunc, PETSC_DECIDE);
    VecSetFromOptions(F);
    VecSet(F, 0.0);
    VecSetOption(F, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

    AssemNonZeroEstimate(locassem, IEN, ID);

    std::vector<int> Dnz, Onz;
    NonZeroCount(K, Dnz, Onz);

    MatDestroy(&K);

    MatCreateAIJ(PETSC_COMM_WORLD, nlocalfunc, nlocalfunc, PETSC_DETERMINE,
        PETSC_DETERMINE, 0, &Dnz[0], 0, &Onz[0], &K);
    MatSetOption(mat, MAT_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
}

GlobalAssembly::~GlobalAssembly()
{
    MatDestroy(&K);
    VecDestroy(&F);
}

void GlobalAssembly::NonZeroCount(const Mat &K, std::vector<int> &dnz, std::vector<int> &onz)
{
    PetscInt rstart, rend;
    PetscInt cstart, cend;
    MatGetOwnershipRange(K, &rstart, &rend);
    MatGetOwnershipRangeColumn(K, &cstart, &cend);

    const PetscInt nlocalfunc = rend - rstart;
    dnz.resize(nlocalfunc, 0);
    onz.resize(nlocalfunc, 0);

    for (PetscInt i = 0; i < nlocalfunc; ++i)
    {
        const PetscInt globalrow = rstart + i;

        PetscInt ncols;
        const PetscInt *cols;
        MatGetRow(K, globalrow, &ncols, &cols, PETSC_NULL);

        for (PetscInt j = 0; j < ncols; ++j)
        {
            const PetscInt globalcol = cols[j];
            if (globalcol >= cstart && globalcol < cend)
            {
                ++dnz[i];
            }
            else
            {
                ++onz[i];
            }
        }
    }
}

void GlobalAssembly::AssemNonZeroEstimate(const LocalAssembly * const &locassem,
    const std::vector<int> &IEN,
    const std::vector<int> &ID)
{
    PetscInt * eID = new PetscInt[nLocBas]

    for (int i = 0; i < nlocalelem; ++i)
    {
        for (int j = 0; j < nLocBas; ++j)
        {
            eID[j] = ID[IEN[i*nLocBas+j]];
        }

        locassem->AssemNonZero();

        MatSetValues(K, nLocBas, eID, nLocBas, eID, locassem->Kloc, ADD_VALUES);
    }

    delete[] eID; eID = nullptr;

    MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
}

void GlobalAssembly::AssemStiffnessLoad(const LocalAssembly * const &locassem,
    const std::vector<int> &IEN,
    const std::vector<int> &ID,
    const std::vector<double> &CP,
    Element * const &elem)
{
    PetscInt * eID = new PetscInt[nLocBas]
    std::vector<double> eCP(2*nLocBas, 0.0);

    for (int i = 0; i < nlocalelem; ++i)
    {
        for (int j = 0; j < nLocBas; ++j)
        {
            eID[j] = ID[IEN[i*nLocBas+j]];
            eCP[2*j] = CP[2*IEN[i*nLocBas+j]];
            eCP[2*j+1] = CP[2*IEN[i*nLocBas+j]+1];
        }

        locassem->AssemLocalStiffnessLoad(elem, eCP);

        MatSetValues(K, nLocBas, eID, nLocBas, eID, locassem->Kloc, ADD_VALUES);

        VecSetValues(F, nLocBas, eID, locassem->Floc, ADD_VALUES);
    }

    delete[] eID; eID = nullptr;

    MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(F);
    VecAssemblyEnd(F);
}