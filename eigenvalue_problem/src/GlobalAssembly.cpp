#include "GlobalAssembly.hpp"

GlobalAssembly::GlobalAssembly(const int &nlocalfunc, const int &nLocBas)
    : nlocalfunc(nlocalfunc)
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