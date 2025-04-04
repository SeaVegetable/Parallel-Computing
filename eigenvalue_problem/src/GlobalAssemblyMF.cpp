#include "GlobalAssemblyMF.hpp"

GlobalAssemblyMF::GlobalAssemblyMF(const std::vector<int> &IEN, const std::vector<int> &ID,
    LocalAssembly * const &locassem, const int &nLocBas,
    const int &nlocalfunc, const int &nlocalelemx, const int &nlocalelemy)
    : nLocBas(nLocBas), nlocalfunc(nlocalfunc),
      nlocalelemx(nlocalelemx), nlocalelemy(nlocalelemy)
{
    VecCreate(PETSC_COMM_WORLD, &F);
    VecSetSizes(F, nlocalfunc, PETSC_DECIDE);
    VecSetFromOptions(F);
    VecSet(F, 0.0);
    VecSetOption(F, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
}

GlobalAssemblyMF::~GlobalAssemblyMF()
{
    VecDestroy(&F);
}

void GlobalAssemblyMF::AssemStiffnessLoad(LocalAssembly * const &locassem,
    const std::vector<int> &IEN,
    const std::vector<int> &ID,
    const std::vector<double> &CP,
    const std::vector<double> &NURBSExtraction1,
    const std::vector<double> &NURBSExtraction2,
    const std::vector<double> &elem_size1,
    const std::vector<double> &elem_size2,
    Element * const &elem)
{
    PetscInt * eID = new PetscInt[nLocBas];
    std::vector<double> eCP(2*nLocBas, 0.0);
    const int pp = elem->GetNumLocalBasis1D(0);
    const int qq = elem->GetNumLocalBasis1D(1);
    std::vector<double> eNURBSExtraction1(pp*pp, 0.0);
    std::vector<double> eNURBSExtraction2(qq*qq, 0.0);

    for (int jj = 0; jj < nlocalelemy; ++jj)
    {
        for (int ii = 0; ii < nlocalelemx; ++ii)
        {
            int elemIndex = jj*nlocalelemx + ii;
            for (int j = 0; j < nLocBas; ++j)
            {
                eID[j] = ID[IEN[elemIndex*nLocBas+j]];
                eCP[2*j] = CP[2*IEN[elemIndex*nLocBas+j]];
                eCP[2*j+1] = CP[2*IEN[elemIndex*nLocBas+j]+1];
            }

            std::copy(NURBSExtraction1.begin() + ii * pp * pp, 
                NURBSExtraction1.begin() + (ii + 1) * pp * pp, 
                eNURBSExtraction1.begin());
            std::copy(NURBSExtraction2.begin() + jj * qq * qq,
                NURBSExtraction2.begin() + (jj + 1) * qq * qq,
                eNURBSExtraction2.begin());
        
            elem->SetElement(eNURBSExtraction1, eNURBSExtraction2, elem_size1[ii], elem_size2[jj]);

            locassem->AssemLocalStiffnessLoad(elem, eCP);

            VecSetValues(F, nLocBas, eID, locassem->Floc, ADD_VALUES);
        }
    }

    delete[] eID; eID = nullptr;

    VecAssemblyBegin(F);
    VecAssemblyEnd(F);
}