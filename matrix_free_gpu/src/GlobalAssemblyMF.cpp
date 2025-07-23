#include "GlobalAssemblyMF.hpp"

GlobalAssemblyMF::GlobalAssemblyMF(const int &nLocBas, const int &nlocalfunc,
    const int &nlocalelemx, const int &nlocalelemy, const std::vector<int> &ghostID)
    : nLocBas(nLocBas), nlocalfunc(nlocalfunc),
      nlocalelemx(nlocalelemx), nlocalelemy(nlocalelemy)
{
    const int nghost = static_cast<int>(ghostID.size());

    PetscInt * ghostIdx = new PetscInt[nghost];

    for (int ii = 0; ii < nghost; ++ii)
    {
        ghostIdx[ii] = ghostID[ii];
    }
    
    VecCreateGhost(PETSC_COMM_WORLD, nlocalfunc, PETSC_DETERMINE,
        nghost, ghostIdx, &F);
    VecSetFromOptions(F);
    VecSet(F, 0.0);
    VecSetOption(F, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

    delete[] ghostIdx; ghostIdx = nullptr;
}

void GlobalAssemblyMF::AssemLoad(LocalAssemblyMF * const &locassem,
    const std::vector<int> &IEN,
    const std::vector<int> &ID,
    const std::vector<int> &Dir,
    const std::vector<double> &CP,
    const std::vector<double> &NURBSExtraction1,
    const std::vector<double> &NURBSExtraction2,
    const std::vector<double> &elem_size1,
    const std::vector<double> &elem_size2,
    ElementMF * const &elemmf)
{
    PetscInt * eID = new PetscInt[nLocBas];
    std::vector<double> eCP(2*nLocBas, 0.0);
    const int pp = elemmf->GetNumLocalBasis1D(0);
    const int qq = elemmf->GetNumLocalBasis1D(1);
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
        
            elemmf->SetElement(eNURBSExtraction1, eNURBSExtraction2, elem_size1[ii], elem_size2[jj]);

            locassem->AssemLocalLoad(elemmf, eCP);

            VecSetValues(F, nLocBas, eID, locassem->Floc, ADD_VALUES);
        }
    }

    VecAssemblyBegin(F);
    VecAssemblyEnd(F);

    delete[] eID; eID = nullptr;

    DirichletBC(Dir);

    VecAssemblyBegin(F);
    VecAssemblyEnd(F);
}

void GlobalAssemblyMF::AssemLoad(LocalAssemblyMFSF * const &locassem,
    const std::vector<int> &IEN,
    const std::vector<int> &ID,
    const std::vector<int> &Dir,
    const std::vector<double> &CP,
    const std::vector<double> &NURBSExtraction1,
    const std::vector<double> &NURBSExtraction2,
    const std::vector<double> &elem_size1,
    const std::vector<double> &elem_size2,
    ElementMFSF * const &elemmf)
{
    PetscInt * eID = new PetscInt[nLocBas];
    std::vector<double> eCP(2*nLocBas, 0.0);
    const int pp = elemmf->GetNumLocalBasis1D(0);
    const int qq = elemmf->GetNumLocalBasis1D(1);
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
        
            elemmf->SetElement(eNURBSExtraction1, eNURBSExtraction2, elem_size1[ii], elem_size2[jj]);

            locassem->AssemLocalLoad(elemmf, eCP);

            // PetscMPIInt rank;
            // MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
            // if (rank == 1)
            VecSetValues(F, nLocBas, eID, locassem->Floc, ADD_VALUES);
        }
    }

    VecAssemblyBegin(F);
    VecAssemblyEnd(F);

    delete[] eID; eID = nullptr;

    DirichletBC(Dir);

    VecAssemblyBegin(F);
    VecAssemblyEnd(F);
}

void GlobalAssemblyMF::MatMulMF(LocalAssemblyMF * const &locassem,
    const std::vector<int> &IEN,
    const std::vector<int> &ID,
    const std::vector<int> &Dir,
    const std::vector<double> &CP,
    const std::vector<double> &NURBSExtraction1,
    const std::vector<double> &NURBSExtraction2,
    const std::vector<double> &elem_size1,
    const std::vector<double> &elem_size2,
    ElementMF * const &elemmf,
    Vec x, Vec y)
{
    PetscInt * eID = new PetscInt[nLocBas];
    PetscInt * eIEN = new PetscInt[nLocBas];
    std::vector<double> eCP(2*nLocBas, 0.0);
    const int pp = elemmf->GetNumLocalBasis1D(0);
    const int qq = elemmf->GetNumLocalBasis1D(1);
    std::vector<double> eNURBSExtraction1(pp*pp, 0.0);
    std::vector<double> eNURBSExtraction2(qq*qq, 0.0);

    VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
    VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);
    Vec localx;
    VecGhostGetLocalForm(x, &localx);

    for (int jj = 0; jj < nlocalelemy; ++jj)
    {
        for (int ii = 0; ii < nlocalelemx; ++ii)
        {
            int elemIndex = jj*nlocalelemx + ii;
            for (int j = 0; j < nLocBas; ++j)
            {
                eID[j] = ID[IEN[elemIndex*nLocBas+j]];
                eIEN[j] = IEN[elemIndex*nLocBas+j];
                eCP[2*j] = CP[2*IEN[elemIndex*nLocBas+j]];
                eCP[2*j+1] = CP[2*IEN[elemIndex*nLocBas+j]+1];
            }

            std::copy(NURBSExtraction1.begin() + ii * pp * pp, 
                NURBSExtraction1.begin() + (ii + 1) * pp * pp, 
                eNURBSExtraction1.begin());
            std::copy(NURBSExtraction2.begin() + jj * qq * qq,
                NURBSExtraction2.begin() + (jj + 1) * qq * qq,
                eNURBSExtraction2.begin());
        
            elemmf->SetElement(eNURBSExtraction1, eNURBSExtraction2, elem_size1[ii], elem_size2[jj]);

            VecGetValues(localx, nLocBas, eIEN, locassem->Floc_in);

            locassem->LocalMatMulMF(elemmf, eCP);

            VecSetValues(y, nLocBas, eID, locassem->Floc_out, ADD_VALUES);
        }
    }
    VecAssemblyBegin(y);
    VecAssemblyEnd(y);

    delete[] eID; eID = nullptr;
    delete[] eIEN; eIEN = nullptr;

    const int nDir = static_cast<int>(Dir.size());
    for (int ii = 0; ii < nDir; ++ii)
    {
        VecSetValue(y, Dir[ii], 0.0, INSERT_VALUES);
    }

    VecAssemblyBegin(y);
    VecAssemblyEnd(y);
    VecGhostRestoreLocalForm(x, &localx);
}

void GlobalAssemblyMF::MatMulMF(LocalAssemblyMFSF * const &locassem,
    const std::vector<int> &IEN,
    const std::vector<int> &ID,
    const std::vector<int> &Dir,
    const std::vector<double> &CP,
    const std::vector<double> &NURBSExtraction1,
    const std::vector<double> &NURBSExtraction2,
    const std::vector<double> &elem_size1,
    const std::vector<double> &elem_size2,
    ElementMFSF * const &elemmf,
    Vec x, Vec y)
{
    PetscInt * eID = new PetscInt[nLocBas];
    PetscInt * eIEN = new PetscInt[nLocBas];
    std::vector<double> eCP(2*nLocBas, 0.0);
    const int pp = elemmf->GetNumLocalBasis1D(0);
    const int qq = elemmf->GetNumLocalBasis1D(1);
    std::vector<double> eNURBSExtraction1(pp*pp, 0.0);
    std::vector<double> eNURBSExtraction2(qq*qq, 0.0);

    VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
    VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);
    Vec localx;
    VecGhostGetLocalForm(x, &localx);

    for (int jj = 0; jj < nlocalelemy; ++jj)
    {
        for (int ii = 0; ii < nlocalelemx; ++ii)
        {
            int elemIndex = jj*nlocalelemx + ii;
            for (int j = 0; j < nLocBas; ++j)
            {
                eID[j] = ID[IEN[elemIndex*nLocBas+j]];
                eIEN[j] = IEN[elemIndex*nLocBas+j];
                eCP[2*j] = CP[2*IEN[elemIndex*nLocBas+j]];
                eCP[2*j+1] = CP[2*IEN[elemIndex*nLocBas+j]+1];
            }

            std::copy(NURBSExtraction1.begin() + ii * pp * pp, 
                NURBSExtraction1.begin() + (ii + 1) * pp * pp, 
                eNURBSExtraction1.begin());
            std::copy(NURBSExtraction2.begin() + jj * qq * qq,
                NURBSExtraction2.begin() + (jj + 1) * qq * qq,
                eNURBSExtraction2.begin());
        
            elemmf->SetElement(eNURBSExtraction1, eNURBSExtraction2, elem_size1[ii], elem_size2[jj]);

            VecGetValues(localx, nLocBas, eIEN, locassem->Floc_in);

            locassem->LocalMatMulMF(elemmf, eCP);

            VecSetValues(y, nLocBas, eID, locassem->Floc_out, ADD_VALUES);
        }
    }
    VecAssemblyBegin(y);
    VecAssemblyEnd(y);

    delete[] eID; eID = nullptr;
    delete[] eIEN; eIEN = nullptr;

    const int nDir = static_cast<int>(Dir.size());
    for (int ii = 0; ii < nDir; ++ii)
    {
        VecSetValue(y, Dir[ii], 0.0, INSERT_VALUES);
    }

    VecAssemblyBegin(y);
    VecAssemblyEnd(y);
    VecGhostRestoreLocalForm(x, &localx);
}

void GlobalAssemblyMF::DirichletBC(const std::vector<int> &Dir)
{
    const int nDir = static_cast<int>(Dir.size());
    for (int ii = 0; ii < nDir; ++ii)
    {
        VecSetValue(F, Dir[ii], 0.0, INSERT_VALUES);
    }
}