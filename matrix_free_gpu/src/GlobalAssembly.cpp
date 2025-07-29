#include "GlobalAssembly.hpp"

__global__ void AssembleKernel(const int nLocBas,
    const double * d_N, const double * d_dN_dxi, const double * d_dN_deta,
    const double * d_weight,
    const int * d_IEN, const int * d_ID, const int * d_Dir,
    const double * d_CP, double * d_val)
{
    extern __shared__ double shared_data[];

    int offset = 0;
    PetscInt * s_eID = (PetscInt*)(shared_data + offset);
    offset += nLocBas * (sizeof(PetscInt)/sizeof(double));
    double * s_eCP = shared_data + offset;

    int elemIndex = blockIdx.y * gridDim.x + blockIdx.x;

    for (int j = 0; j < nLocBas; ++j)
    {
        s_eID[j] = d_ID[d_IEN[elemIndex * nLocBas + j]];
        s_eCP[2 * j] = d_CP[2 * d_IEN[elemIndex * nLocBas + j]];
        s_eCP[2 * j + 1] = d_CP[2 * d_IEN[elemIndex * nLocBas + j] + 1];
    }

    __syncthreads();


}

__device__ void ElementKernel(const int nLocBas,
    const double * d_N, const double * d_dN_dxi, const double * d_dN_deta,
    const double * d_weight,
    const double * s_eCP)
{
    
}

GlobalAssembly::GlobalAssembly(const int &nlocalfunc,
    const int &nnz, const std::vector<int> &rows,
    const std::vector<int> &cols)
: nLocBas(nlocalfunc), nnz(nnz)
{
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
}

GlobalAssembly::~GlobalAssembly()
{
    MatDestroy(&K);
}

void GlobalAssembly::AssemStiffnessLoad(LocalAssembly * const &locassem,
    const std::vector<int> &IEN,
    const std::vector<int> &ID,
    const std::vector<int> &Dir,
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

    int * d_IEN, * d_ID, * d_Dir;
    double * d_CP, * d_NURBSExtraction1, * d_NURBSExtraction2;
    double * d_elem_size1, * d_elem_size2;
    cudaMalloc((void**)&d_IEN, IEN.size() * sizeof(int));
    cudaMalloc((void**)&d_ID, ID.size() * sizeof(int));
    cudaMalloc((void**)&d_Dir, Dir.size() * sizeof(int));
    cudaMalloc((void**)&d_CP, CP.size() * sizeof(double));
    cudaMalloc((void**)&d_NURBSExtraction1, NURBSExtraction1.size() * sizeof(double));
    cudaMalloc((void**)&d_NURBSExtraction2, NURBSExtraction2.size() * sizeof(double));
    cudaMalloc((void**)&d_elem_size1, elem_size1.size() * sizeof(double));
    cudaMalloc((void**)&d_elem_size2, elem_size2.size() * sizeof(double));

    cudaMemcpy(d_IEN, IEN.data(), IEN.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ID, ID.data(), ID.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Dir, Dir.data(), Dir.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_CP, CP.data(), CP.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_NURBSExtraction1, NURBSExtraction1.data(), NURBSExtraction1.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_NURBSExtraction2, NURBSExtraction2.data(), NURBSExtraction2.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_elem_size1, elem_size1.data(), elem_size1.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_elem_size2, elem_size2.data(), elem_size2.size() * sizeof(double), cudaMemcpyHostToDevice);

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

            MatSetValues(K, nLocBas, eID, nLocBas, eID, locassem->Kloc, ADD_VALUES);

            VecSetValues(F, nLocBas, eID, locassem->Floc, ADD_VALUES);
        }
    }
    VecAssemblyBegin(F);
    VecAssemblyEnd(F);

    delete[] eID; eID = nullptr;

    DirichletBC(Dir);

    MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(F);
    VecAssemblyEnd(F);

    cudaFree(d_IEN);
    cudaFree(d_ID);
    cudaFree(d_Dir);
    cudaFree(d_CP);
    cudaFree(d_NURBSExtraction1);
    cudaFree(d_NURBSExtraction2);
    cudaFree(d_elem_size1);
    cudaFree(d_elem_size2);
}

void GlobalAssembly::AssemStiffnessLoad(QuadraturePoint * const &quad1,
    QuadraturePoint * const &quad2,
    const std::vector<int> &IEN,
    const std::vector<int> &ID,
    const std::vector<int> &Dir,
    const std::vector<double> &CP,
    ElementFEM * const &elem)
{
    PetscInt * eID = new PetscInt[nLocBas];
    std::vector<double> eCP(2*nLocBas, 0.0);

    const int nqp1 = quad1->GetNumQuadraturePoint();
    const int nqp2 = quad2->GetNumQuadraturePoint();
    const std::vector<double> qp1 = quad1->GetQuadraturePoint();
    const std::vector<double> qp2 = quad2->GetQuadraturePoint();
    const std::vector<double> w1 = quad1->GetWeight();
    const std::vector<double> w2 = quad2->GetWeight();
    const int n = elem->GetNumLocalBasis();

    std::vector<double> sqp_N{};
    std::vector<double> sqp_dN_dxi{};
    std::vector<double> sqp_dN_deta{};

    std::vector<double> N{};
    std::vector<double> dN_dxi{};
    std::vector<double> dN_deta{};
    std::vector<double> weight{};

    for (int jj = 0; jj < nqp2; ++jj)
    {
        for (int ii = 0; ii < nqp1; ++ii)
        {
            elem->GenerateRefElementSingleQP(qp1[ii], qp2[jj], sqp_N, sqp_dN_dxi, sqp_dN_deta);
            N.push_back(sqp_N);
            dN_dxi.push_back(sqp_dN_dxi);
            dN_deta.push_back(sqp_dN_deta);
            weight.push_back(w1[ii] * w2[jj]);
        }
    }

    double * d_N, * d_dN_dxi, * d_dN_deta, * d_weight;
    cudaMalloc((void**)&d_N, N.size() * sizeof(double));
    cudaMalloc((void**)&d_dN_dxi, dN_dxi.size() * sizeof(double));
    cudaMalloc((void**)&d_dN_deta, dN_deta.size() * sizeof(double));
    cudaMalloc((void**)&d_weight, weight.size() * sizeof(double));
    cudaMemcpy(d_N, N.data(), N.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_dN_dxi, dN_dxi.data(), dN_dxi.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_dN_deta, dN_deta.data(), dN_deta.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_weight, weight.data(), weight.size() * sizeof(double), cudaMemcpyHostToDevice);

    int * d_IEN, * d_ID, * d_Dir;
    double * d_CP;
    cudaMalloc((void**)&d_IEN, IEN.size() * sizeof(int));
    cudaMalloc((void**)&d_ID, ID.size() * sizeof(int));
    cudaMalloc((void**)&d_Dir, Dir.size() * sizeof(int));
    cudaMalloc((void**)&d_CP, CP.size() * sizeof(double));
    cudaMemcpy(d_IEN, IEN.data(), IEN.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ID, ID.data(), ID.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Dir, Dir.data(), Dir.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_CP, CP.data(), CP.size() * sizeof(double), cudaMemcpyHostToDevice);

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

            locassem->AssemLocalStiffnessLoad(elem, eCP);

            MatSetValues(K, nLocBas, eID, nLocBas, eID, locassem->Kloc, ADD_VALUES);

            VecSetValues(F, nLocBas, eID, locassem->Floc, ADD_VALUES);
        }
    }
    VecAssemblyBegin(F);
    VecAssemblyEnd(F);

    delete[] eID; eID = nullptr;

    DirichletBC(Dir);

    MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(F);
    VecAssemblyEnd(F);

    cudaFree(d_N);
    cudaFree(d_dN_dxi);
    cudaFree(d_dN_deta);
    cudaFree(d_weight);

    cudaFree(d_IEN);
    cudaFree(d_ID);
    cudaFree(d_Dir);
    cudaFree(d_CP);
}

void GlobalAssembly::DirichletBCK(const std::vector<int> &Dir)
{
    for (int i = 0; i < static_cast<int>(Dir.size()); ++i)
        MatSetValue(K, Dir[i], Dir[i], 1.0, ADD_VALUES);
}

void GlobalAssembly::DirichletBC(const std::vector<int> &Dir)
{
    for (int i = 0; i < static_cast<int>(Dir.size()); ++i)
    {
        MatSetValue(K, Dir[i], Dir[i], 1.0, ADD_VALUES);
        VecSetValue(F, Dir[i], 0.0, INSERT_VALUES);
    }
}