#include "GlobalAssembly.hpp"

__global__ void AssembleKernel(const int nLocBas, const int nqp,
    const double * d_N, const double * d_dN_dxi, const double * d_dN_deta,
    const double * d_weight,
    const int * d_IEN, const int * d_ID,
    const double * d_CP, const int * elem2coo,
    double * d_val)
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

    int qp = threadIdx.y * blockDim.x + threadIdx.x;

    double R[nLocBas];
    double dN_dxi_q[nLocBas];
    double dN_deta_q[nLocBas];

    if (qp < nqp)
    {
        for (int i = 0; i < nLocBas; ++i)
        {
            R[i] = d_N[qp * nLocBas + i];
            dN_dxi_q[i] = d_dN_dxi[qp * nLocBas + i];
            dN_deta_q[i] = d_dN_deta[qp * nLocBas + i];
        }

        double jacobian;
        double dR_dx[nLocBas];
        double dR_dy[nLocBas];
        compute_jacobian_basis_derivative(nLocBas, dN_dxi_q, dN_deta_q, s_eCP, jacobian, dR_dx, dR_dy);

        for (int i = 0; i < nLocBas; ++i)
        {
            for (int j = 0; j < nLocBas; ++j)
            {
                int elem2coo_index = elem2coo[elemIndex * nLocBas * nLocBas + i * nLocBas + j];
                double val = (dR_dx[i] * dR_dx[j] + dR_dy[i] * dR_dy[j]) * jacobian * d_weight[qp];
                atomicAdd(&d_val[elem2coo_index], val);
            }
        }
    }
}

__global__ void DirichletBCKernel(const int * dir2coo, const int dirsize, double * d_val)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < dirsize)
    {
        int coo_index = dir2coo[idx];
        if (coo_index >= 0)
        {
            d_val[coo_index] = 1.0;
        }
    }
}

__device__ void compute_jacobian_basis_derivative(
    const int nLocBas, const double * dN_dxi, const double * dN_deta,
    const double * eCP, double &jacobian, double * dR_dx, double * dR_dy)
{
    double dx_dxi = 0.0;
    double dx_deta = 0.0;
    double dy_dxi = 0.0;
    double dy_deta = 0.0;

    jacobian = 0.0;

    double dxi_dx = 0.0;
    double dxi_dy = 0.0;
    double deta_dx = 0.0;
    double deta_dy = 0.0;

    for (int ii = 0; ii < nLocBas; ++ii)
    {
        dx_dxi += eCP[ii * 2] * dN_dxi[ii];
        dx_deta += eCP[ii * 2] * dN_deta[ii];
        dy_dxi += eCP[ii * 2 + 1] * dN_dxi[ii];
        dy_deta += eCP[ii * 2 + 1] * dN_deta[ii];
    }

    jacobian = dx_dxi*dy_deta - dx_deta*dy_dxi;

    dxi_dx = dy_deta/jacobian;
    dxi_dy = -dy_dxi/jacobian;
    deta_dx = -dx_deta/jacobian;
    deta_dy = dx_dxi/jacobian;

    for (int ii = 0; ii < nLocBas; ++ii)
    {
        dR_dx[ii] = dxi_dx * dN_dxi[ii] + dxi_dy * dN_deta[ii];
        dR_dy[ii] = deta_dx * dN_dxi[ii] + deta_dy * dN_deta[ii];
    }
}

GlobalAssembly::GlobalAssembly(const int &nLocBas,
    const int &nnz, const int &nlocalfunc,
    const int &nlocalelemx, const int &nlocalelemy,
    const std::vector<int> &rows,
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

void GlobalAssembly::AssemStiffness(QuadraturePoint * const &quad1,
    QuadraturePoint * const &quad2,
    const std::vector<int> &IEN,
    const std::vector<int> &ID,
    const std::vector<int> &dir2coo,
    const std::vector<double> &CP,
    const std::vector<int> &elem2coo,
    ElementFEM * const &elem)
{
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

    int * d_IEN, * d_ID, * d_dir2coo, * d_elem2coo;
    double * d_CP, * d_val;
    cudaMalloc((void**)&d_IEN, IEN.size() * sizeof(int));
    cudaMalloc((void**)&d_ID, ID.size() * sizeof(int));
    cudaMalloc((void**)&d_dir2coo, dir2coo.size() * sizeof(int));
    cudaMalloc((void**)&d_CP, CP.size() * sizeof(double));
    cudaMalloc((void**)&d_elem2coo, elem2coo.size() * sizeof(int));
    cudaMalloc((void**)&d_val, nnz * sizeof(double));
    cudaMemcpy(d_IEN, IEN.data(), IEN.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ID, ID.data(), ID.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Dir, dir2coo.data(), dir2coo.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_CP, CP.data(), CP.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_elem2coo, elem2coo.data(), elem2coo.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemset(d_val, 0, nnz * sizeof(double));

    AssembleKernel<<<dim3(nlocalelemx, nlocalelemy), dim3(nqp1, nqp2), nLocBas * sizeof(PetscInt) + 2 * nLocBas * sizeof(double)>>>(
        nLocBas, nqp1 * nqp2, d_N, d_dN_dxi, d_dN_deta, d_weight,
        d_IEN, d_ID, d_Dir, d_CP, d_elem2coo, d_val);

    DirichletBCKernel<<<(Dir.size() + 255) / 256, 256>>>(d_dir2coo, dir2coo.size(), d_val);

    MarSetValuesCOO(K, d_val, INSERT_VALUES);

    cudaFree(d_N);
    cudaFree(d_dN_dxi);
    cudaFree(d_dN_deta);
    cudaFree(d_weight);

    cudaFree(d_IEN);
    cudaFree(d_ID);
    cudaFree(d_Dir);
    cudaFree(d_CP);
    cudaFree(d_elem2coo);
    cudaFree(d_val);
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