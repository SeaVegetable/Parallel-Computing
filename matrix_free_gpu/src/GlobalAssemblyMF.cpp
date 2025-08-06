#include "GlobalAssemblyMF.hpp"

__global__ void AssembleKernel(const int p, const int q,
    double *d_B1, double *d_B2,
    double *d_dB1, double *d_dB2,
    double *d_nurbs_extraction1, double *d_nurbs_extraction2,
    double *d_elem_size1, double *d_elem_size2,
    int *d_IEN, int *d_ID,
    double *d_CP,
    double *qw1, double *qw2,
    double *d_x_array
    )
{
    extern __shared__ double shared_data[];

    int offset = 0;
    PetscInt *s_eID = (PetscInt*)(shared_data + offset);
    offset += (p + 1) * (q + 1) * (sizeof(PetscInt) / sizeof(double));
    double *s_eCP = shared_data + offset;
    offset += (p + 1) * (p + 1);
    double *s_eNURBSExtraction1 = shared_data + offset;
    offset += (q + 1) * (q + 1);
    double *s_eNURBSExtraction2 = shared_data + offset;
    offset += (p + 1) * (q + 1);
    double *s_qw = shared_data + offset;
    
    int elemIndex = blockIdx.y * gridDim.x + blockIdx.x;
    int nLocBas = (p + 1) * (q + 1);

    for (int j = 0; j < nLocBas; ++j)
    {
        s_eID[j] = d_ID[d_IEN[elemIndex * nLocBas + j]];
        s_eCP[2 * j] = d_CP[2 * d_IEN[elemIndex * nLocBas + j]];
        s_eCP[2 * j + 1] = d_CP[2 * d_IEN[elemIndex * nLocBas + j] + 1];
    }

    for (int j = 0; j < q + 1; ++j)
    {
        for (int i = 0; i < p + 1; ++i)
        {
            s_qw[j * (p + 1) + i] = qw1[i] * qw2[j];
        }
    }

    for (int i = 0; i < (p + 1) * (p + 1); ++i)
        s_eNURBSExtraction1[i] = d_nurbs_extraction1[elemIndex * (p + 1) * (p + 1) + i];
    for (int i = 0; i < (q + 1) * (q + 1); ++i)
        s_eNURBSExtraction2[i] = d_nurbs_extraction2[elemIndex * (q + 1) * (q + 1) + i];

    double h1 = d_elem_size1[blockIdx.x];
    double h2 = d_elem_size2[blockIdx.y];

    __syncthreads();

    int qpx = threadIdx.x;
    int qpy = threadIdx.y;
    int qp = threadIdx.y * blockDim.x + threadIdx.x;

    double B1[p + 1];
    double dB1[p + 1];
    double B2[q + 1];
    double dB2[q + 1];

    if (qp < nqp)
    {
        for (int i = 0; i < p + 1; ++i)
        {
            B1[i] = d_B1[qpx * (p + 1) + i];
            dB1[i] = d_dB1[qpx * (p + 1) + i];
        }
        for (int j = 0; j < q + 1; ++j)
        {
            B2[j] = d_B2[qpy * (q + 1) + j];
            dB2[j] = d_dB2[qpy * (q + 1) + j];
        }

        double jacobian;
        double R[nLocBas];

        compute_jacobian_basis(p, q, h1, h2, B1, B2, dB1, dB2,
            s_eNURBSExtraction1, s_eNURBSExtraction2, s_eCP, jacobian, R);
        
        double x = 0.0;
        double y = 0.0;

        for (int ii = 0; ii < nLocBas; ++ii)
        {
            x += s_eCP[2 * ii] * R[ii];
            y += s_eCP[2 * ii + 1] * R[ii];
        }

        for (int ii = 0; ii < nLocBas; ++ii)
        {
            double val = R[ii] * jacobian * s_qw[qp];
            int coo_index = s_eID[ii];
            if (coo_index >= 0)
            {
                atomicAdd(&d_x_array[coo_index], val);
            }
        }
    }
}

__global__ void MatrixFreeMatMultKernel(const int p, const int q,
    double *d_B1, double *d_B2,
    double *d_dB1, double *d_dB2,
    double *d_nurbs_extraction1, double *d_nurbs_extraction2,
    double *d_elem_size1, double *d_elem_size2,
    int *d_IEN, int *d_ID,
    double *d_CP,
    double *qw1, double *qw2,
    double *d_F_array_in,
    double *d_F_array_out
    )
{
    extern __shared__ double shared_data[];

    int nLocBas = (p + 1) * (q + 1);

    int offset = 0;
    PetscInt *s_eID = (PetscInt*)(shared_data + offset);
    offset += nLocBas * (sizeof(PetscInt) / sizeof(double));
    double *s_eCP = shared_data + offset;
    offset += (p + 1) * (p + 1);
    double *s_eNURBSExtraction1 = shared_data + offset;
    offset += (q + 1) * (q + 1);
    double *s_eNURBSExtraction2 = shared_data + offset;
    offset += nLocBas;
    double *s_qw = shared_data + offset;
    offset += nLocBas;
    double *Floc_in = shared_data + offset;
    offset += nLocBas;
    double *Floc_out = shared_data + offset;
    
    int elemIndex = blockIdx.y * gridDim.x + blockIdx.x;

    for (int j = 0; j < nLocBas; ++j)
    {
        s_eID[j] = d_ID[d_IEN[elemIndex * nLocBas + j]];
        s_eCP[2 * j] = d_CP[2 * d_IEN[elemIndex * nLocBas + j]];
        s_eCP[2 * j + 1] = d_CP[2 * d_IEN[elemIndex * nLocBas + j] + 1];
    }

    for (int j = 0; j < q + 1; ++j)
    {
        for (int i = 0; i < p + 1; ++i)
        {
            s_qw[j * (p + 1) + i] = qw1[i] * qw2[j];
        }
    }

    for (int i = 0; i < (p + 1) * (p + 1); ++i)
        s_eNURBSExtraction1[i] = d_nurbs_extraction1[elemIndex * (p + 1) * (p + 1) + i];
    for (int i = 0; i < (q + 1) * (q + 1); ++i)
        s_eNURBSExtraction2[i] = d_nurbs_extraction2[elemIndex * (q + 1) * (q + 1) + i];

    for (int i = 0; i < nLocBas; ++i)
    {
        int coo_index = d_IEN[elemIndex * nLocBas + i];
        Floc_in[i] = d_F_array_in[coo_index];
        Floc_out[i] = 0.0;
    }

    double h1 = d_elem_size1[blockIdx.x];
    double h2 = d_elem_size2[blockIdx.y];

    __syncthreads();

    int qpx = threadIdx.x;
    int qpy = threadIdx.y;
    int qp = threadIdx.y * blockDim.x + threadIdx.x;

    double B1[p + 1];
    double dB1[p + 1];
    double B2[q + 1];
    double dB2[q + 1];

    if (qp < nqp)
    {
        for (int i = 0; i < p + 1; ++i)
        {
            B1[i] = d_B1[qpx * (p + 1) + i];
            dB1[i] = d_dB1[qpx * (p + 1) + i];
        }
        for (int j = 0; j < q + 1; ++j)
        {
            B2[j] = d_B2[qpy * (q + 1) + j];
            dB2[j] = d_dB2[qpy * (q + 1) + j];
        }

        double jacobian;
        double dR_dx[nLocBas];
        double dR_dy[nLocBas];

        compute_jacobian_derivative(p, q, h1, h2, B1, B2, dB1, dB2,
            s_eNURBSExtraction1, s_eNURBSExtraction2, s_eCP, jacobian, dR_dx, dR_dy);
        
        double temp_x = 0.0;
        double temp_y = 0.0;

        for (int jj = 0; jj < nLocBas; ++jj)
        {
            temp_x += dR_dx(jj) * Floc_in[jj];
            temp_y += dR_dy(jj) * Floc_in[jj];
        }

        temp_x *= -s_qw(qp) * jacobian;
        temp_y *= -s_qw(qp) * jacobian;

        for (int ii = 0; ii < nLocBas; ++ii)
        {
            Floc_out[ii] += (dR_dx[ii] * temp_x + dR_dy[ii] * temp_y);
        }

        __syncthreads();

        for (int ii = 0; ii < nLocBas; ++ii)
        {
            int coo_index = s_eID[ii];
            if (coo_index >= 0)
            {
                atomicAdd(&d_F_array_out[coo_index], Floc_out[ii]);
            }
        }
    }
}

__global__ void DirichletBCKernel(int * d_Dir, const int dirsize, double * d_val)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < dirsize)
    {
        int coo_index = d_Dir[idx];
        if (coo_index >= 0)
        {
            d_val[coo_index] = 0.0;
        }
    }
}

__device__ void compute_jacobian_basis(
    const int p, const int q, const double h1, const double h2,
    const double *d_B1, const double *d_B2,
    const double *d_dB1, const double *d_dB2,
    const double *s_nurbs_extraction1, const double *s_nurbs_extraction2,
    const double *eCP,
    double &jacobian,
    double *R)
{
    double N1[p + 1];
    double dN1[p + 1];
    for (int i = 0; i < p + 1; ++i)
    {
        N1[i] = 0.0;
        dN1[i] = 0.0;
    }
    double N2[q + 1];
    double dN2[q + 1];
    for (int j = 0; j < q + 1; ++j)
    {
        N2[j] = 0.0;
        dN2[j] = 0.0;
    }

    for (int jj = 0; jj < p + 1; ++jj)
    {
        for (int kk = 0; kk < p + 1; ++kk)
        {
            N1[jj] += s_nurbs_extraction1[jj * (p + 1) + kk] * d_B1[kk];
            dN1[jj] += s_nurbs_extraction1[jj * (p + 1) + kk] * d_dB1[kk];
        }
        dN1[jj] /= h1;
    }
    for (int jj = 0; jj < q + 1; ++jj)
    {
        for (int kk = 0; kk < q + 1; ++kk)
        {
            N2[jj] += s_nurbs_extraction2[jj * (q + 1) + kk] * d_B2[kk];
            dN2[jj] += s_nurbs_extraction2[jj * (q + 1) + kk] * d_dB2[kk];
        }
        dN2[jj] /= h2;
    }

    int nLocBas = (p + 1) * (q + 1);
    double N[nLocBas];
    double dN_dxi[nLocBas];
    double dN_deta[nLocBas];
    double w = 0.0;
    double dw_dxi = 0.0;
    double dw_deta = 0.0;
    for (int j = 0; j < q + 1; ++j)
    {
        for (int i = 0; i < p + 1; ++i)
        {
            N[j * (p + 1) + i] = N1[i] * N2[j];
            w += N[j * (p + 1) + i];
            dN_dxi[j * (p + 1) + i] = dN1[i] * N2[j];
            dw_dxi += dN_dxi[j * (p + 1) + i];
            dN_deta[j * (p + 1) + i] = N1[i] * dN2[j];
            dw_deta += dN_deta[j * (p + 1) + i];
        }
    }

    double dR_dxi[nLocBas];
    double dR_deta[nLocBas];
    for (int j = 0; j < q + 1; ++j)
    {
        for (int i = 0; i < p + 1; ++i)
        {
            R[j * (p + 1) + i] = N[j * (p + 1) + i] / w;
            dR_dxi[j * (p + 1) + i] = (dN_dxi[j * (p + 1) + i] - dw_dxi * R[j * (p + 1) + i]) / w;
            dR_deta[j * (p + 1) + i] = (dN_deta[j * (p + 1) + i] - dw_deta * R[j * (p + 1) + i]) / w;
        }
    }

    double dx_dxi = 0.0;
    double dx_deta = 0.0;
    double dy_dxi = 0.0;
    double dy_deta = 0.0;

    for (int j = 0; j < q + 1; ++j)
    {
        for (int i = 0; i < p + 1; ++i)
        {
            dx_dxi += eCP[2 * (j * (p + 1) + i)] * dR_dxi[j * (p + 1) + i];
            dx_deta += eCP[2 * (j * (p + 1) + i)] * dR_deta[j * (p + 1) + i];
            dy_dxi += eCP[2 * (j * (p + 1) + i) + 1] * dR_dxi[j * (p + 1) + i];
            dy_deta += eCP[2 * (j * (p + 1) + i) + 1] * dR_deta[j * (p + 1) + i];
        }
    }
    jacobian = dx_dxi * dy_deta - dx_deta * dy_dxi;
}

__device__ void compute_jacobian_derivative(
    const int p, const int q, const double h1, const double h2,
    const double *d_B1, const double *d_B2,
    const double *d_dB1, const double *d_dB2,
    const double *s_nurbs_extraction1, const double *s_nurbs_extraction2,
    const double *eCP, 
    double &jacobian,
    double *dR_dx, 
    double *dR_dy)
{
    double N1[p + 1];
    double dN1[p + 1];
    for (int i = 0; i < p + 1; ++i)
    {
        N1[i] = 0.0;
        dN1[i] = 0.0;
    }
    double N2[q + 1];
    double dN2[q + 1];
    for (int j = 0; j < q + 1; ++j)
    {
        N2[j] = 0.0;
        dN2[j] = 0.0;
    }

    for (int jj = 0; jj < p + 1; ++jj)
    {
        for (int kk = 0; kk < p + 1; ++kk)
        {
            N1[jj] +=  s_nurbs_extraction1[jj * (p + 1) + kk] * B1[kk];
            dN1[jj] += s_nurbs_extraction1[jj * (p + 1) + kk] * dB1[kk];
        }
        dN1[jj] /= h1;
    }
    for (int jj = 0; jj < q + 1; ++jj)
    {
        for (int kk = 0; kk < q + 1; ++kk)
        {
            N2[jj] += s_nurbs_extraction2[jj * (q + 1) + kk] * B2[kk];
            dN2[jj] += s_nurbs_extraction2[jj * (q + 1) + kk] * dB2[kk];
        }
        dN2[jj] /= h2;
    }

    int nLocBas = (p + 1) * (q + 1);
    double N[nLocBas];
    double dN_dxi[nLocBas];
    double dN_deta[nLocBas];
    double w = 0.0;
    double dw_dxi = 0.0;
    double dw_deta = 0.0;

    for (int j = 0; j < q + 1; ++j)
    {
        for (int i = 0; i < p + 1; ++i)
        {
            N[j * (p + 1) + i] = N1[i] * N2[j];
            w += N[j * (p + 1) + i];
            dN_dxi[j * (p + 1) + i] = dN1[i] * N2[j];
            dw_dxi += dN_dxi[j * (p + 1) + i];
            dN_deta[j * (p + 1) + i] = N1[i] * dN2[j];
            dw_deta += dN_deta[j * (p + 1) + i];
        }
    }

    double dR_dxi[nLocBas];
    double dR_deta[nLocBas];

    for (int j = 0; j < q + 1; ++j)
    {
        for (int i = 0; i < p + 1; ++i)
        {
            R[j * (p + 1) + i] = N[j * (p + 1) + i] / w;
            dR_dxi[j * (p + 1) + i] = (dN_dxi[j * (p + 1) + i] - dw_dxi * R[j * (p + 1) + i]) / w;
            dR_deta[j * (p + 1) + i] = (dN_deta[j * (p + 1) + i] - dw_deta * R[j * (p + 1) + i]) / w;
        }
    }

    double dx_dxi = 0.0;
    double dx_deta = 0.0;
    double dy_dxi = 0.0;
    double dy_deta = 0.0;
    double dxi_dx = 0.0;
    double dxi_dy = 0.0;
    double deta_dx = 0.0;
    double deta_dy = 0.0;

    for (int j = 0; j < q + 1; ++j)
    {
        for (int i = 0; i < p + 1; ++i)
        {
            dx_dxi += eCP[2 * (j * (p + 1) + i)] * dR_dxi[j * (p + 1) + i];
            dx_deta += eCP[2 * (j * (p + 1) + i)] * dR_deta[j * (p + 1) + i];
            dy_dxi += eCP[2 * (j * (p + 1) + i) + 1] * dR_dxi[j * (p + 1) + i];
            dy_deta += eCP[2 * (j * (p + 1) + i) + 1] * dR_deta[j * (p + 1) + i];
        }
    }

    jacobian = dx_dxi * dy_deta - dx_deta * dy_dxi;

    dxi_dx = dy_deta / jacobian;
    dxi_dy = -dx_deta / jacobian;
    deta_dx = -dy_dxi / jacobian;
    deta_dy = dx_dxi / jacobian;

    for (int j = 0; j < q + 1; ++j)
    {
        for (int i = 0; i < p + 1; ++i)
        {
            dR_dx[j * (p + 1) + i] = dxi_dx * dR_dxi[j * (p + 1) + i] + deta_dx * dR_deta[j * (p + 1) + i];
            dR_dy[j * (p + 1) + i] = dxi_dy * dR_dxi[j * (p + 1) + i] + deta_dy * dR_deta[j * (p + 1) + i];
        }
    }

    jacobian *= h1 * h2;
}

GlobalAssemblyMF::GlobalAssemblyMF(const int &nLocBas, const int &nlocalfunc,
    const int &nlocalelemx, const int &nlocalelemy)
    : nLocBas(nLocBas), nlocalfunc(nlocalfunc),
      nlocalelemx(nlocalelemx), nlocalelemy(nlocalelemy)
{   
    VecCreate(PETSC_COMM_WORLD, &F);
    VecSetSizes(F, nlocalfunc, PETSC_DETERMINE);
    VecSetType(F, VECCUDA);
    VecSet(F, 0.0);
}

void GlobalAssemblyMF::AssemLoad(QuadraturePoint * const &quad1,
    QuadraturePoint * const &quad2,
    const std::vector<int> &IEN,
    const std::vector<int> &ID,
    const std::vector<int> &Dir,
    const std::vector<double> &CP,
    const std::vector<double> &NURBSExtraction1,
    const std::vector<double> &NURBSExtraction2,
    const std::vector<double> &elem_size1,
    const std::vector<double> &elem_size2,
    ElementMF * const &elemmf,
    BernsteinBasis * const &bernstein)
{
    const int nqp1 = quad1->GetNumQuadraturePoint();
    const int nqp2 = quad2->GetNumQuadraturePoint();
    const std::vector<double> qp1 = quad1->GetQuadraturePoint();
    const std::vector<double> qp2 = quad2->GetQuadraturePoint();
    const std::vector<double> w1 = quad1->GetWeight();
    const std::vector<double> w2 = quad2->GetWeight();
    const int nLocBas = elemmf->GetNumLocalBasis();

    std::vector<double> B1{};
    std::vector<double> B2{};
    std::vector<double> dB1{};
    std::vector<double> dB2{};

    for (int i = 0; i < nqp1; ++i)
    {
        B1.push_back(bernstein->GetBernsteinBasisSingleQP(qp1[i]));
        dB1.push_back(bernstein->GetBernsteinBasisDerivativeSingleQP(qp1[i]));
    }
    for (int j = 0; j < nqp2; ++j)
    {
        B2.push_back(bernstein->GetBernsteinBasisSingleQP(qp2[j]));
        dB2.push_back(bernstein->GetBernsteinBasisDerivativeSingleQP(qp2[j]));
    }

    double * d_B1, * d_B2, * d_dB1, * d_dB2;
    cudaMalloc((void**)&d_B1, B1.size() * sizeof(double));
    cudaMalloc((void**)&d_B2, B2.size() * sizeof(double));
    cudaMalloc((void**)&d_dB1, dB1.size() * sizeof(double));
    cudaMalloc((void**)&d_dB2, dB2.size() * sizeof(double));
    cudaMemcpy(d_B1, B1.data(), B1.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_B2, B2.data(), B2.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_dB1, dB1.data(), dB1.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_dB2, dB2.data(), dB2.size() * sizeof(double), cudaMemcpyHostToDevice);

    double * d_NURBSExtraction1, * d_NURBSExtraction2;
    cudaMalloc((void**)&d_NURBSExtraction1, NURBSExtraction1.size() * sizeof(double));
    cudaMalloc((void**)&d_NURBSExtraction2, NURBSExtraction2.size() * sizeof(double));
    cudaMemcpy(d_NURBSExtraction1, NURBSExtraction1.data(),
        NURBSExtraction1.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_NURBSExtraction2, NURBSExtraction2.data(),
        NURBSExtraction2.size() * sizeof(double), cudaMemcpyHostToDevice);
    
    int * d_IEN, * d_ID;
    double * d_CP;
    cudaMalloc((void**)&d_IEN, IEN.size() * sizeof(int));
    cudaMalloc((void**)&d_ID, ID.size() * sizeof(int));
    cudaMalloc((void**)&d_CP, CP.size() * sizeof(double));
    cudaMemcpy(d_IEN, IEN.data(), IEN.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ID, ID.data(), ID.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_CP, CP.data(), CP.size() * sizeof(double), cudaMemcpyHostToDevice);

    double * d_elem_size1, * d_elem_size2;
    cudaMalloc((void**)&d_elem_size1, elem_size1.size() * sizeof(double));
    cudaMalloc((void**)&d_elem_size2, elem_size2.size() * sizeof(double));
    cudaMemcpy(d_elem_size1, elem_size1.data(), elem_size1.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_elem_size2, elem_size2.data(), elem_size2.size() * sizeof(double), cudaMemcpyHostToDevice);

    double * qw1, * qw2;
    cudaMalloc((void**)&qw1, nqp1 * sizeof(double));
    cudaMalloc((void**)&qw2, nqp2 * sizeof(double));
    cudaMemcpy(qw1, w1.data(), nqp1 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(qw2, w2.data(), nqp2 * sizeof(double), cudaMemcpyHostToDevice);

    double *d_F_array;
    VecCUDAGetArray(F, &d_F_array);

    AssembleKernel<<<dim3(nlocalelemx, nlocalelemy), dim3(nqp1, nqp2), 
        (nLocBas + 1) * sizeof(double)>>>(pp, qq,
        d_B1, d_B2, d_dB1, d_dB2,
        d_NURBSExtraction1, d_NURBSExtraction2,
        d_elem_size1, d_elem_size2,
        d_IEN, d_ID, d_CP,
        qw1, qw2, d_F_array);
    
    DirichletBCKernel<<<(Dir.size() + 255) / 256, 256>>>(Dir.data(), Dir.size(), d_F_array);

    VecCUDARestoreArray(F, &d_F_array);

    cudaFree(d_B1);
    cudaFree(d_B2);
    cudaFree(d_dB1);
    cudaFree(d_dB2);

    cudaFree(d_NURBSExtraction1);
    cudaFree(d_NURBSExtraction2);

    cudaFree(d_IEN);
    cudaFree(d_ID);
    cudaFree(d_CP);

    cudaFree(d_elem_size1);
    cudaFree(d_elem_size2);

    cudaFree(qw1);
    cudaFree(qw2);
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

void GlobalAssemblyMF::MatMulMF(QuadraturePoint * const &quad1,
    QuadraturePoint * const &quad2,
    const std::vector<int> &IEN,
    const std::vector<int> &ID,
    const std::vector<int> &Dir,
    const std::vector<double> &CP,
    const std::vector<double> &NURBSExtraction1,
    const std::vector<double> &NURBSExtraction2,
    const std::vector<double> &elem_size1,
    const std::vector<double> &elem_size2,
    ElementMF * const &elemmf,
    BernsteinBasis * const &bernstein,
    Vec x, Vec y)
{
    const int nqp1 = quad1->GetNumQuadraturePoint();
    const int nqp2 = quad2->GetNumQuadraturePoint();
    const std::vector<double> qp1 = quad1->GetQuadraturePoint();
    const std::vector<double> qp2 = quad2->GetQuadraturePoint();
    const std::vector<double> w1 = quad1->GetWeight();
    const std::vector<double> w2 = quad2->GetWeight();
    const int nLocBas = elemmf->GetNumLocalBasis();

    std::vector<double> B1{};
    std::vector<double> B2{};
    std::vector<double> dB1{};
    std::vector<double> dB2{};

    for (int i = 0; i < nqp1; ++i)
    {
        B1.push_back(bernstein->GetBernsteinBasisSingleQP(qp1[i]));
        dB1.push_back(bernstein->GetBernsteinBasisDerivativeSingleQP(qp1[i]));
    }
    for (int j = 0; j < nqp2; ++j)
    {
        B2.push_back(bernstein->GetBernsteinBasisSingleQP(qp2[j]));
        dB2.push_back(bernstein->GetBernsteinBasisDerivativeSingleQP(qp2[j]));
    }

    double * d_B1, * d_B2, * d_dB1, * d_dB2;
    cudaMalloc((void**)&d_B1, B1.size() * sizeof(double));
    cudaMalloc((void**)&d_B2, B2.size() * sizeof(double));
    cudaMalloc((void**)&d_dB1, dB1.size() * sizeof(double));
    cudaMalloc((void**)&d_dB2, dB2.size() * sizeof(double));
    cudaMemcpy(d_B1, B1.data(), B1.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_B2, B2.data(), B2.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_dB1, dB1.data(), dB1.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_dB2, dB2.data(), dB2.size() * sizeof(double), cudaMemcpyHostToDevice);

    double * d_NURBSExtraction1, * d_NURBSExtraction2;
    cudaMalloc((void**)&d_NURBSExtraction1, NURBSExtraction1.size() * sizeof(double));
    cudaMalloc((void**)&d_NURBSExtraction2, NURBSExtraction2.size() * sizeof(double));
    cudaMemcpy(d_NURBSExtraction1, NURBSExtraction1.data(),
        NURBSExtraction1.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_NURBSExtraction2, NURBSExtraction2.data(),
        NURBSExtraction2.size() * sizeof(double), cudaMemcpyHostToDevice);
    
    int * d_IEN, * d_ID;
    double * d_CP;
    cudaMalloc((void**)&d_IEN, IEN.size() * sizeof(int));
    cudaMalloc((void**)&d_ID, ID.size() * sizeof(int));
    cudaMalloc((void**)&d_CP, CP.size() * sizeof(double));
    cudaMemcpy(d_IEN, IEN.data(), IEN.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ID, ID.data(), ID.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_CP, CP.data(), CP.size() * sizeof(double), cudaMemcpyHostToDevice);

    double * d_elem_size1, * d_elem_size2;
    cudaMalloc((void**)&d_elem_size1, elem_size1.size() * sizeof(double));
    cudaMalloc((void**)&d_elem_size2, elem_size2.size() * sizeof(double));
    cudaMemcpy(d_elem_size1, elem_size1.data(), elem_size1.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_elem_size2, elem_size2.data(), elem_size2.size() * sizeof(double), cudaMemcpyHostToDevice);

    double * qw1, * qw2;
    cudaMalloc((void**)&qw1, nqp1 * sizeof(double));
    cudaMalloc((void**)&qw2, nqp2 * sizeof(double));
    cudaMemcpy(qw1, w1.data(), nqp1 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(qw2, w2.data(), nqp2 * sizeof(double), cudaMemcpyHostToDevice);

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

    cudaFree(d_B1);
    cudaFree(d_B2);
    cudaFree(d_dB1);
    cudaFree(d_dB2);

    cudaFree(d_NURBSExtraction1);
    cudaFree(d_NURBSExtraction2);

    cudaFree(d_IEN);
    cudaFree(d_ID);
    cudaFree(d_CP);

    cudaFree(d_elem_size1);
    cudaFree(d_elem_size2);

    cudaFree(qw1);
    cudaFree(qw2);
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