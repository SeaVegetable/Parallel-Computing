#include "mult.cuh"

__device__ void compute_jacobian_basis(
    const int in_p, const int in_q, const double h1, const double h2,
    const double *d_B1, const double *d_B2,
    const double *d_dB1, const double *d_dB2,
    const double *s_nurbs_extraction1, const double *s_nurbs_extraction2,
    const double *eCP,
    double &jacobian,
    double *R)
{
    const int p = 3;
    const int q = 3;
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

    const int nLocBas = (p + 1) * (q + 1);
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
    jacobian *= h1*h2;
}

__device__ void compute_jacobian_derivative(
    const int in_p, const int in_q, const double h1, const double h2,
    const double *d_B1, const double *d_B2,
    const double *d_dB1, const double *d_dB2,
    const double *s_nurbs_extraction1, const double *s_nurbs_extraction2,
    const double *eCP, 
    double &jacobian,
    double *dR_dx, 
    double *dR_dy)
{
    const int p = 3;
    const int q = 3;
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
            N1[jj] +=  s_nurbs_extraction1[jj * (p + 1) + kk] * d_B1[kk];
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

    const int nLocBas = (p + 1) * (q + 1);
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

    double R[nLocBas];
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

__device__ double get_force(double x, double y)
{
    return x * (1.0 - x) * y * (1.0 - y);
}

__global__ void AssembleKernel(const int in_p, const int in_q,
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
    const int p = 3;
    const int q = 3;
    extern __shared__ char shared_data[];

    int offset = 0;
    int *s_eID = (int*)(shared_data + offset);
    offset += (p + 1) * (q + 1) * sizeof(int);
    double *s_eCP = (double*)(shared_data + offset);
    offset += 2 * (p + 1) * (q + 1) * sizeof(double);
    double *s_eNURBSExtraction1 = (double*)(shared_data + offset);
    offset += (p + 1) * (p + 1) * sizeof(double);
    double *s_eNURBSExtraction2 = (double*)(shared_data + offset);
    offset += (q + 1) * (q + 1) * sizeof(double);
    double *s_qw = (double*)(shared_data + offset);

    int elemIndex = blockIdx.y * gridDim.x + blockIdx.x;
    const int nLocBas = (p + 1) * (q + 1);

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
        s_eNURBSExtraction1[i] = d_nurbs_extraction1[blockIdx.x * (p + 1) * (p + 1) + i];
    for (int i = 0; i < (q + 1) * (q + 1); ++i)
        s_eNURBSExtraction2[i] = d_nurbs_extraction2[blockIdx.y * (q + 1) * (q + 1) + i];

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

    if (qp < (p+1)*(q+1))
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

        double force = get_force(x, y);

        for (int ii = 0; ii < nLocBas; ++ii)
        {
            double val = R[ii] * force * jacobian * s_qw[qp];
            int coo_index = s_eID[ii];
            if (coo_index >= 0)
            {
                atomicAdd(&d_x_array[coo_index], val);
            }
        }
    }
}

__global__ void MatrixFreeMatMultKernel(const int in_p, const int in_q,
    double *d_B1, double *d_B2,
    double *d_dB1, double *d_dB2,
    double *d_nurbs_extraction1, double *d_nurbs_extraction2,
    double *d_elem_size1, double *d_elem_size2,
    int *d_IEN, int *d_ID,
    double *d_CP,
    double *qw1, double *qw2,
    const double *d_F_array_in,
    double *d_F_array_out
    )
{
    extern __shared__ char shared_data[];

    const int p = 3;
    const int q = 3;
    const int nLocBas = (p + 1) * (q + 1);

    int offset = 0;
    int *s_eID = (int*)(shared_data + offset);
    offset += nLocBas * sizeof(int);
    double *s_eCP = (double*)(shared_data + offset);
    offset += 2 * nLocBas * sizeof(double);
    double *s_eNURBSExtraction1 = (double*)(shared_data + offset);
    offset += (p + 1) * (p + 1) * sizeof(double);
    double *s_eNURBSExtraction2 = (double*)(shared_data + offset);
    offset += (q + 1) * (q + 1) * sizeof(double);
    double *s_qw = (double*)(shared_data + offset);
    offset += (p + 1) * (q + 1) * sizeof(double);
    double *Floc_in = (double*)(shared_data + offset);
    offset += nLocBas * sizeof(double);
    double *Floc_out = (double*)(shared_data + offset);
    
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
        s_eNURBSExtraction1[i] = d_nurbs_extraction1[blockIdx.x * (p + 1) * (p + 1) + i];
    for (int i = 0; i < (q + 1) * (q + 1); ++i)
        s_eNURBSExtraction2[i] = d_nurbs_extraction2[blockIdx.y * (q + 1) * (q + 1) + i];

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

    if (qp < (p+1)*(q+1))
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
            temp_x += dR_dx[jj] * Floc_in[jj];
            temp_y += dR_dy[jj] * Floc_in[jj];
        }

        temp_x *= -s_qw[qp]*jacobian;
        temp_y *= -s_qw[qp]*jacobian;

        for (int ii = 0; ii < nLocBas; ++ii)
        {
            Floc_out[ii] += (dR_dx[ii] * temp_x + dR_dy[ii] * temp_y);
        }

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

__global__ void DirichletBCKernel(const int * d_Dir, const int dirsize, double * d_val, double value)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < dirsize)
    {
        int coo_index = d_Dir[idx];
        if (coo_index >= 0)
        {
            d_val[coo_index] = value;
        }
    }
}

void AssembleLoadCUDA(const int p, const int q,
    const int nlocalelemx, const int nlocalelemy,
    double * d_B1, double * d_B2,
    double * d_dB1, double * d_dB2,
    double * d_nurbs_extraction1, double * d_nurbs_extraction2,
    double * d_elem_size1, double * d_elem_size2,
    int * d_IEN, int * d_ID, double * d_CP,
    double * qw1, double * qw2, double * d_F_array)
{
    int shared_size = (p + 1) * (q + 1) * sizeof(int)
                + 2 * (p + 1) * (q + 1) * sizeof(double)
                + (p + 1) * (p + 1) * sizeof(double)
                + (q + 1) * (q + 1) * sizeof(double)
                + (p + 1) * (q + 1) * sizeof(double);

    AssembleKernel<<<dim3(nlocalelemx, nlocalelemy), dim3(p+1, q+1), shared_size>>>(
        p, q, d_B1, d_B2, d_dB1, d_dB2,
        d_nurbs_extraction1, d_nurbs_extraction2,
        d_elem_size1, d_elem_size2,
        d_IEN, d_ID, d_CP,
        qw1, qw2, d_F_array);

    cudaDeviceSynchronize();

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        printf("Error in AssembleLoadCUDA: %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

void MatrixFreeMatMultCUDA(const int p, const int q,
    const int nlocalelemx, const int nlocalelemy,
    double * d_B1, double * d_B2,
    double * d_dB1, double * d_dB2,
    double * d_nurbs_extraction1, double * d_nurbs_extraction2,
    double * d_elem_size1, double * d_elem_size2,
    int * d_IEN, int * d_ID, double * d_CP,
    double * qw1, double * qw2,
    const double * d_F_array_in, double * d_F_array_out)
{
    int shared_size = (p + 1) * (q + 1) * sizeof(int)
                + 2 * (p + 1) * (q + 1) * sizeof(double)
                + (p + 1) * (p + 1) * sizeof(double)
                + (q + 1) * (q + 1) * sizeof(double)
                + (p + 1) * (q + 1) * sizeof(double)
                + (p + 1) * (q + 1) * sizeof(double)
                + (p + 1) * (q + 1) * sizeof(double);

    MatrixFreeMatMultKernel<<<dim3(nlocalelemx, nlocalelemy), dim3(p+1, q+1), shared_size>>>(
        p, q, d_B1, d_B2, d_dB1, d_dB2,
        d_nurbs_extraction1, d_nurbs_extraction2,
        d_elem_size1, d_elem_size2,
        d_IEN, d_ID, d_CP,
        qw1, qw2,
        d_F_array_in, d_F_array_out);

    cudaDeviceSynchronize();

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        printf("Error in MatrixFreeMatMultCUDA: %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

void DirichletBCCUDA(const int * d_Dir, const int dirsize, double * d_x_array, double value)
{
    int blocksize = 256;
    int nblocks = (dirsize + blocksize - 1) / blocksize;

    DirichletBCKernel<<<nblocks, blocksize>>>(d_Dir, dirsize, d_x_array, value);

    cudaDeviceSynchronize();

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        printf("Error in DirichletBCCUDA: %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}