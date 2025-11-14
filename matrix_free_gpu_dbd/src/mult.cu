#ifndef P_ORDER
#define P_ORDER 3
#endif
#ifndef Q_ORDER
#define Q_ORDER 3
#endif
#define NLOC_MOD2 (((P_ORDER + 1) * (Q_ORDER + 1)) % 2)
#define MYOFFSET (NLOC_MOD2 ? 4 : 0)

#include "mult.cuh"

__device__ void compute_jacobian_basis(
    const double h1, const double h2,
    const double *d_B1, const double *d_B2,
    const double *d_dB1, const double *d_dB2,
    const double *s_nurbs_extraction1, const double *s_nurbs_extraction2,
    const double *eCP,
    double &jacobian,
    double *R)
{
    const int nx = P_ORDER + 1;
    const int ny = Q_ORDER + 1;
    double N1[nx];
    double dN1[nx];
    for (int i = 0; i < nx; ++i)
    {
        N1[i] = 0.0;
        dN1[i] = 0.0;
    }
    double N2[ny];
    double dN2[ny];
    for (int j = 0; j < ny; ++j)
    {
        N2[j] = 0.0;
        dN2[j] = 0.0;
    }

    for (int jj = 0; jj < nx; ++jj)
    {
        for (int kk = 0; kk < nx; ++kk)
        {
            N1[jj] += s_nurbs_extraction1[jj * (nx) + kk] * d_B1[kk];
            dN1[jj] += s_nurbs_extraction1[jj * (nx) + kk] * d_dB1[kk];
        }
        dN1[jj] /= h1;
    }
    for (int jj = 0; jj < ny; ++jj)
    {
        for (int kk = 0; kk < ny; ++kk)
        {
            N2[jj] += s_nurbs_extraction2[jj * (ny) + kk] * d_B2[kk];
            dN2[jj] += s_nurbs_extraction2[jj * (ny) + kk] * d_dB2[kk];
        }
        dN2[jj] /= h2;
    }

    const int nLocBas = nx * ny;
    double N[nLocBas];
    double dN_dxi[nLocBas];
    double dN_deta[nLocBas];
    double w = 0.0;
    double dw_dxi = 0.0;
    double dw_deta = 0.0;
    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            N[j * (nx) + i] = N1[i] * N2[j];
            w += N[j * (nx) + i];
            dN_dxi[j * (nx) + i] = dN1[i] * N2[j];
            dw_dxi += dN_dxi[j * (nx) + i];
            dN_deta[j * (nx) + i] = N1[i] * dN2[j];
            dw_deta += dN_deta[j * (nx) + i];
        }
    }

    double dR_dxi[nLocBas];
    double dR_deta[nLocBas];
    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            R[j * (nx) + i] = N[j * (nx) + i] / w;
            dR_dxi[j * (nx) + i] = (dN_dxi[j * (nx) + i] - dw_dxi * R[j * (nx) + i]) / w;
            dR_deta[j * (nx) + i] = (dN_deta[j * (nx) + i] - dw_deta * R[j * (nx) + i]) / w;
        }
    }

    double dx_dxi = 0.0;
    double dx_deta = 0.0;
    double dy_dxi = 0.0;
    double dy_deta = 0.0;

    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            dx_dxi += eCP[2 * (j * (nx) + i)] * dR_dxi[j * (nx) + i];
            dx_deta += eCP[2 * (j * (nx) + i)] * dR_deta[j * (nx) + i];
            dy_dxi += eCP[2 * (j * (nx) + i) + 1] * dR_dxi[j * (nx) + i];
            dy_deta += eCP[2 * (j * (nx) + i) + 1] * dR_deta[j * (nx) + i];
        }
    }
    jacobian = dx_dxi * dy_deta - dx_deta * dy_dxi;
    jacobian *= h1*h2;
}

__device__ void compute_jacobian_derivative(
    const double h1, const double h2,
    const double *d_B1, const double *d_B2,
    const double *d_dB1, const double *d_dB2,
    const double *s_nurbs_extraction1, const double *s_nurbs_extraction2,
    const double *eCP, 
    double &jacobian,
    double *dR_dx, 
    double *dR_dy)
{
    const int nx = P_ORDER + 1;
    const int ny = Q_ORDER + 1;
    double N1[nx];
    double dN1[nx];
    for (int i = 0; i < nx; ++i)
    {
        N1[i] = 0.0;
        dN1[i] = 0.0;
    }
    double N2[ny];
    double dN2[ny];
    for (int j = 0; j < ny; ++j)
    {
        N2[j] = 0.0;
        dN2[j] = 0.0;
    }

    for (int jj = 0; jj < nx; ++jj)
    {
        for (int kk = 0; kk < nx; ++kk)
        {
            N1[jj] +=  s_nurbs_extraction1[jj * nx + kk] * d_B1[kk];
            dN1[jj] += s_nurbs_extraction1[jj * nx + kk] * d_dB1[kk];
        }
        dN1[jj] /= h1;
    }
    for (int jj = 0; jj < ny; ++jj)
    {
        for (int kk = 0; kk < ny; ++kk)
        {
            N2[jj] += s_nurbs_extraction2[jj * ny + kk] * d_B2[kk];
            dN2[jj] += s_nurbs_extraction2[jj * ny + kk] * d_dB2[kk];
        }
        dN2[jj] /= h2;
    }

    const int nLocBas = nx * ny;
    double N[nLocBas];
    double dN_dxi[nLocBas];
    double dN_deta[nLocBas];
    double w = 0.0;
    double dw_dxi = 0.0;
    double dw_deta = 0.0;

    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            N[j * nx + i] = N1[i] * N2[j];
            w += N[j * nx + i];
            dN_dxi[j * nx + i] = dN1[i] * N2[j];
            dw_dxi += dN_dxi[j * nx + i];
            dN_deta[j * nx + i] = N1[i] * dN2[j];
            dw_deta += dN_deta[j * nx + i];
        }
    }

    double R[nLocBas];
    double dR_dxi[nLocBas];
    double dR_deta[nLocBas];

    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            R[j * nx + i] = N[j * nx + i] / w;
            dR_dxi[j * nx + i] = (dN_dxi[j * nx + i] - dw_dxi * R[j * nx + i]) / w;
            dR_deta[j * nx + i] = (dN_deta[j * nx + i] - dw_deta * R[j * nx + i]) / w;
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

    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            dx_dxi += eCP[2 * (j * nx + i)] * dR_dxi[j * nx + i];
            dx_deta += eCP[2 * (j * nx + i)] * dR_deta[j * nx + i];
            dy_dxi += eCP[2 * (j * nx + i) + 1] * dR_dxi[j * nx + i];
            dy_deta += eCP[2 * (j * nx + i) + 1] * dR_deta[j * nx + i];
        }
    }

    jacobian = dx_dxi * dy_deta - dx_deta * dy_dxi;

    dxi_dx = dy_deta / jacobian;
    dxi_dy = -dx_deta / jacobian;
    deta_dx = -dy_dxi / jacobian;
    deta_dy = dx_dxi / jacobian;

    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            dR_dx[j * nx + i] = dxi_dx * dR_dxi[j * nx + i] + deta_dx * dR_deta[j * nx + i];
            dR_dy[j * nx + i] = dxi_dy * dR_dxi[j * nx + i] + deta_dy * dR_deta[j * nx + i];
        }
    }

    jacobian *= h1 * h2;
}

__device__ double get_force(double x, double y)
{
    return x * (1.0 - x) * y * (1.0 - y);
}

__global__ void AssembleKernel(
    int nfunc,
    double *d_B1, double *d_B2,
    double *d_dB1, double *d_dB2,
    double *d_nurbs_extraction1, double *d_nurbs_extraction2,
    double *d_elem_size1, double *d_elem_size2,
    int *d_IEN, int *d_ID,
    double *d_CP,
    double *qw1, double *qw2,
    int *d_invlm_elemNum, int *d_invlm_offset,
    int *d_invlm_elemIdx, int *d_invlm_baseIdx,
    int *d_xelemIdx, int *d_yelemIdx,
    double *d_x_array
    )
{
    const int nx = P_ORDER + 1;
    const int ny = Q_ORDER + 1;
    const int nLocBas = nx * ny;

    double s_qw[nLocBas];

    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            s_qw[j * nx + i] = qw1[i] * qw2[j];
        }
    }

    int nodeIndex = blockIdx.x * blockDim.x + threadIdx.x;

    if (nodeIndex < nfunc)
    {
        int invlm_elemNum = d_invlm_elemNum[nodeIndex];
        int invlm_offset = d_invlm_offset[nodeIndex];

        for (int ee = 0; ee < invlm_elemNum; ++ee)
        {
            int elemIndex = d_invlm_elemIdx[invlm_offset + ee];
            int baseIndex = d_invlm_baseIdx[invlm_offset + ee];

            int eIEN[nLocBas];
            double eCP[2 * nLocBas];
            for (int j = 0; j < nLocBas; ++j)
            {
                eIEN[j] = d_IEN[elemIndex * nLocBas + j];
                eCP[2 * j] = d_CP[2 * d_IEN[elemIndex * nLocBas + j]];
                eCP[2 * j + 1] = d_CP[2 * d_IEN[elemIndex * nLocBas + j] + 1];
            }

            double eNURBSExtraction1[nx * nx];
            double eNURBSExtraction2[ny * ny];
            for (int i = 0; i < nx * nx; ++i)
                eNURBSExtraction1[i] = d_nurbs_extraction1[d_xelemIdx[elemIndex] * nx * nx + i];
            for (int i = 0; i < ny * ny; ++i)
                eNURBSExtraction2[i] = d_nurbs_extraction2[d_yelemIdx[elemIndex] * ny * ny + i];

            double h1 = d_elem_size1[d_xelemIdx[elemIndex]];
            double h2 = d_elem_size2[d_yelemIdx[elemIndex]];

            int coo_index = eIEN[baseIndex];

            for (int qpy = 0; qpy < ny; ++qpy)
            {
                for (int qpx = 0; qpx < nx; ++qpx)
                {
                    double B1[nx];
                    double dB1[nx];
                    double B2[ny];
                    double dB2[ny];

                    for (int i = 0; i < nx; ++i)
                    {
                        B1[i] = d_B1[qpx * nx + i];
                        dB1[i] = d_dB1[qpx * nx + i];
                    }
                    for (int j = 0; j < ny; ++j)
                    {
                        B2[j] = d_B2[qpy * ny + j];
                        dB2[j] = d_dB2[qpy * ny + j];
                    }

                    double jacobian;
                    double R[nLocBas];

                    compute_jacobian_basis(h1, h2, B1, B2, dB1, dB2,
                        eNURBSExtraction1, eNURBSExtraction2, eCP, jacobian, R);
                    
                    double x = 0.0;
                    double y = 0.0;

                    for (int ii = 0; ii < nLocBas; ++ii)
                    {
                        x += eCP[2 * ii] * R[ii];
                        y += eCP[2 * ii + 1] * R[ii];
                    }

                    double force = get_force(x, y);

                    double val = R[baseIndex] * force * jacobian * s_qw[qpy * nx + qpx];

                    d_x_array[coo_index] += val;
                }
            }
        }
    }
}

__global__ void MatrixFreeMatMultKernel(
    int nfunc,
    double *d_B1, double *d_B2,
    double *d_dB1, double *d_dB2,
    double *d_nurbs_extraction1, double *d_nurbs_extraction2,
    double *d_elem_size1, double *d_elem_size2,
    int *d_IEN, int *d_ID,
    double *d_CP,
    double *qw1, double *qw2,
    int *d_invlm_elemNum, int *d_invlm_offset,
    int *d_invlm_elemIdx, int *d_invlm_baseIdx,
    int *d_xelemIdx, int *d_yelemIdx,
    const double *d_F_array_in,
    double *d_F_array_out
    )
{
    const int nx = P_ORDER + 1;
    const int ny = Q_ORDER + 1;
    const int nLocBas = nx * ny;

    double s_qw[nx * ny];

    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            s_qw[j * nx + i] = qw1[i] * qw2[j];
        }
    }
    
    int nodeIndex = blockIdx.x * blockDim.x + threadIdx.x;

    if (nodeIndex < nfunc)
    {
        int invlm_elemNum = d_invlm_elemNum[nodeIndex];
        int invlm_offset = d_invlm_offset[nodeIndex];

        for (int ee = 0; ee < invlm_elemNum; ++ee)
        {
            int elemIndex = d_invlm_elemIdx[invlm_offset + ee];
            int baseIndex = d_invlm_baseIdx[invlm_offset + ee];

            int eIEN[nLocBas];
            double eCP[2 * nLocBas];
            for (int j = 0; j < nLocBas; ++j)
            {
                eIEN[j] = d_IEN[elemIndex * nLocBas + j];
                eCP[2 * j] = d_CP[2 * d_IEN[elemIndex * nLocBas + j]];
                eCP[2 * j + 1] = d_CP[2 * d_IEN[elemIndex * nLocBas + j] + 1];
            }

            double eNURBSExtraction1[nx * nx];
            double eNURBSExtraction2[ny * ny];
            for (int i = 0; i < nx * nx; ++i)
                eNURBSExtraction1[i] = d_nurbs_extraction1[d_xelemIdx[elemIndex] * nx * nx + i];
            for (int i = 0; i < ny * ny; ++i)
                eNURBSExtraction2[i] = d_nurbs_extraction2[d_yelemIdx[elemIndex] * ny * ny + i];

            double h1 = d_elem_size1[d_xelemIdx[elemIndex]];
            double h2 = d_elem_size2[d_yelemIdx[elemIndex]];

            int coo_index = eIEN[baseIndex];

            for (int qpy = 0; qpy < ny; ++qpy)
            {
                for (int qpx = 0; qpx < nx; ++qpx)
                {
                    double B1[nx];
                    double dB1[nx];
                    double B2[ny];
                    double dB2[ny];

                    for (int i = 0; i < nx; ++i)
                    {
                        B1[i] = d_B1[qpx * nx + i];
                        dB1[i] = d_dB1[qpx * nx + i];
                    }
                    for (int j = 0; j < ny; ++j)
                    {
                        B2[j] = d_B2[qpy * ny + j];
                        dB2[j] = d_dB2[qpy * ny + j];
                    }

                    double jacobian;
                    double dR_dx[nLocBas];
                    double dR_dy[nLocBas];

                    compute_jacobian_derivative(h1, h2, B1, B2, dB1, dB2,
                        eNURBSExtraction1, eNURBSExtraction2, eCP, jacobian, dR_dx, dR_dy);

                    double temp_x = 0.0;
                    double temp_y = 0.0;

                    for (int jj = 0; jj < nLocBas; ++jj)
                    {
                        temp_x += dR_dx[jj] * d_F_array_in[eIEN[jj]];
                        temp_y += dR_dy[jj] * d_F_array_in[eIEN[jj]];
                    }

                    temp_x *= -s_qw[qpy * nx + qpx] * jacobian;
                    temp_y *= -s_qw[qpy * nx + qpx] * jacobian;

                    double val = (dR_dx[baseIndex] * temp_x + dR_dy[baseIndex] * temp_y);

                    d_F_array_out[coo_index] += val;
                }
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

void AssembleLoadCUDA(
    const int nfunc,
    const int p, const int q,
    const int nlocalelemx, const int nlocalelemy,
    double * d_B1, double * d_B2,
    double * d_dB1, double * d_dB2,
    double * d_nurbs_extraction1, double * d_nurbs_extraction2,
    double * d_elem_size1, double * d_elem_size2,
    int * d_IEN, int * d_ID, double * d_CP,
    double * qw1, double * qw2,
    int * d_invlm_elemNum, int * d_invlm_offset,
    int * d_invlm_elemIdx, int * d_invlm_baseIdx,
    int * d_xelemIdx, int * d_yelemIdx,
    double * d_F_array)
{
    if (p != P_ORDER || q != Q_ORDER)
    {
        printf("Error: p and q must match the defined P_ORDER and Q_ORDER.\n");
        exit(EXIT_FAILURE);
    }

    int blocksize = 256;
    int nblocks = (nfunc + blocksize - 1) / blocksize;

    AssembleKernel<<<nblocks, blocksize>>>(
        nfunc,
        d_B1, d_B2, d_dB1, d_dB2,
        d_nurbs_extraction1, d_nurbs_extraction2,
        d_elem_size1, d_elem_size2,
        d_IEN, d_ID, d_CP,
        qw1, qw2,
        d_invlm_elemNum, d_invlm_offset,
        d_invlm_elemIdx, d_invlm_baseIdx,
        d_xelemIdx, d_yelemIdx,
        d_F_array);

    cudaDeviceSynchronize();

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        printf("Error in AssembleLoadCUDA: %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

void MatrixFreeMatMultCUDA(
    const int nfunc,
    const int p, const int q,
    const int nlocalelemx, const int nlocalelemy,
    double * d_B1, double * d_B2,
    double * d_dB1, double * d_dB2,
    double * d_nurbs_extraction1, double * d_nurbs_extraction2,
    double * d_elem_size1, double * d_elem_size2,
    int * d_IEN, int * d_ID, double * d_CP,
    double * qw1, double * qw2,
    int * d_invlm_elemNum, int * d_invlm_offset,
    int * d_invlm_elemIdx, int * d_invlm_baseIdx,
    int * d_xelemIdx, int * d_yelemIdx,
    const double * d_F_array_in, double * d_F_array_out)
{
    if (p != P_ORDER || q != Q_ORDER)
    {
        printf("Error: p and q must match the defined P_ORDER and Q_ORDER.\n");
        exit(EXIT_FAILURE);
    }

    int blocksize = 256;
    int nblocks = (nfunc + blocksize - 1) / blocksize;

    MatrixFreeMatMultKernel<<<nblocks, blocksize>>>(
        nfunc,
        d_B1, d_B2, d_dB1, d_dB2,
        d_nurbs_extraction1, d_nurbs_extraction2,
        d_elem_size1, d_elem_size2,
        d_IEN, d_ID, d_CP,
        qw1, qw2,
        d_invlm_elemNum, d_invlm_offset,
        d_invlm_elemIdx, d_invlm_baseIdx,
        d_xelemIdx, d_yelemIdx,
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