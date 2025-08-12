#include <coda_runtime.h>
#include "assembly.hpp"

void AssemStiffness(const int nLocBas, const int nqp1, const int nqp2,
    const int nlocalelemx, const int nlocalelemy,
    const double * d_N, const double * d_dN_dxi, const double * d_dN_deta,
    const double * d_weight, const int * d_IEN, const int * d_ID,
    const int * d_Dir, const double * d_CP, const int * d_elem2coo,
    double * d_val)
{
    int nelem = nlocalelemx * nlocalelemy;
    int blocksize = 256;
    int gridsize = (nelem + blocksize - 1) / blocksize;

    AssembleStiffnessKernel<<<gridsize, blocksize>>>(nLocBas, nqp1*nqp2,
        d_N, d_dN_dxi, d_dN_deta,
        d_weight, d_IEN, d_ID,
        d_Dir, d_CP, d_elem2coo,
        d_val);

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "Error launching kernel: %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

void DirichletBCK(const int * d_dir2coo, const int dir_size, double * d_val)
{
    int blockSize = 256;
    int numBlocks = (dir_size + blockSize - 1) / blockSize;

    DirichletBCKKernel<<<numBlocks, blockSize>>>(d_dir2coo, dir_size, d_val);

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "Error launching DirichletBC kernel: %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

__global__ void AssembleKernel(const int nLocBas, const int nqp,
    const double * d_N, const double * d_dN_dxi, const double * d_dN_deta,
    const double * d_weight,
    const int * d_IEN, const int * d_ID,
    const double * d_CP, const int * elem2coo,
    double * d_val)
{
    int elemIndex = blockIdx.x * blockDim.x + threadIdx.x;

    if (elemIndex >= gridDim.x * blockDim.x) return;

    int eID[nLocBas];
    double eCP[2 * nLocBas];

    for (int j = 0; j < nLocBas; ++j)
    {
        int ien = d_IEN[elemIndex * nLocBas + j];
        eID[j] = d_ID[ien];
        eCP[2 * j]     = d_CP[2 * ien];
        eCP[2 * j + 1] = d_CP[2 * ien + 1];
    }

    for (int qp = 0; qp < nqp; ++qp)
    {
        double R[nLocBas];
        double dN_dxi_q[nLocBas];
        double dN_deta_q[nLocBas];

        for (int i = 0; i < nLocBas; ++i)
        {
            R[i] = d_N[qp * nLocBas + i];
            dN_dxi_q[i] = d_dN_dxi[qp * nLocBas + i];
            dN_deta_q[i] = d_dN_deta[qp * nLocBas + i];
        }

        double jacobian;
        double dR_dx[nLocBas];
        double dR_dy[nLocBas];
        compute_jacobian_basis_derivative(nLocBas, dN_dxi_q, dN_deta_q, eCP, jacobian, dR_dx, dR_dy);

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

    if (idx >= dirsize) return;

    int coo_index = dir2coo[idx];
    d_val[coo_index] = 1.0;
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