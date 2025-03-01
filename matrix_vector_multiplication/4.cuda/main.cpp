#include <iostream>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <cuda_runtime.h>

__global__ void matrixVectorKernel(const double* matrix, const double* vec, double* result, int N) {
    int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row < N) {
        double sum = 0.0;
        for (int j = 0; j < N; ++j) {
            sum += matrix[row * N + j] * vec[j];
        }
        result[row] = sum;
    }
}

void matrixVectorMultiplyCUDA(double* matrix, double* vec, double* result, int N) {
    double *d_matrix = nullptr, *d_vec = nullptr, *d_result = nullptr;

    size_t sizeMatrix = sizeof(double) * N * N;
    size_t sizeVec = sizeof(double) * N;

    cudaMalloc((void**)&d_matrix, sizeMatrix);
    cudaMalloc((void**)&d_vec, sizeVec);
    cudaMalloc((void**)&d_result, sizeVec);
    
    cudaMemcpy(d_matrix, matrix, sizeMatrix, cudaMemcpyHostToDevice);
    cudaMemcpy(d_vec, vec, sizeVec, cudaMemcpyHostToDevice);
    cudaMemset(d_result, 0, sizeVec);

    int blockSize = 256;
    int gridSize = (N + blockSize - 1) / blockSize;
    matrixVectorKernel<<<gridSize, blockSize>>>(d_matrix, d_vec, d_result, N);

    cudaDeviceSynchronize();

    cudaMemcpy(result, d_result, sizeVec, cudaMemcpyDeviceToHost);

    cudaFree(d_matrix);
    cudaFree(d_vec);
    cudaFree(d_result);
}

double* generateRandomMatrix(int N)
{
    double* matrix = new double[N * N];
    std::srand(std::time(0));
    for (int i = 0; i < N * N; ++i) {
        matrix[i] = static_cast<double>(std::rand()) / RAND_MAX * 2.0;
    }
    return matrix;
}

double* generateRandomVector(int N)
{
    double* vec = new double[N];
    for (int i = 0; i < N; ++i) {
        vec[i] = static_cast<double>(std::rand()) / RAND_MAX * 2.0;
    }
    return vec;
}


int main()
{
    int N;
    std::cout << "Enter the size of the matrix: ";
    std::cin >> N;

    auto start = std::chrono::high_resolution_clock::now();

    double* matrix = generateRandomMatrix(N);
    double* vector = generateRandomVector(N);
    double* result = new double[N];

    matrixVectorMultiplyCUDA(matrix, vec, result, N);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    delete[] matrix;
    delete[] vector;
    delete[] result;

    return 0;
}
