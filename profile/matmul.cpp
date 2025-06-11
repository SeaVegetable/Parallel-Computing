#include <iostream>
#include <vector>
#include <cstdlib>
#include <omp.h>

const int N = 1024;

void init_matrix(std::vector<std::vector<double>>& mat, bool randomize = true) {
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            mat[i][j] = randomize ? static_cast<double>(rand()) / RAND_MAX : 0.0;
}

void matmul(const std::vector<std::vector<double>>& A,
             const std::vector<std::vector<double>>& B,
             std::vector<std::vector<double>>& C,
             bool use_openmp) {
    #pragma omp parallel for collapse(2) if(use_openmp)
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            double sum = 0.0;
            for (int k = 0; k < N; ++k)
                sum += A[i][k] * B[k][j];
            C[i][j] = sum;
        }
}

int main() {
    std::vector<std::vector<double>> A(N, std::vector<double>(N));
    std::vector<std::vector<double>> B(N, std::vector<double>(N));
    std::vector<std::vector<double>> C(N, std::vector<double>(N));

    init_matrix(A);
    init_matrix(B, true);

    bool use_openmp = true;
    double start = omp_get_wtime();
    matmul(A, B, C, use_openmp);
    double end = omp_get_wtime();

    std::cout << "Matrix multiplication completed in " << (end - start) << " seconds." << std::endl;
    std::cout << "Use OpenMP: " << (use_openmp ? "yes" : "no") << std::endl;

    return 0;
}
