#include <iostream>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <arm_neon.h>

double* generateRandomMatrix(int N)
{
    double * matrix = new double[N * N];
    std::srand(std::time(0));
    
    for (int i = 0; i < N * N; ++i) {
        matrix[i] = static_cast<double>(std::rand()) / RAND_MAX * 2;
    }
    
    return matrix;
}

double* generateRandomVector(int N)
{
    double * vector = new double[N];
    std::srand(std::time(0));
    
    for (int i = 0; i < N; ++i) {
        vector[i] = static_cast<double>(std::rand()) / RAND_MAX * 2;
    }
    
    return vector;
}

double* matrixVectorMultiplicationNEON(double* matrix, double* vector, int N)
{
    double * result = new double[N];

    for (int i = 0; i < N; ++i)
    {
        result[i] = 0;
        float64x2_t sum = vdupq_n_f64(0.0);

        int j = 0;
        for (; j < N; j += 4)
        {
            float64x2_t matrix_0 = vld1q_f64(matrix + i * N + j);
            float64x2_t vector_0 = vld1q_f64(vector + j);
            sum = vmlaq_f64(sum, matrix_0, vector_0);

            float64x2_t matrix_1 = vld1q_f64(matrix + i * N + j + 2);
            float64x2_t vector_1 = vld1q_f64(vector + j + 2);
            sum = vmlaq_f64(sum, matrix_1, vector_1);
        }

        result[i] = vgetq_lane_f64(sum, 0) + vgetq_lane_f64(sum, 1);

        for (; j < N; ++j)
        {
            result[i] += matrix[i * N + j] * vector[j];
        }
    }
    
    return result;
}

int main()
{
    int N;
    std::cout << "Enter the size of the matrix: ";
    std::cin >> N;

    auto start = std::chrono::high_resolution_clock::now();

    double * matrix = generateRandomMatrix(N);
    double * vector = generateRandomVector(N);
    double * result = matrixVectorMultiplicationNEON(matrix, vector, N);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;

    std::cout << "Time: " << diff.count() << " s" << std::endl;

    delete[] matrix;
    delete[] vector;
    delete[] result;
}