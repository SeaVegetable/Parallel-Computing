#include <iostream>
#include <cstdlib>
#include <ctime>
#include <chrono>

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

double* matrixVectorMultiplication(double* matrix, double* vector, int N)
{
    double * result = new double[N];

    for (int i = 0; i < N; ++i) {
        result[i] = 0;
        for (int j = 0; j < N; ++j) {
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
    double * result = matrixVectorMultiplication(matrix, vector, N);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;

    std::cout << "Time: " << diff.count() << " s" << std::endl;

    delete[] matrix;
    delete[] vector;
    delete[] result;
}