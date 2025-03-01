#include <iostream>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <thread>
#include <vector>

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

void matrixVectorMultiplication(double* matrix, double* vector, double* result, int N, int start, int end)
{
    for (int i = start; i < end; ++i) {
        result[i] = 0;
        for (int j = 0; j < N; ++j) {
            result[i] += matrix[i * N + j] * vector[j];
        }
    }
}

int main()
{
    int N;
    std::cout << "Enter the size of the matrix: ";
    std::cin >> N;

    auto start = std::chrono::high_resolution_clock::now();

    double * matrix = generateRandomMatrix(N);
    double * vector = generateRandomVector(N);
    double * result = new double[N];

    int numThreads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads;
    int chunkSize = N / numThreads;

    for (int i = 0; i < numThreads; ++i) {
        int start = i * chunkSize;
        int end = (i == numThreads - 1) ? N : (i + 1) * chunkSize;
        threads.push_back(std::thread(matrixVectorMultiplication, matrix, vector, result, N, start, end));
    }

    for (auto& thread : threads) {
        thread.join();
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;

    std::cout << "Time: " << diff.count() << " s" << std::endl;

    delete[] matrix;
    delete[] vector;
    delete[] result;
}