#include <mpi.h>
#include <iostream>

double* generateRandomMatrix(int row, int col)
{
    double * matrix = new double[row * col];
    std::srand(std::time(0));
    
    for (int i = 0; i < row*col; ++i) {
        matrix[i] = static_cast<double>(std::rand()) / RAND_MAX * 2;
    }
    
    return matrix;
}

double* generateRandomVector(int row)
{
    double * vector = new double[row];
    std::srand(std::time(0));
    
    for (int i = 0; i < row; ++i) {
        vector[i] = static_cast<double>(std::rand()) / RAND_MAX * 2;
    }
    
    return vector;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int total_row, total_col;
    if (rank == 0)
    {
        std::cout << "Enter the size of the matrix: ";
        std::cin >> total_row >> total_col;
    }

    MPI_Bcast(&total_row, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&total_col, 1, MPI_INT, 0, MPI_COMM_WORLD);

    double * matrix = nullptr;
    double * vector = nullptr;

    if (rank == 0)
    {
        matrix = generateRandomMatrix(total_row, total_col);
        vector = generateRandomVector(total_col);
    }
    else
    {
        matrix = new double[total_col];
        vector = new double[total_col];
    }

    double start_time = MPI_Wtime();

    if (rank == 0)
    {
        double * result = new double[total_row];
        MPI_Bcast(vector, total_col, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        int numsent = 0;

        for (int i = 1; i < size; ++i)
        {
            MPI_Send(matrix + numsent * total_col, total_col, MPI_DOUBLE, i, numsent, MPI_COMM_WORLD);
        }
        numsent++;

        for (int i = 0; i < total_row; ++i)
        {
            MPI_Status status;
            double ans = 0;
            MPI_Recv(&ans, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            int sender = status.MPI_SOURCE;
            int tag = status.MPI_TAG;
            result[tag] = ans;

            if (numsent < total_row)
            {
                MPI_Send(matrix + numsent * total_col, total_col, MPI_DOUBLE, sender, numsent, MPI_COMM_WORLD);
                numsent++;
            }
        }

        for (int i = 1; i< size; ++i)
        {
            MPI_Send(MPI_BOTTOM, 0, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
        delete[] result;
    }
    else
    {
        double * local_matrix = new double[total_col];
        while(true)
        {
            MPI_Status status;
            MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            int count;
            MPI_Get_count(&status, MPI_DOUBLE, &count);

            if (count == 0)
            {
                MPI_Recv(nullptr, 0, MPI_DOUBLE, 0, status.MPI_TAG, MPI_COMM_WORLD, &status);
                break;
            }

            int row = status.MPI_TAG;
            MPI_Recv(local_matrix, total_col, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            double local_result = 0;
            for (int j = 0; j < total_col; ++j)
            {
                local_result += local_matrix[j] * vector[j];
            }

            MPI_Send(&local_result, 1, MPI_DOUBLE, 0, row, MPI_COMM_WORLD);
        }
        delete[] local_matrix;
    }

    double end_time = MPI_Wtime();

    if (rank == 0)
    {
        std::cout << "Time: " << end_time - start_time << std::endl;
    }

    MPI_Finalize();

    delete[] matrix;
    delete[] vector;
    return 0;
}