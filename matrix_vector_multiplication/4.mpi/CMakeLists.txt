cmake_minimum_required(VERSION 3.10)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(matmul)

# Find MPI package
find_package(MPI REQUIRED)

# Include MPI headers
include_directories(${MPI_CXX_INCLUDE_PATH})

# Add the executable
add_executable(matmul main.cpp)

# Link MPI libraries
target_link_libraries(matmul MPI::MPI_CXX)