cmake_minimum_required(VERSION 3.10)
project(matrix_vector_multiplication)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(CUDA REQUIRED)

add_executable(main main.cpp)
target_link_libraries(main ${CUDA_LIBRARIES})