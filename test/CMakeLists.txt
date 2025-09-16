cmake_minimum_required(VERSION 3.18)

project(EigenvalueProblem LANGUAGES CXX C CUDA)

enable_language(CUDA)

IF(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
  cmake_policy(SET CMP0074 NEW)
ENDIF(COMMAND cmake_policy)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Add your source files here
set(CUDA_INC "/public/software/compiler/nvidia/cuda/12.1.0/include")

# Include directories
include_directories(${EIGEN3_INCLUDE_DIR})

# Find PETSc
set(PETSC_DIR "/work/home/mae_yuexh/lib/petsc-3.20.6-opt")
set(PETSC_ARCH .)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/module")
find_package(PETSc REQUIRED)
message(STATUS "PETSC_INC: ${PETSC_INC}")
message(STATUS "PETSC_LIB: ${PETSC_LIB}")

set(CMAKE_C_COMPILER /public/software/mpi/openmpi/4.1.5/bin/mpicc)
set(CMAKE_CXX_COMPILER /public/software/mpi/openmpi/4.1.5/bin/mpicxx)

include_directories(${PETSC_INC})

add_executable(test ghost_with_cuda.cpp)
target_link_libraries(test ${PETSC_LIB})
target_include_directories(test PRIVATE ${PETSC_INC} ${CUDA_INC})
set_target_properties(test PROPERTIES
  CUDA_STANDARD 14
  CUDA_STANDARD_REQUIRED YES
  CUDA_SEPARABLE_COMPILATION ON
  CUDA_ARCHITECTURES 86
)

add_executable(test_func cuda_function_test.cu)
target_link_libraries(test_func ${PETSC_LIB})
target_include_directories(test_func PRIVATE ${PETSC_INC} ${CUDA_INC})
set_target_properties(test_func PROPERTIES
  CUDA_STANDARD 14
  CUDA_STANDARD_REQUIRED YES
  CUDA_SEPARABLE_COMPILATION ON
  CUDA_ARCHITECTURES 86
)