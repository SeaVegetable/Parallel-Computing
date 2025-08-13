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

# OpenMP
find_package(OpenMP REQUIRED)

# MPI
find_package(MPI REQUIRED)

# Add your source files here
set(CUDA_INC "/usr/local/cuda/include")
set(EIGEN3_INCLUDE_DIR "$ENV{HOME}/lib/eigen-3.4.0")
set(INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

message(STATUS "EIGEN3_INCLUDE_DIR: ${EIGEN3_INCLUDE_DIR}")
message(STATUS "INCLUDE_DIR: ${INCLUDE_DIR}")
message(STATUS "SOURCE_DIR: ${SOURCE_DIR}")

# Include directories
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${INCLUDE_DIR})

# Add source files
set(SOURCES ${SOURCE_DIR}/BSplineBasis.cpp
            ${SOURCE_DIR}/ControlPointGenerator.cpp
            ${SOURCE_DIR}/IENGenerator.cpp
            ${SOURCE_DIR}/IDGenerator.cpp
            ${SOURCE_DIR}/NURBSExtractionGenerator.cpp
            ${SOURCE_DIR}/Partition.cpp
            ${SOURCE_DIR}/FileManager.cpp
            ${SOURCE_DIR}/BernsteinBasis.cpp
            ${SOURCE_DIR}/QuadraturePoint.cpp
            ${SOURCE_DIR}/LocalAssembly.cpp
            ${SOURCE_DIR}/GlobalAssemblyDR.cpp
            ${SOURCE_DIR}/AbscissaeGenerator.cpp
            ${SOURCE_DIR}/ElementFEM.cpp
            ${SOURCE_DIR}/Elem2COOGenerator.cpp)

set(CUDA_SOURCES ${SOURCE_DIR}/assembly.cu
                 ${SOURCE_DIR}/GlobalAssembly.cu)

# Find PETSc
set(PETSC_DIR $ENV{HOME}/lib/petsc-3.23.4-opt)
set(PETSC_ARCH .)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/module")
find_package(PETSc REQUIRED)
message(STATUS "PETSC_INC: ${PETSC_INC}")
message(STATUS "PETSC_LIB: ${PETSC_LIB}")

set(CMAKE_C_COMPILER $ENV{HOME}/lib/openmpi-5.0.8/bin/mpicc)
set(CMAKE_CXX_COMPILER $ENV{HOME}/lib/openmpi-5.0.8/bin/mpicxx)

include_directories(${PETSC_INC})

add_library(eigenvalue_problem_lib ${SOURCES})
add_library(eigenvalue_problem_cuda_lib ${CUDA_SOURCES})
target_include_directories(eigenvalue_problem_cuda_lib PRIVATE ${MPI_C_INCLUDE_PATH})
set_target_properties(eigenvalue_problem_cuda_lib PROPERTIES
  CUDA_STANDARD 11
  CUDA_STANDARD_REQUIRED YES
  CUDA_SEPARABLE_COMPILATION ON
  CUDA_ARCHITECTURES 75
)

add_executable(preprocess_fem preprocess_fem.cpp)
target_link_libraries(preprocess_fem eigenvalue_problem_lib)

add_executable(dry_run dry_run.cpp)
target_link_libraries(dry_run eigenvalue_problem_lib ${PETSC_LIB})

add_executable(coo2aij coo2aij.cu)
target_link_libraries(coo2aij eigenvalue_problem_lib ${PETSC_LIB})
target_include_directories(coo2aij PRIVATE ${PETSC_INC} ${CUDA_INC} ${MPI_C_INCLUDE_PATH})

add_executable(memorytest memorytest.cu)
target_link_libraries(memorytest eigenvalue_problem_lib ${PETSC_LIB})
target_include_directories(memorytest PRIVATE ${PETSC_INC} ${CUDA_INC} ${MPI_C_INCLUDE_PATH})

add_executable(assemblytest assemblytest.cu)
target_link_libraries(assemblytest eigenvalue_problem_lib eigenvalue_problem_cuda_lib ${PETSC_LIB})
target_include_directories(assemblytest PRIVATE ${PETSC_INC} ${CUDA_INC} ${MPI_C_INCLUDE_PATH})

if(OPENMP_CXX_FOUND)
  set_target_properties( eigenvalue_problem_lib PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DUSE_OPENMP" )
  set_target_properties( preprocess_fem PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DUSE_OPENMP" )
  target_include_directories( eigenvalue_problem_lib PRIVATE ${OpenMP_CXX_INCLUDE_DIR} )
  target_include_directories( preprocess_fem PRIVATE ${OpenMP_CXX_INCLUDE_DIR} )
  target_link_libraries( preprocess_fem ${OpenMP_CXX_LIBRARIES} )
endif()

set_target_properties(coo2aij PROPERTIES
  CUDA_STANDARD 11
  CUDA_STANDARD_REQUIRED YES
  CUDA_SEPARABLE_COMPILATION ON
  CUDA_ARCHITECTURES 75
)

set_target_properties(memorytest PROPERTIES
  CUDA_STANDARD 11
  CUDA_STANDARD_REQUIRED YES
  CUDA_SEPARABLE_COMPILATION ON
  CUDA_ARCHITECTURES 75
)

set_target_properties(assemblytest PROPERTIES
  CUDA_STANDARD 11
  CUDA_STANDARD_REQUIRED YES
  CUDA_SEPARABLE_COMPILATION ON
  CUDA_ARCHITECTURES 75
)