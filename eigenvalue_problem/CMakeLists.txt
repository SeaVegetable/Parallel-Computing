cmake_minimum_required(VERSION 3.10)

project(EigenvalueProblem)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add your source files here
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
            ${SOURCE_DIR}/Partition.cpp)

# Find PETSc
set(PETSC_DIR $ENV{HOME}/lib/petsc-3.18.6-opt)
set(PETSC_ARCH .)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/module")
find_package(PETSc REQUIRED)

set(CMAKE_C_COMPILER  /opt/homebrew/opt/mpich/bin/mpicc)
set(CMAKE_CXX_COMPILER /opt/homebrew/opt/mpich/bin/mpicxx)

include_directories(${PETSC_INC})

add_executable(preprocess ${SOURCES} preprocess.cpp)
add_executable(solve driver.cpp)

target_link_libraries(solve ${PETSC_LIB})