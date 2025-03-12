cmake_minimum_required(VERSION 3.10)

project(EigenvalueProblem)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add your source files here
add_executable(preprocess preprocess.cpp)