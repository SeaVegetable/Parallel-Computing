cmake_minimum_required(VERSION 3.10)

project(matmul)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(matmul main.cpp)
#target_link_libraries(main pthread)