cmake_minimum_required(VERSION 3.3)
project(parallelClustering)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.c em.c kmeans.c em.h kmeans.h common.h)
add_executable(parallelClustering ${SOURCE_FILES})