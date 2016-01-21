cmake_minimum_required(VERSION 3.3)
project(parallelClustering)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.c em.c kmeans.c)
set(SOURCE_K main.c em.c kmeans.c)

add_executable(parallelClustering ${SOURCE_FILES})
add_executable(kmeans ${SOURCE_K})

target_link_libraries(parallelClustering pthread m)
target_link_libraries(kmeans pthread m)
