cmake_minimum_required(VERSION 3.10)

project(SphereMeshBlendShape VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(GLM_DIR "${CMAKE_SOURCE_DIR}/external/glm")

add_library(SphereMeshBlendShape STATIC
        include/sphere_mesh.h
        src/sphere_mesh.cpp
        include/bumper_graph.h
        src/bumper_graph.cpp
        include/bumper_grid.h
        src/bumper_grid.cpp
)

target_include_directories(SphereMeshBlendShape PUBLIC include ${GLM_DIR})