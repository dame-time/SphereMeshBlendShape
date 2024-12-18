cmake_minimum_required(VERSION 3.10)

# Project name
project(SphereMeshBlendShape VERSION 1.0.0)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the GLM submodule
set(GLM_DIR "${CMAKE_SOURCE_DIR}/external/glm")

# Add the library
add_library(SphereMeshBlendShape STATIC
        src/sphere_mesh.cpp
        include/sphere_mesh.h
        include/bumper_graph.h
        src/bumer_grid.cpp
        src/bumper_graph.cpp
)

# Include the header files
target_include_directories(SphereMeshBlendShape PUBLIC include ${GLM_DIR})

# Add an executable for testing the library (optional)
add_executable(SphereMeshBlendShapeTest main.cpp)

# Link the library to the test executable
target_link_libraries(SphereMeshBlendShapeTest PRIVATE SphereMeshBlendShape)
target_compile_options(SphereMeshBlendShapeTest PRIVATE
        -O3        # Maximum optimization
        -march=native # Optimize for the host machine
        -flto      # Link time optimization
        -ffast-math # Aggressive math optimizations
        -funroll-loops # Unroll loops for performance
        -DNDEBUG   # Disable asserts
)