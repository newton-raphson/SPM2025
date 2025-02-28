cmake_minimum_required(VERSION 3.5.0)
project(Ray_tracing)

set(CMAKE_CXX_STANDARD 11)

add_library(raytracer
        include/raytracer/stl_reader.h
        src/stl_reader.cpp
        include/raytracer/ray_tracer.h
        src/ray_tracer.cpp
        include/raytracer/ray_tracer_2d.h
        src/ray_tracer_2d.cpp
        include/raytracer/math_types.h
        include/testpoints.h
        include/raytracer/gmsh_2d_reader.h
        src/gmsh_2d_reader.cpp
        )
target_include_directories(raytracer PUBLIC include)

option(RAYTRACER_BUILD_EXAMPLE "Build the Raytracer example program" OFF)
option(BENCHMARK "Enabling Raytracing Benchmark" OFF)
option(BENCHMARKLOG "Enabling Benchmark Log" OFF)
if(BENCHMARK)
    add_definitions(-DBENCHMARK)
    message("Enabling Raytracing Benchmark")
    set(BENCHMARK,ON)
    if(BENCHMARKLOG)
        add_definitions(-DBENCHMARKLOG)
        message("Enabling Benchmark Log")
        set(BENCHMARKLOG,ON)
    endif()
endif()



set(CMAKE_BUILD_TYPE "Release")
if (RAYTRACER_BUILD_EXAMPLE)
    add_executable(raytracing-example main.cpp)
    target_link_libraries(raytracing-example raytracer)
endif()
