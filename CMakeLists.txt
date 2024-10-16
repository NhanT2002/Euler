cmake_minimum_required(VERSION 3.22.1)
project(Euler)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

add_executable(Euler main.cpp
        read_PLOT3D.h
        read_PLOT3D.cpp
        cell.cpp
        cell.h
        SpatialDiscretization.cpp
        SpatialDiscretization.h
        vector_helper.h
        vector_helper.cpp
        TemporalDiscretization.cpp
        TemporalDiscretization.h)
