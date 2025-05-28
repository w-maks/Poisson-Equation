cmake_minimum_required(VERSION 3.29)
project(Poisson_Equation)

set(CMAKE_CXX_STANDARD 17)

include_directories(.)

add_executable(Poisson_Equation
        grid.cpp
        grid.hpp
        poissonSolver.cpp
        poissonSolver.hpp
        main.cpp
        solverSOR.cpp
        solverSOR.hpp
        solverParabolic.cpp
        solverParabolic.hpp
        solverGradient.cpp
        solverGradient.hpp
        solverRandom.cpp
        solverRandom.hpp
)
