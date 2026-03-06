cmake_minimum_required(VERSION 3.30)

set(PROJECT_NAME Spin2DLattice)
project(
    ${PROJECT_NAME}   
)

add_executable(${PROJECT_NAME}
    include/main.cpp
    include/core/data.cpp 
    include/core/simulation.cpp 
    include/lattice/atoms.cpp 
    include/lattice/geometry.cpp 
    include/physics/calculator.cpp 
)