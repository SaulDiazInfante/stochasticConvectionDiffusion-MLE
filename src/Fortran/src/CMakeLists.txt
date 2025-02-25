cmake_minimum_required(VERSION 3.15)
project(stoCD-MLE Fortran)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# Set the Fortran compiler to Intel ifx
set(CMAKE_Fortran_COMPILER ifx)
set(CMAKE_Fortran_FLAGS_DEBUG "-g")

# Add source files
add_executable(a.out main.f90 mod_par_generators.f90)


