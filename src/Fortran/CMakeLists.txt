cmake_minimum_required(VERSION 3.15)
project(stoCD-MLE Fortran)

# Set the Fortran compiler to Intel ifx
set(CMAKE_Fortran_COMPILER ifx)

# Add source files
add_executable(a.out main.f90 mod_generators.f90)

