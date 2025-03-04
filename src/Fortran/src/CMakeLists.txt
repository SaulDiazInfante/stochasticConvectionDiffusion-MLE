cmake_minimum_required(VERSION 3.14)
project(FortranMKLExample LANGUAGES Fortran C)
# Enable verbose output for debugging
set(CMAKE_VERBOSE_MAKEFILE ON)

# Set CMAKE_SHARED_LIBRARY_LINK_C_FLAGS to empty to avoid linker warnings
set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")

# Find Intel OneAPI MKL
find_package(MKL REQUIRED)

# Set Fortran compiler (Intel ifx)
set(CMAKE_Fortran_COMPILER ifx)  

# Compiler flags for optimization
set(CMAKE_Fortran_FLAGS "-O3 -xHost -qopenmp -m64 -fpp -fPIC")

# Set MKL include and link flags
include_directories(${MKL_INCLUDE_DIRS})

# C compiler flags
set(CMAKE_C_FLAGS "-O3 -m64 -fopenmp -fPIC")
link_directories(${MKL_LIBRARY_DIRS} /opt/intel/oneapi/compiler/2025.0/lib)

# Source files
set(SOURCES main.f90 mod_par_generators.f90 mod_sde_coefficients.f90 c_main.c)

# Create the executable
add_executable(main ${SOURCES})

# Link with MKL and Intel Fortran runtime libraries
target_link_libraries(main MKL::MKL ifcore ifport /opt/intel/oneapi/compiler/2025.0/lib/libsvml.so /opt/intel/oneapi/compiler/2025.0/lib/libintlc.so.5)

# Ensure correct runtime paths
set_target_properties(main PROPERTIES 
    LINK_FLAGS "-Wl,-rpath,${MKL_LIBRARY_DIRS}" 
    LINKER_LANGUAGE C)

# Mark Fortran modules for proper dependency tracking
set_source_files_properties(
    mod_par_generators.f90 
    mod_sde_coefficients.f90
    PROPERTIES COMPILE_FLAGS "-module ${CMAKE_CURRENT_BINARY_DIR}")

# Add include directory for Fortran modules
include_directories(${CMAKE_CURRENT_BINARY_DIR})
