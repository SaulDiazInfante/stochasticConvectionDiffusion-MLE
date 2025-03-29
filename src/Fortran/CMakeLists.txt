cmake_minimum_required(VERSION 3.19)
project(intel_mkl Fortran)
set(ENV{MKLROOT} "/opt/intel/oneapi/mkl/latest")
message(STATUS "ENV{MKLROOT}: $ENV{MKLROOT}")
set(CMAKE_Fortran_COMPILER /opt/intel/oneapi/compiler/2025.1/bin/ifx)
set(CMAKE_Fortran_FLAGS "-qmkl") 
set(MKL_ARCH intel64)
set(MKL_THREADING sequential)
set(MKL_INTERFACE lp64)
find_package(MKL CONFIG REQUIRED PATHS $ENV{MKLROOT})
message(STATUS "Imported oneMKL targets: ${MKL_IMPORTED_TARGETS}")
#set(ENABLE_BLAS95 true)
#set(ENABLE_LAPACK95 true)
# set(MKL_THREADING "tbb_thread")
set(SRC_FILES
    ./src/mod_data_io.f90
    ./src/mod_alloc.f90
    ./src/mod_global_parameters_and_shared_data.f90
    ./src/mod_par_generators.f90
    ./src/mod_sde_coefficients.f90
    ./src/mkl_vsl.f90
    ./src/mod_random_number_generator.f90
    ./src/mod_sde_solver
    ./src/mod_setup_and_logs.f90
    ./src/main.f90
)
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

enable_testing()
include_directories(${MKL_INCLUDE_DIRS})
include_directories(${MKLROOT}/include)
link_directories($ENV{MKLROOT}/lib/intel64)

# Add executable
add_executable(my_fortran_program ${SRC_FILES})


# Add dependencies
add_dependencies(my_fortran_program 
    mod_data_io
    mod_alloc
    mod_global_parameters_and_shared_data
    mod_par_generators
    mod_sde_coefficients    
    mkl_vsl
    mod_random_number_generator
    mod_sde_solver
    mod_setup_and_logs
)
add_library(mod_data_io OBJECT ./src/mod_data_io.f90)
add_library(mod_alloc OBJECT ./src/mod_alloc.f90)
add_library(
    mod_global_parameters_and_shared_data OBJECT ./src/mod_global_parameters_and_shared_data.f90
)
add_library(mod_par_generators OBJECT ./src/mod_par_generators.f90)
add_library(mod_sde_coefficients OBJECT ./src/mod_sde_coefficients.f90)
add_library(mkl_vsl OBJECT ./src/mkl_vsl.f90)
add_library(mod_random_number_generator OBJECT ./src/mod_random_number_generator.f90)
add_library(mod_sde_solver OBJECT ./src/mod_sde_solver.f90)
add_library(mod_setup_and_logs OBJECT ./src/mod_setup_and_logs.f90)
add_dependencies(
    mod_sde_solver 
    mod_random_number_generator 
    mod_par_generators 
    mod_sde_coefficients
    mod_data_io
    mod_alloc
    mod_global_parameters_and_shared_data
    mkl_vsl
    mod_setup_and_logs
    mod_random_number_generator
    mod_sde_coefficients
)
# Include the module directory for compilation
target_include_directories(my_fortran_program PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

# Set compiler options to specify module directory
target_compile_options(my_fortran_program PRIVATE -module ${CMAKE_Fortran_MODULE_DIRECTORY})

target_include_directories(mod_random_number_generator PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})
target_include_directories(mod_par_generators PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})
target_include_directories(mod_sde_coefficients PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})