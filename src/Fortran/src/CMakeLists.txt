cmake_minimum_required(VERSION 3.19)
project(intel_mkl Fortran)
set(ENV{MKLROOT} "/opt/intel/oneapi/mkl/latest")
message(STATUS "ENV{MKLROOT}: $ENV{MKLROOT}")
#set(MKL_LINK "static")
# variables set before find_package(MKL)
# https://www.intel.com/content/www/us/en/docs/onemkl/developer-guide-linux/2023-2/cmake-config-for-onemkl.html
set(CMAKE_Fortran_COMPILER ifx)
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
    mod_sde_coefficients.f90
    mod_par_generators.f90
    mkl_vsl.f90
    mod_random_number_generator.f90
    mod_sde_solver
    main.f90
)
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}")
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

enable_testing()
include_directories(${MKL_INCLUDE_DIRS})
include_directories(${MKLROOT}/include)

# Add executable
add_executable(my_fortran_program ${SRC_FILES})

# Link MKL libraries
target_link_libraries(my_fortran_program ${MKL_LIBRARIES})

# Add dependencies
add_dependencies(my_fortran_program 
    mod_sde_coefficients    
    mod_par_generators
    mkl_vsl
    mod_random_number_generator
    mod_sde_solver
)
add_library(mod_sde_coefficients OBJECT mod_sde_coefficients.f90)
add_library(mod_par_generators OBJECT mod_par_generators.f90)
add_library(mkl_vsl OBJECT mkl_vsl.f90)
add_library(mod_random_number_generator OBJECT mod_random_number_generator.f90)
add_library(mod_sde_solver OBJECT mod_sde_solver.f90)
# Include the module directory for compilation
target_include_directories(my_fortran_program PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

# Set compiler options to specify module directory
target_compile_options(my_fortran_program PRIVATE -module ${CMAKE_Fortran_MODULE_DIRECTORY})