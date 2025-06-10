cmake_minimum_required(VERSION 3.20)
project(intel_mkl LANGUAGES Fortran)

# Set Intel oneAPI compiler
set(CMAKE_Fortran_COMPILER /opt/intel/oneapi/compiler/2025.1/bin/ifx)
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# Set MKLROOT
set(ENV{MKLROOT} "/opt/intel/oneapi/mkl/latest")
message(STATUS "MKLROOT = $ENV{MKLROOT}")

# Enable Debug build type
set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)

# Compiler flags for debug mode
# -g: debugging info, -O0: no optimization, -fpp: preprocess, -qmkl: link MKL
set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -fpp -qmkl")

# Source directory
set(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(TEST_DIR "${CMAKE_CURRENT_SOURCE_DIR}/test")

# Define all module libraries
add_library(mod_alloc STATIC ${SRC_DIR}/mod_alloc.f90)
add_library(mod_global_parameters_and_shared_data STATIC ${SRC_DIR}/mod_global_parameters_and_shared_data.f90)
add_library(mod_par_generators STATIC ${SRC_DIR}/mod_par_generators.f90)
add_library(mod_sde_coefficients STATIC ${SRC_DIR}/mod_sde_coefficients.f90)
add_library(mod_data_io STATIC ${SRC_DIR}/mod_data_io.f90)
add_library(mkl_vsl STATIC ${SRC_DIR}/mkl_vsl.f90)
add_library(mod_random_number_generator STATIC ${SRC_DIR}/mod_random_number_generator.f90)
add_library(mod_sde_solver STATIC ${SRC_DIR}/mod_sde_solver.f90)
add_library(mod_setup_and_logs STATIC ${SRC_DIR}/mod_setup_and_logs.f90)

# Add module dependencies
target_include_directories(mod_alloc PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

target_link_libraries(mod_global_parameters_and_shared_data PRIVATE mod_alloc)
target_include_directories(mod_global_parameters_and_shared_data PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

target_link_libraries(mod_par_generators PRIVATE mod_global_parameters_and_shared_data)
target_include_directories(mod_par_generators PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

target_link_libraries(mod_sde_coefficients PRIVATE mod_alloc mod_global_parameters_and_shared_data)
target_include_directories(mod_sde_coefficients PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

target_include_directories(mod_data_io PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

target_include_directories(mkl_vsl PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY} $ENV{MKLROOT}/include)

target_link_libraries(mod_random_number_generator PRIVATE mkl_vsl)
target_include_directories(mod_random_number_generator PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

target_link_libraries(mod_sde_solver PRIVATE 
    mod_random_number_generator 
    mod_par_generators 
    mod_sde_coefficients
    mod_global_parameters_and_shared_data
)
target_include_directories(mod_sde_solver PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

target_link_libraries(mod_setup_and_logs PRIVATE 
    mod_data_io
    mod_global_parameters_and_shared_data
    mod_par_generators
    mod_sde_coefficients
    mod_random_number_generator
    mod_sde_solver
)
target_include_directories(mod_setup_and_logs PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

# Ensure correct dependency order
set(MY_LIBS
        mod_alloc
        mod_global_parameters_and_shared_data
        mod_par_generators
        mod_sde_coefficients
        mod_data_io
        mkl_vsl
        mod_random_number_generator
        mod_sde_solver
        mod_setup_and_logs
)

# Add main program
add_executable(simulate_sde_path ${SRC_DIR}/main.f90)
target_link_libraries(simulate_sde_path PRIVATE ${MY_LIBS})
target_include_directories(simulate_sde_path PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY} $ENV{MKLROOT}/include)

# Enable testing infrastructure
enable_testing()

# Add test program for mod_data_io
add_executable(test_mod_data_io ${TEST_DIR}/test_mod_data_io.f90)
target_link_libraries(test_mod_data_io PRIVATE mod_data_io)
target_include_directories(test_mod_data_io PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})
add_test(NAME TestModDataIO COMMAND test_mod_data_io)

# Add test program for mod_alloc
add_executable(test_mod_alloc ${TEST_DIR}/test_mod_alloc.f90)
target_link_libraries(test_mod_alloc PRIVATE mod_alloc)
target_include_directories(test_mod_alloc PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})
add_test(NAME TestModAlloc COMMAND test_mod_alloc)

# Add test program for mod_sde_solver
add_executable(test_mod_sde_solver ${TEST_DIR}/test_mod_sde_solver.f90)
target_link_libraries(test_mod_sde_solver PRIVATE
        mod_alloc
        mod_global_parameters_and_shared_data
        mod_par_generators
        mod_sde_coefficients
        mod_data_io
        mkl_vsl
        mod_random_number_generator
        mod_sde_solver
        mod_setup_and_logs
)
target_include_directories(test_mod_sde_solver PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})
add_test(NAME TestModSDESolver COMMAND test_mod_sde_solver)

# Add test program for mod_random_number_generator
add_executable(test_mod_random_number_generator ${TEST_DIR}/test_mod_random_number_generator.f90)
target_link_libraries(test_mod_random_number_generator PRIVATE mod_random_number_generator mkl_vsl)
target_include_directories(test_mod_random_number_generator PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})
add_test(NAME TestModRandomNumberGenerator COMMAND test_mod_random_number_generator)

# Add test program for mod_sde_coefficients
add_executable(test_mod_sde_coefficients ${TEST_DIR}/test_mod_sde_coefficients.f90)
target_link_libraries(test_mod_sde_coefficients PRIVATE
        mod_alloc
        mod_global_parameters_and_shared_data
        mod_par_generators
        mod_sde_coefficients
        mod_data_io
        mkl_vsl
        mod_random_number_generator
        mod_sde_solver
        mod_setup_and_logs
)
target_include_directories(test_mod_sde_coefficients PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})
add_test(NAME TestModSDECoefficients COMMAND test_mod_sde_coefficients)

