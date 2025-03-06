cmake_minimum_required(VERSION 3.15)
project(FortranModules LANGUAGES Fortran)

# Define Fortran Compiler
set(CMAKE_Fortran_COMPILER ifx)  # or "ifx" for Intel Fortran Compiler

# Define source files (list modules first)
set(SRC_FILES
    src/mod_par_generators.f90
    src/mod_sde_coefficients.f90
    src/mod_random_number_generator.f90
    src/mkl_vsl.f90
    src/main.f90
)

# Set module output directory
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# Add executable
add_executable(my_fortran_program ${SRC_FILES})
add_dependencies(my_fortran_program 
    mod_par_generators
    mod_sde_coefficients
    mod_random_number_generator
)
add_library(mod_par_generators OBJECT src/mod_par_generators.f90)
add_library(mod_sde_coefficients OBJECT src/mod_sde_coefficients.f90)
add_library(mod_random_number_generator OBJECT 
    src/mod_random_number_generator.f90
#    src/mkl_vsl.f90
)

target_link_libraries(my_fortran_program 
    mod_par_generators
    mod_sde_coefficients
    mod_random_number_generator
)

# Include the module directory for compilation
target_include_directories(my_fortran_program PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})
