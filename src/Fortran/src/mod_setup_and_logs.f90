module mod_setup_and_logs
  use iso_fortran_env, only: int32, real64
  use mod_data_io
  use mod_global_parameters_and_shared_data
  use mod_par_generators
  use mod_sde_coefficients
  use mod_random_number_generator
  use mod_sde_solver
  implicit none

  !! This module wraps the main Fortran functionality to be called from C
  !! It contains the main function and the setup for the SDE solver.
  !! It also includes the logging functionality to write the results to a file.
  !! The module uses the MKL library for random number generation and Gaussian sampling.
  !! The module also includes the setup for the SDE solver, including the generation of
  !! observation times, eigenvalues, and the drift and diffusion matrices.                              
  contains
  !!> @brief Initializes the SDE solver and sets up the parameters.
  !!>
  !!> This subroutine initializes the SDE solver by allocating memory for the
  !!> required matrices and vectors. It also generates the observation times,
  !!> eigenvalues, and the drift and diffusion matrices. The subroutine uses
  !!> the MKL library for random number generation and Gaussian sampling.
  !!>
  !!> @param file_name The name of the file to write the results to.
  !!> @param header The header information to be written to the file.
  !!> @param indexed_times The array to store the observation times.
  !!> @param array The array to store the results of the SDE solver.
  !!> @param vector_drift The vector to store the drift coefficients.
  !!> @param vector_diffusion The vector to store the diffusion coefficients.
  !!> 

  !!> @note This subroutine is called from the main function to set up the SDE solver.  
  subroutine build_sde()
  
    call allocate_dynamic_memory()
    call gen_observation_times()
    call gen_eigen_values()
    call build_matrix_B(eigen_values, B)
    call gen_matrix_diag_B()
    call gen_lambda_matrix()
    call assemble_matrix_A()
    call gen_drift_matrix()
    call gen_diffusion_matrix()

  end subroutine  build_sde

! TODO: code a routine to display constants parameters
  !!> @brief Displays the domain problem arrays.


  subroutine display_domain_problem_arrays()
    call print_vector_with_indices("times", times(1:10), 10)
    call print_vector_with_indices("eigen values", eigen_values(1:10), 10)
    call print_matrix_with_indices("B", B(1:5, 1:5), 5, 5)
    call print_vector_with_indices("diag(B)", B_(1:5), 5)
    call print_matrix_with_indices("Lambda", lambda_matrix(1:5, 1:5), 5, 5)  
    call print_matrix_with_indices("A", A(1:5, 1:5), 5, 5)
    call print_matrix_with_indices("Lambda_diag", B_(1:5), 5, 5)  
    call print_matrix_with_indices("Drift matrix", drift_mat(1:5, 1:5), 5, 5)
    call print_matrix_with_indices(&
      &"Diffusion matrix", &
      &diffusion_mat(1:5, 1:5), &
      &5, &
      &5 &
    &)
  end subroutine display_domain_problem_arrays
end module