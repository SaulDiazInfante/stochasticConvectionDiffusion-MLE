!! ifx -qmkl mod_alloc.f90 mod_data_io.f90 mod_global_parameters_and_shared_data.f90 mod_par_generators.f90 mod_sde_coefficients.f90 mkl_vsl.f90 mod_random_number_generator.f90 mod_sde_solver.f90 main.f90
program main
  !! This module wraps the main Fortran functionality to be called from C
  use iso_fortran_env, only: int32, real64
  use mod_data_io
  use mod_global_parameters_and_shared_data
  use mod_par_generators
  use mod_sde_coefficients
  use mod_random_number_generator
  use mod_setup_and_logs
  use mod_sde_solver
  
  implicit none
  
  character(len=20), dimension(2) :: header 
  character(len=50) :: file_name 
  real(real64), dimension(DIM, 2) :: indexed_times
  real(real64), allocatable :: array(:,:), vector_drift(:), vector_diffusion(:)
  real(real64), allocatable :: vectorial_winner_delta(:), initial_vector_winner(:)
  print *, "DIM: ", DIM
  print *, "Nx: ", Nx
  print *, "Ny: ", Ny
  print *, "SEED: ", SEED
  print *, "nobs: ", nobs
  print *, "PI: ", PI   
  print *, "theta: ", theta
  print *, "beta: ", beta
  print *, "gamma: ", gamma
  print *, "sigma: ", sigma
  print *, "delta: ", delta
  print *, "L1: ", L1
  print *, "L2: ", L2 

  call build_sde()
  call display_domain_problem_arrays()
  
  ! Allocate arrays for calculations
  allocate(vector_drift(DIM))
  allocate(vector_diffusion(DIM))
  
  u(:) = 1.0_real64
  call eval_whole_drift(u, vector_drift)
  call print_vector_with_indices("drift(par, u)", vector_drift(1:5), 5)
  
  call eval_drift_at_u(u, vector_drift)
  call print_vector_with_indices("drift(u)", vector_drift(1:5), 5)
  
  call eval_diffusion_at_u(u, vector_diffusion)
  call print_vector_with_indices("diffusion(u)", vector_diffusion(1:5), 5)
  call eval_diagonal_diffusion_at_u(u, vector_diffusion)
  call print_vector_with_indices(&
    &"diffusion(U) from diag(B)", &
    & vector_diffusion(1:5), 5 &
  &)
  
  ! mean_a = 0.0
  ! std_a = 1.0
  ! call mkl_gaussian_sampler(10000, mean_a, std_a, gaussian_sample, SEED)
  ! call print_vector_with_indices(&
  !   &"Gaussian(mu, std)",&
  !   &gaussian_sample(9000:9010), &
  !   &10)

  call scalar_winner_increment(0.0_real64, winner_delta)
  ! print *, "Winner delta: ", winner_delta 
  ! call vectorial_winner_increment(&
  !   & 0.1_real64, &
  !   & DIM, &
  !   & 100, & 
  !   &initial_vector_winner, &
  !   & vectorial_winner_delta &
  ! &)
  ! call print_vector_with_indices(&
  !   &"Winner delta", &
  !   & vectorial_winner_delta(90:95), &
  !   & 5 &
  ! &)

  ! call milstein_step(&
  !   &DIM, &
  !   &delta, &
  !   &beta, &
  !   &theta, &
  !   &driftmat, &
  !   &sigma, &
  !   &vector_diffusion, &
  !   &u, &
  !   &vectorial_winner_delta, &
  !   &u &  
  ! &)
  ! call print_vector_with_indices("U_milstein", u(1:DIM), DIM)
  ! Deallocate the arrays we allocated
  deallocate(vector_drift)
  deallocate(vector_diffusion)
end program main
