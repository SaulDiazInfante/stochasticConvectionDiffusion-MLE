program main
  !! This module wraps the main Fortran functionality to be called from C
  use iso_fortran_env, only: int32, real64
  
  use mod_par_generators
  use mod_sde_coefficients
  use mod_random_number_generator
  implicit none
  integer(int32), parameter :: Nx = 10
  integer(int32), parameter :: Ny = 10
  integer(int32), parameter :: DIM = Nx * Ny
  integer(int32), parameter :: SEED = 765431
  integer(int32), parameter :: nobs= 1000

  real(real64), parameter :: PI = 2.D0 * DASIN(1.D0)
  real(real64),  parameter :: theta = 0.5_real64
  real(real64),  parameter :: beta = 0.5_real64
  real(real64),  parameter :: gamma = 1.0_real64
  real(real64),  parameter :: sigma = 0.2_real64
  real(real64),  parameter :: delta = 0.0001_real64
  real(real64),  parameter :: L1 = 5.0_real64
  real(real64),  parameter :: L2 = 5.0_real64

  real(real64) x, lambda_matrix(DIM,DIM), A(DIM,DIM), path(0:nobs,DIM)
  real(real64) eigen_values(DIM)
  real(real64) lambdas(DIM), B(DIM,DIM), B_(DIM)
  real(real64) hs(DIM), startx(DIM), Ls(DIM), AM(DIM*DIM)
  real(real64) drift_mat(DIM,DIM), diffusion_mat(DIM,DIM)
  real(real64) U(DIM), vector_drift(DIM), vector_diffusion(DIM)
  real(real64) brownian(nobs,DIM), HT
  real(real64) :: times(0:nobs)
  real(kind=8) :: gaussian_sample

! load matrix A entries
  open(99, file="../src/MatrixA.dat")
    read(99,*) AM
  close(99)

  U(:) = 1.0D0
  ! generate times
  call gen_observation_times(NOBS, DELTA, times)
  call print_vector_with_indices("times", times(1:10), 10)

  call gen_lambdas(DIM, Nx, Ny, L1, L2, eigen_values)
  call print_vector_with_indices("eigen values", eigen_values(1:10),10)
  
  call MB(DIM, eigen_values, gamma, B)
  call print_matrix_with_indices("B", B(1:5, 1:5) ,5 ,5)

  call gen_matrix_diag_B(DIM, eigen_values, gamma, B_)
  call print_vector_with_indices("diag(B)", B_(1:5), 5)

  
  call gen_lambda_matrix(DIM, eigen_values, lambda_matrix)
  call print_matrix_with_indices("Lambda", lambda_matrix(1:5, 1:5) ,5 ,5)

  call MA(DIM, Nx, Ny, AM, A)
  call print_matrix_with_indices("A", lambda_matrix(1:5, 1:5) ,5 ,5)

  call gen_drift_matrix(DIM, theta, beta, eigen_values, A, drift_mat)
  call print_matrix_with_indices("Drift matrix", drift_mat(1:5, 1:5) ,5 ,5)

  call  gen_diffusion_matrix(DIM, 1.0_real64, B, diffusion_mat)
  call print_matrix_with_indices("Diffusion matrix", diffusion_mat(1:5, 1:5) ,5 ,5)

  call eval_whole_drift(DIM, beta, theta, eigen_values, A, U, vector_drift)
  call print_vector_with_indices("drift(par, U)", vector_drift(1:5), 5)
  
  call eval_drift_at_u(DIM, beta, theta, drift_mat, U, vector_drift)
  call print_vector_with_indices("drift(U)", vector_drift(1:5), 5)
  
  call eval_diffusion_at_u(DIM, sigma, diffusion_mat, U, vector_diffusion)
  call print_vector_with_indices("diffusion(U)", vector_diffusion(1:5), 5)
  
  mkl_gaussian_sampler(&
      &100, &
      &0.0, &
      &1.0, &
      &123456, &
      &gaussian_sample, &
      &FALSE&
  &)
end program main
