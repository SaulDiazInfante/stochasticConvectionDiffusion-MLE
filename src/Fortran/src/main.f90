program main
  !! This program illustrate the use of all subroutines in the module
  !! mode_generators
  !! for the model presented in [ref:]
  use iso_fortran_env, only: int32, real64
  use mod_par_generators
  use mod_sde_coefficients

  implicit none
  integer(int32), parameter :: Nx = 10
  integer(int32), parameter :: Ny = 10
  integer(int32), parameter :: DIM = Nx * Ny
  integer(int32), parameter :: SEED = 765431
  integer(int32), parameter :: nobs= 1000

  real(real64), parameter :: PI = 2.D0 * DASIN(1.D0)
  real(real64),  parameter :: theta = 0.5
  real(real64),  parameter :: beta = 0.5
  real(real64),  parameter :: gamma = 1.0
  real(real64),  parameter :: sigma = 0.2
  real(real64),  parameter :: delta = 0.0001
  real(real64),  parameter :: L1 = 5.0
  real(real64),  parameter :: L2 = 5.0

  real(real64) x, lambda_matrix(DIM,DIM), A(DIM,DIM), path(0:nobs,DIM)
  real(real64) lambda_numbers(DIM)
  real(real64) lambdas(DIM), B(DIM,DIM)
  real(real64) hs(DIM), startx(DIM), Ls(DIM), AM(DIM*DIM)
  real(real64) drift_mat(DIM,DIM), diffusion_mat(DIM,DIM)
  real(real64) brownian(nobs,DIM), HT
  real(real64) :: times(0:nobs)

! load matrix A entries
!print*,"Fucking PI: ", PI 
  open(99, file="../src/MatrixA.dat")
  read(99,*) AM
  close(99)

! generate times
  call gen_observation_times(NOBS, DELTA, times)
  print*,"times :)"
 
  call gen_lambdas(DIM, Nx, Ny, L1, L2, lambda_numbers)
  print*, "lambdas :)", lambda_numbers(50:70)
 
  call MB(DIM, lambda_numbers, gamma, B)
  print*, "B :)", B(1:5, 1:5)
 
  call gen_lambda_matrix(DIM, lambdas, lambda_matrix)
  print*,"lambda_matrix :)"

  call MA(DIM, Nx, Ny, AM, A)
  print*,"A :)"

  call gen_drift_matrix(DIM, theta, beta, lambda_matrix, A, drift_mat)
  print*,"Drift_matrix :)", drift_mat(1:5, 1:5)
  ! TODO: diffusion
  call  gen_diffusion_matrix(DIM, sigma, B, diffusion_mat)
  print*,"Diffusion_matrix :)", diffusion_mat(1:5, 1:5)
end program main

! ifx -o winner main.f90 -I"${MKLROOT}/include/mkl" -I"${MKLROOT}/include/mkl/intel64/lp64" -qmkl=sequential