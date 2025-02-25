program main
  !! This program illustrate the use of all subroutines in the module
  !! mode_generators
  !! for the model presented in [ref:]
  use iso_fortran_env, only: int32, real64
  use mod_par_generators
  !!use mod_sde_coefficients

  implicit none
  integer(int32), parameter :: Nx = 10
  integer(int32), parameter :: Ny = 10
  integer(int32), parameter :: DIM = Nx * Ny
  integer(int32), parameter :: SEED = 765431
  integer(int32), parameter :: nobs= 1000

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
  real(real64) hs(DIM), startx(DIM), Ls(DIM)
  real(real64) Talpha(DIM,DIM), Tsigma(DIM,DIM), brownian(nobs,DIM), HT
  real(real64) :: times(0:nobs)

! load matrix A entries

open(99,file="MatrixA.dat")
  read(99,*) AM
!    print*,'Dimension=',AM

close(99)


! generate times
  call gen_observation_times(NOBS, DELTA, times)
  print*,"times:", times
  write(*,*)

  call gen_lambdas(DIM, L1, L2, lambda_numbers)
  print*, "lambdas", lambdas
  write(*,*)

  call MB(DIM, lambdas, gamma, B)
  print*, "B", B
  write(*,*)

  call gen_lambda_matrix(DIM, lambdas, lambda_matrix)
  print*,"lambda", lambda_matrix
  write(*,*)

!!  call MA(DIM, hs, A)
!!  print*,"A",A

!!  call MDrift(DIM, theta, beta, lambda_matrix, A, Talpha)
  ! print*,"Talpha", Talpha
  ! TODO: diffusion
end program main

! gfortran main.f90 mod_par_generators.f90 mod_sde_coefficients.f90 -o a.out
