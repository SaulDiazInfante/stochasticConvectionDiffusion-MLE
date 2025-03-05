!> To compile with intel fortran and mkl use
!>  
!> ifx -qmkl mod_sde_coefficients.f90 mod_par_generators.f90 main_.f90


program main
    !! This module wraps the main Fortran functionality to be called from C
    use iso_fortran_env, only: int32, real64
    use iso_c_binding, only: c_int
    use mod_par_generators
    use mod_sde_coefficients
  
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
    real(real64) lambda_numbers(DIM)
    real(real64) lambdas(DIM), B(DIM,DIM), B_(DIM)
    real(real64) hs(DIM), startx(DIM), Ls(DIM), AM(DIM*DIM)
    real(real64) drift_mat(DIM,DIM), diffusion_mat(DIM,DIM)
    real(real64) U(DIM), vector_drift(DIM), vector_diffusion(DIM)
    real(real64) brownian(nobs,DIM), HT
    real(real64) :: times(0:nobs)
  
  ! load matrix A entries
    open(99, file="../src/MatrixA.dat")
    read(99,*) AM
    close(99)
  
    U(:) = 1.0D0
    ! generate times
    call gen_observation_times(NOBS, DELTA, times)
    print*,"(+++) times:)"
    call print_vector_with_indices(times(1:10), 10)

    call gen_lambdas(DIM, Nx, Ny, L1, L2, lambda_numbers)
    print*,"(+++) lambda eigenvalues :)"
    call print_vector_with_indices(lambda_numbers(1:10),10)
   
    call MB(DIM, lambda_numbers, gamma, B)
    print*, "(+++) B :)" 
    call print_matrix_with_indices(B(1:5, 1:5) ,5 ,5)

    call gen_matrix_diag_B(DIM, lambda_numbers, gamma, B_)
    print*, "(++++) diag(B) :)"
    call print_vector_with_indices(B_(1:5), 5)

   
    call gen_lambda_matrix(DIM, lambdas, lambda_matrix)
    print*,"lambda_matrix :)"
  
    call MA(DIM, Nx, Ny, AM, A)
    print*,"A :)"
  
    call gen_drift_matrix(DIM, theta, beta, lambda_numbers, A, drift_mat)
    print*,"Drift_matrix :)", drift_mat(1:5, 1:5)
  
    call  gen_diffusion_matrix(DIM, 1.0_real64, B, diffusion_mat)
    print*,"Diffusion_matrix :)", diffusion_mat(1:5, 1:5)
  
    call eval_whole_drift(DIM, beta, theta, lambda_numbers, A, U, vector_drift)
    print*, "whole drift :)", vector_drift(1:5)
    
    call eval_drift_at_u(DIM, beta, theta, drift_mat, U, vector_drift)
    print*, "drift :)", vector_drift(1:5)
    
    call eval_diffusion_at_u(DIM, sigma, diffusion_mat, U, vector_diffusion)
    print*, "diffusion :)", vector_diffusion(1:5)
  end program main
  