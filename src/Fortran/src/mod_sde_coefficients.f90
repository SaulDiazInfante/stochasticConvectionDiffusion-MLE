!> @ingroup modules
!> @author F. Baltazar-Larios, S. Diaz-Infante
!> @brief This module computes the coefficients of a SDE regarding
!> the PDE model.
module mod_sde_coefficients
  use iso_fortran_env, only: int32, real64
  implicit none

  ! Include the MKL module
  include 'mkl_blas.fi'
contains
  subroutine gen_drift_matrix(DIM, theta, beta, lambda_mat, A, drift_mat)
    implicit none
    integer(int32), intent(in) ::  DIM
    !! dimension $(N_x * N_y)$
    real(real64), intent(in) :: theta
    !! meaning
    real(real64), intent(in) :: beta
    !! meaning
    real(real64), intent(in) :: A(DIM, DIM)
    !! meaning
    real(real64), intent(in) :: lambda_mat(DIM, DIM)
    !!
    real(real64), intent(out) :: drift_mat(DIM,DIM)
    
    drift_mat = beta * lambda_mat(:,:) + theta * A(:, :)
    drift_mat = -1.0_real64 * drift_mat
    return
  end subroutine gen_drift_matrix

  subroutine gen_diffusion_matrix(DIM, sigma, B, diffusion_mat)
    implicit none
    integer(int32), intent(in) :: DIM
    !!>  dimension
    real(real64), intent(in) :: sigma
    !! Noise intensity
    real(real64), intent(in) :: B(DIM, DIM)

    real(real64), intent(out) :: diffusion_mat(DIM, DIM)

    diffusion_mat(:, :) = sigma * B(:, :)
    return
  end subroutine gen_diffusion_matrix

!> @brief compute the drift coefficient of SDE equation
!!  \f $dU = -(\beta \Lambda  + \theta A) U dt + \sigma B U dW(t)\f$. 
!!
!! @param[in]   DIM           int32 matrix dimension
!! @param[in]   beta          real64
!! @param[in]   theta         real64
!! @param[in]   drift_matrix  real64(DIM, DIM) current values in the array
!! @param[in]   U
!! @param[out]  vector_drift  \f $-(\beta \Lambda  + \theta A) U dt
  subroutine eval_drift(DIM, beta, theta, Lambda_diagonal, A_matrix, U, vector_drift)
    implicit none
    integer(int32), intent(in) :: DIM
    real(real64), intent(in) :: beta
    real(real64), intent(in) :: theta
    real(real64), intent(inout) :: Lambda_diagonal(DIM)
    real(real64), intent(inout) :: A_matrix(DIM, DIM)
    real(real64), intent(in) :: U(DIM) 
    real(real64), intent(out) :: vector_drift(DIM)
    
    character(1) :: trans
    integer(int32) :: i

    trans = 'N'  ! No transpose
    
    Lambda_diagonal(:) = beta * Lambda_diagonal(:)
    do i = 1, DIM
      A_matrix(i, i) = A_matrix(i, i) + Lambda_diagonal(i)
    end do
    
    call dgemv(trans, DIM, DIM, -1.0_real64 * theta , A_matrix, DIM, U, 1, 0.0_real64, vector_drift, 1)
  
  end subroutine eval_drift
  
!> @brief compute the diffusion coefficient of SDE equation
!!  \f $dU = -(\beta \Lambda  + \theta A) U dt + \sigma B U dW(t)\f$. 
!!
!! @param[in]   DIM           int32 matrix dimension
!! @param[in]   sigma real64
!! @param[in]   diffusion_matrix  real64(DIM, DIM) current values in the array
!! @param[in]   U                 real64(DIM)
!! @param[out]  vector_diffusion  real64(DIM  \f $\sigma B U \f$. 

  subroutine eval_diffusion(DIM, sigma, diffusion_matrix, U, vector_diffusion)
    implicit none
    integer(int32), intent(in) :: DIM
    real(real64), intent(in) :: sigma
    real(real64), intent(in) :: diffusion_matrix(DIM, DIM)
    real(real64), intent(in) :: U(DIM) 
    real(real64), intent(out) :: vector_diffusion(DIM)

    real(real64), allocatable :: A(:,:), x(:), y(:)
    real(real64) :: alpha_, beta_
    character(1) :: trans
    integer(int32) :: info
    allocate(A(DIM,DIM), x(DIM), y(DIM), stat=info)
    if (info /= 0) then
        write(*,*) "Error: Memory allocation failed!"
        stop
    end if
    A(:,:) = diffusion_matrix(:,:)
    x(:) = U(:)
    y(:) = 0.0_real64    
    
    ! Set parameters for DGEMV
    ! y = alpha*A*x + beta*y
    
    alpha_ = sigma
    beta_ = 0.0_real64
    trans = 'N'  ! No transpose
    
    call dgemv(trans, DIM, DIM, alpha_, A, DIM, x, 1, beta_, y, 1)
    vector_diffusion(:) = y(:)
    deallocate(A, x, y)
  end subroutine eval_diffusion
end module mod_sde_coefficients

