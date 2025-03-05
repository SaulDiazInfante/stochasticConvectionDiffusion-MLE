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

!> @brief Generates the drift matrix for a stochastic differential equation.
!>
!> The drift matrix is computed as:
!> @f[
!> \text{drift_mat} = - \left( \beta \Lambda + \theta A \right)
!> @f]
!> where:
!> - @f$ \Lambda @f$ is a diagonal matrix.
!> - @f$ A @f$ is a general interaction matrix.
!> - @f$ \beta, \theta @f$ are scalar parameters.
!>
!> @param[in]  DIM             Dimension of the square matrices.
!> @param[in]  theta           Parameter controlling the influence of matrix A.
!> @param[in]  beta            Parameter controlling the influence of Lambda_diagonal.
!> @param[in]  Lambda_diagonal A vector containing the diagonal elements of the Lambda matrix.
!> @param[in]  A               The interaction matrix of size (DIM, DIM).
!> @param[out] drift_mat       The resulting drift matrix of size (DIM, DIM).
subroutine gen_drift_matrix(DIM, theta, beta, Lambda_diagonal, A, drift_mat)
    implicit none
    integer(int32), intent(in) ::  DIM
    real(real64), intent(in) :: theta
    real(real64), intent(in) :: beta
    real(real64), intent(in) :: A(DIM, DIM)
    real(real64), intent(in) :: Lambda_diagonal(DIM)
    real(real64), intent(out) :: drift_mat(DIM,DIM)
    
    integer(int32) :: i
    
    drift_mat(:, :) = theta * A(:, :)
    
    do i=1, DIM
      drift_mat(i, i) = drift_mat(i, i) + Lambda_diagonal(i)
    end do
    drift_mat = -1.0_real64 * drift_mat
    return
  end subroutine gen_drift_matrix

!> @brief Generates the diffusion matrix for a stochastic differential equation.
!> @param[in]  DIM           Dimension of the square matrices.
!> @param[in]  sigma         Noise intensity parameter.
!> @param[in]  B             A matrix of size (DIM, DIM) influencing the 
!>                           diffusion term.
!> @param[out] diffusion_mat The resulting diffusion matrix of size (DIM, DIM).

  subroutine gen_diffusion_matrix(DIM, sigma, B, diffusion_mat)
    implicit none
    integer(int32), intent(in) :: DIM
    real(real64), intent(in) :: sigma
    real(real64), intent(in) :: B(DIM, DIM)
    real(real64), intent(out) :: diffusion_mat(DIM, DIM)
    diffusion_mat(:, :) = sigma * B(:, :)
    return
  end subroutine gen_diffusion_matrix

!> @brief 
!> Given Matrix A, vector Lambda and parameters beta, theta, this   
!> subroutine computes the evaluation of the drift coefficient in the SDE
!> at vector U and conforming to the SDE
!>  @f$
!>    dU = -(\beta \Lambda  + \theta A) \ U dt + \sigma B\ U dW(t) 
!>  @f$
!>
!> @param[in]   DIM           int32 matrix dimension
!> @param[in]   beta          real64
!> @param[in]   theta         real64
!> @param[in]   Lambda        real64(DIM) Diagonal of matrix Lambda 
!> @param[in]   A             real64(DIM, DIM) Matrix from spectral 
!>
!> decomposition
!> @param[in]   U             real64(DIM)
!> @param[out]  vector_drift  @f$ -(\beta \Lambda  + \theta A) U @f$
  subroutine eval_whole_drift(&
    & DIM, & 
    & beta, &
    & theta, &
    & Lambda_diagonal, &
    & A_matrix, &
    & U, &
    & vector_drift &
  & )
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
  
  end subroutine eval_whole_drift
!> @brief 
!> Given Matrix A, vector Lambda and parameters beta, theta, this   
!> subroutine computes the evaluation of the drift coefficient in the SDE
!> at vector U and conforming to the following SDE
!!  @f$ dU = -(\beta \Lambda  + \theta A) U dt + \sigma B U dW(t) @f$. 
!!
!! @param[in]   DIM           int32 matrix dimension
!! @param[in]   beta          real64
!! @param[in]   theta         real64
!! @param[in]   Lambda        real64(DIM) Diagonal of matrix Lambda 
!! @param[in]   A             real64(DIM, DIM) Matrix from spectral 
!! decomposition
!! @param[in]   U             real64(DIM)
!! @param[out]  vector_drift  @f$-(\beta \Lambda  + \theta A) U dt @f$

subroutine eval_drift_at_u(&
    & DIM, & 
    & beta, &
    & theta, &
    & drift_matrix, &
    & U, &
    & vector_drift &
  & )
    implicit none
    integer(int32), intent(in) :: DIM
    real(real64), intent(in) :: beta
    real(real64), intent(in) :: theta
    real(real64), intent(out) :: drift_matrix(DIM, DIM)
    real(real64), intent(in) :: U(DIM) 
    real(real64), intent(out) :: vector_drift(DIM)
    
    character(1) :: trans
    
    trans = 'N'  ! No transpose
            
    call dgemv(&
      & trans, DIM, DIM, &
      & 1.0_real64 , drift_matrix, DIM,&
      & U, 1, 0.0_real64, vector_drift, 1)
  end subroutine eval_drift_at_u



!> @brief compute the diffusion coefficient of SDE equation
!!  @f$ dU = -(\beta \Lambda  + \theta A) U dt + \sigma B U dW(t) @f$. 
!!
!! @param[in]   DIM           int32 matrix dimension
!! @param[in]   sigma real64
!! @param[in]   diffusion_matrix  real64(DIM, DIM) current values in the array
!! @param[in]   U                 real64(DIM)
!! @param[out]  vector_diffusion  real64(DIM  \f $\sigma B U \f$. 

  subroutine eval_diffusion_at_u(DIM, sigma, diffusion_matrix, U, vector_diffusion)
    implicit none
    integer(int32), intent(in) :: DIM
    real(real64), intent(in) :: sigma
    real(real64), intent(in) :: diffusion_matrix(DIM, DIM)
    real(real64), intent(in) :: U(DIM) 
    real(real64), intent(out) :: vector_diffusion(DIM)

    character(1) :: trans
    integer(int32) :: info 
    ! Set parameters for DGEMV
    ! y = alpha*A*x + beta*y
    
    trans = 'N'  ! No transpose
    call dgemv(&
      & trans, DIM, DIM, &
      & 1.0_real64, diffusion_matrix, &
      DIM, U, 1, 0.0_real64, vector_diffusion, 1)
  end subroutine eval_diffusion_at_u
end module mod_sde_coefficients

