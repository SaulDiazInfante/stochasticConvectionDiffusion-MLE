!> @ingroup modules
!> @author F. Baltazar-Larios, S. Diaz-Infante
!> @brief This module computes the coefficients of a SDE regarding
!> the PDE model.
module mod_sde_coefficients
    use iso_fortran_env, only : int32, real64
    use mod_alloc
    use mod_global_parameters_and_shared_data
    implicit none
    ! Include the MKL module
    !include 'mkl_blas.fi'
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
    subroutine gen_drift_matrix()
        implicit none
        integer(int32) :: i
        driftmat(:, :) = theta * a(:, :)
        do i = 1, DIM
            driftmat(i, i) = driftmat(i, i) + beta * lambdas(i)
        end do
        driftmat = -1.0_real64 * driftmat
        return
    end subroutine gen_drift_matrix
    
    !> @brief Generates the diffusion matrix for a stochastic differential equation.
    !> @param[in]  DIM           Dimension of the square matrices.
    !> @param[in]  sigma         Noise intensity parameter.
    !> @param[in]  B             A matrix of size (DIM, DIM) influencing the
    !>                           diffusion term.
    !> @param[out] diffusion_mat The resulting diffusion matrix of size (DIM, DIM).
    
    subroutine gen_diffusion_matrix()
        implicit none
        diffusionmat(:, :) = sigma * bmat(:, :)
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
    subroutine eval_whole_drift(vector_U, vector_drift)
        implicit none
        real(real64), intent(in) :: vector_U(DIM)
        real(real64), allocatable, intent(out) :: vector_drift(:)
        character(1) :: trans
        real(real64), allocatable :: temp_diagonal(:), temp_matrix(:, :)
        integer(int32) :: i, j
        call alloc_vector(vector_drift, DIM)
        call alloc_vector(temp_diagonal, DIM)
        call alloc_array(temp_matrix, DIM, DIM)
        trans = 'N'  ! No transpose
        
        temp_diagonal(:) = beta * lambdas(:)
        temp_matrix(:, :) = theta * a(:, :)
        
        do i = 1, DIM
            temp_matrix(i, i) = temp_matrix(i, i) + temp_diagonal(i)
        end do
        vector_drift = vector_drift - MATMUL(temp_matrix, vector_U)
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
            & vector_U, &
            & vector_drift &
            &)
        implicit none
        real(real64), intent(in) :: vector_U(DIM)
        real(real64), allocatable, intent(out) :: vector_drift(:)
        
        integer(int32) :: i, j
        
        call alloc_vector(vector_drift, DIM)
        vector_drift = vector_drift + MATMUL(driftmat, vector_U)
    end subroutine eval_drift_at_u
    
    !> @brief compute the diffusion coefficient of SDE equation
    !!  @f$ dU = -(\beta \Lambda  + \theta A) U dt + \sigma B U dW(t) @f$.
    !!
    !! @param[in]   DIM           int32 matrix dimension
    !! @param[in]   sigma real64
    !! @param[in]   diffusion_matrix  real64(DIM, DIM) current values in the array
    !! @param[in]   U                 real64(DIM)
    !! @param[out]  vector_diffusion  real64(DIM  \f $\sigma B U \f$.
    
    subroutine eval_diffusion_at_u(vector_U, vector_diffusion)
        implicit none
        real(real64), intent(in) :: vector_U(DIM)
        real(real64), allocatable, intent(out) :: vector_diffusion(:)
        call alloc_vector(vector_diffusion, DIM)
        
        vector_diffusion = vector_diffusion + MATMUL(diffusionmat, vector_U)
    end subroutine eval_diffusion_at_u
    
    
    !> @brief compute the diffusion coefficient of SDE equation
    !!  @f$ dU = -(\beta \Lambda  + \theta A) U dt + \sigma B U dW(t) @f$.
    !!
    !! @param[in]   DIM           int32 matrix dimension
    !! @param[in]   sigma real64
    !! @param[in]   diffusion_matrix  real64(DIM, DIM) current values in the array
    !! @param[in]   U                 real64(DIM)
    !! @param[out]  vector_diffusion  real64(DIM  \f $\sigma B U \f$.
    
    subroutine eval_diagonal_diffusion_at_u(vector_U, vector_diffusion)
        implicit none
        real(real64), intent(in) :: vector_U(DIM)
        real(real64), allocatable, intent(out) :: vector_diffusion(:)
        call alloc_vector(vector_diffusion, DIM)
        vector_diffusion = sigma * b(:) * vector_U(:)
    end subroutine eval_diagonal_diffusion_at_u
    
    !> @brief Compute the derivative of the diffusion coefficient for the SDE equation
    !!  @f$ dU = -(\beta \Lambda  + \theta A) U dt + \sigma B U dW(t) @f$
    !!
    !! @param[out] diffusion_derivative real64(DIM) The derivative of diffusion term @f$ \sigma B @f$
    subroutine compute_diffusion_derivative()
        implicit none
        diffusion_derivative = sigma * b(:)
    end subroutine compute_diffusion_derivative
    
    !> @brief Evaluate the derivative of the diffusion coefficient for the SDE equation
    !!  @f$ dU = -(\beta \Lambda  + \theta A) U dt + \sigma B U dW(t) @f$
    !!
    !! @param[out] vector_diffusion_derivative real64(DIM) The derivative of diffusion term @f$ \sigma B @f$
    subroutine eval_diffusion_derivative(vector_diffusion_derivative)
        implicit none
        real(real64), allocatable, intent(out) :: vector_diffusion_derivative(:)
        call alloc_vector(vector_diffusion_derivative, DIM)
        vector_diffusion_derivative = sigma * b(:)
    end subroutine eval_diffusion_derivative
end module mod_sde_coefficients