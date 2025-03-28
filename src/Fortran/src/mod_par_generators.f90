!> @ingroup modules
!> @author F. Baltazar-Larios, S. Diaz-Infante
!> @brief This module generates the parameters of sde
!> 
!>  @f$
!>    dU = -(\beta \Lambda  + \theta A) \ U dt + \sigma B\ U dW(t) 
!>  @f$
!> conforming to a spectral decomposition.
module mod_par_generators
  use iso_fortran_env, only: int32, real64
  use mod_global_parameters_and_shared_data
  implicit none
contains

!> @brief Generates an array of observation times.
  !> 
  !> Creates a stencil of length `nobs`, with step-size `delta`.
  !> The observation times are computed as:
  !> @f[
  !> \text{times}(i) = i \cdot \delta, \quad i \in \{0, ..., nobs\}
  !> @f]
  !>
  !> @param[in]  nobs  The number of observations in the time interval.
  !> @param[in]  delta Step-size, such that @f$ nobs \cdot \delta = T @f$.
  !> @param[out] times Array of observation times.
  subroutine gen_observation_times()
    implicit none
    integer(int32) :: i
    do i=0, nobs
      times(i) = i * delta
    end do
    return
  end subroutine gen_observation_times

  !> @brief Generates eigenvalues for a differential operator.
  !> 
  !> Computes the diagonal elements of the matrix @f$ \Lambda @f$ using:
  !> @f[
  !> \lambda_{ij} = \pi^2 \left( \frac{i}{L_1} \right)^2 + \pi^2 \left( \frac{j}{L_2} \right)^2
  !> @f]
  !> where @f$ i @f$ and @f$ j @f$ are grid indices.
  !>
  !> @param[in]  DIM  Number of elements.
  !> @param[in]  Nx   Number of eigen basic vectors in the x-direction.
  !> @param[in]  Ny   Number of eigen basic vectors in the y-direction.
  !> @param[in]  L1   Length in the x-direction.
  !> @param[in]  L2   Length in the y-direction.
  !> @param[out] eigen_values Computed eigenvalues.
  subroutine gen_eigen_values()
    implicit none
    real(real64) pi_square, L1_res, L2_res, lambda_ij
    integer(int32) i, j, k, m, l, n

    pi_square = PI ** 2
    L1_res = (L1) ** (-1)
    L2_res = (L2) ** (-1)

    k=0
    do i=1, Nx
      do j=1, Ny
        m = i + (j - 1) * Ny
        do k=1, Nx
          do l=1, Ny
            n  = k + (l - 1) * Ny
            if (m == n) then
              lambda_ij = pi_square * ((i/L1)**2 + (j/L2)**2)
              eigen_values(m) = lambda_ij
            endif
          end do
        end do
      end do
    end do
    return
  end subroutine gen_eigen_values

  !> @brief Generates the matrix @f$ B @f$ based on eigenvalues.
  !> 
  !> Computes the diagonal matrix @f$ B @f$ using:
  !> @f[
  !> B_{ii} = \lambda_i^{-\gamma}
  !> @f]
  !> except for the first entry, which is set to 1.
  !>
  !> @param[in]  DIM      Dimension of the matrix.
  !> @param[in]  lambdas  Array of eigenvalues.
  !> @param[in]  gamma    Power exponent.
  !> @param[out] B        The computed diagonal matrix.
  subroutine build_matrix_B(eigen_values, matrix_B)
    implicit none
    real(real64), intent(in) :: eigen_values(DIM)
    real(real64), intent(inout) :: matrix_B(DIM, DIM)
    integer(int32) i
    do i=1, DIM
      matrix_B(i, i) = eigen_values(i) ** (-gamma)
    end do
    return
  end subroutine build_matrix_B

  !> @brief Generates the elements of the diagonal matrix @f$ B @f$ based on eigenvalues @f$ \lambda_{\mathbf{k}} @f$.
  !> 
  !> Computes the diagonal matrix @f$ B @f$ using:
  !> @f[
  !> B_{ii} = \lambda_i^{-\gamma}
  !> @f]
  !> except for the first entry, which is set to 1.
  !>
  !> @param[in]  DIM      Dimension of the matrix.
  !> @param[in]  lambdas  Array of eigenvalues.
  !> @param[in]  gamma    Power exponent.
  !> @param[out] B        The computed diagonal matrix.
  subroutine gen_matrix_diag_B()
    implicit none
    integer(int32) i

    B_(:) = 0.0
    do i=1, DIM
      B_(i) = eigen_values(i) ** (-gamma)
    end do
    return
  end subroutine gen_matrix_diag_B

  !> @brief Generates a diagonal matrix from a vector of eigenvalues.
  !>
  !> Converts the eigenvalue vector `lambdas` into a diagonal matrix.
  !> 
  !> @param[in]  DIM          Number of elements.
  !> @param[in]  lambdas      Eigenvalue vector.
  !> @param[out] lambda_matrix Resulting diagonal matrix.
  subroutine gen_lambda_matrix()
    implicit none
  
    integer(int32) i

    do i=1,DIM
      lambda_matrix(i, i) = eigen_values(i) 
    end do
    return
  end subroutine gen_lambda_matrix

  !> @brief Constructs the interaction matrix A.
  !>
  !> Rearranges the vectorized matrix `AM` into a 2D matrix `A` with 
  !> dimensions determined by `Nx` and `Ny`.
  !>
  !> @param[in]  DIM  Dimension of the square matrix.
  !> @param[in]  Nx   Number of eigen base vectors in x-direction.
  !> @param[in]  Ny   Number of eigen base vectors in y-direction.
  !> @param[in]  AM   Flattened matrix data.
  !> @param[out] A    Reshaped 2D matrix.
  
  subroutine assemble_matrix_A()
    implicit none

    integer(int32) :: i, j, k, l, m, n, idx
    idx=1

    do i = 1, Nx
      do j = 1, Ny
        m = i + (j - 1) * Ny
          do k = 1, Nx
            do l = 1, Ny
                n = k + (l - 1) * Ny
                A(m, n) = AM(idx)
                idx = idx + 1
            end do
          end do
      end do
    end do
    return
  end subroutine assemble_matrix_A

end module mod_par_generators
