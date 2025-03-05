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
  implicit none
contains

  !> @brief Prints a matrix to the console.
  !> 
  !> Prints a `rows x cols` matrix `A` with formatted floating-point numbers.
  !> Each row is printed separately.
  !>
  !> @param[in]  A     The matrix to print.
  !> @param[in]  rows  Number of rows in the matrix.
  !> @param[in]  cols  Number of columns in the matrix.
  subroutine print_matrix(A, rows, cols)
    implicit none
    integer, intent(in) :: rows, cols
    real, intent(in) :: A(rows, cols)
    integer :: i

    print *, "Matrix:"
    do i = 1, rows
        print '(100F8.3)', A(i, 1:cols)  ! Adjust format based on max column size
    end do

  end subroutine print_matrix


  !> @brief Prints a matrix with row and column indices.
  !> 
  !> This function prints a matrix `A` of size `rows x cols`, including the column indices at the top 
  !> and row indices on the left for better readability.
  !>
  !> @param[in]  A     The matrix to print.
  !> @param[in]  rows  Number of rows in the matrix.
  !> @param[in]  cols  Number of columns in the matrix.
  subroutine print_matrix_with_indices(A, rows, cols)
    implicit none
    integer(int32), intent(in) :: rows, cols
    real(real64), intent(in) :: A(rows, cols)
    integer :: i, j

    print *, "------------------------------------------------"
    print *, ""

    ! Print column headers
    write(*, '(A, *(I8))') " ", (j, j=1, cols)

    ! Print matrix with row indices
    do i = 1, rows
        write(*, '(I3, 100F8.4)') i, A(i, :)
    end do
    print *, "------------------------------------------------"
    print *, ""
  end subroutine print_matrix_with_indices


  !> @brief Writes a matrix to a file.
  !> 
  !> This subroutine writes a `rows x cols` matrix `A` to a specified file.
  !> The matrix is saved in a structured format where each row is written as a line.
  !>
  !> @param[in]  A        The matrix to write.
  !> @param[in]  rows     Number of rows in the matrix.
  !> @param[in]  cols     Number of columns in the matrix.
  !> @param[in]  filename Name of the file to save the matrix.
  subroutine write_matrix_to_file(A, rows, cols, filename)
    implicit none
    integer, intent(in) :: rows, cols
    real, intent(in) :: A(rows, cols)
    character(len=*), intent(in) :: filename
    integer :: i, unit_number

    ! Open a file for writing
    open(newunit=unit_number, file=filename, status='replace')

    ! Write the matrix
    do i = 1, rows
        write(unit_number, '(100F8.3)') A(i, :)
    end do

    ! Close the file
    close(unit_number)
  end subroutine write_matrix_to_file
  !> @brief Prints a vector to the console.
  !> 
  !> This subroutine prints a vector `V` of size `N`, with each element on a new line.
  !>
  !> @param[in] V The vector to print.
  !> @param[in] N The size of the vector.
  subroutine print_vector(V, N)
    implicit none
    integer(int32), intent(in) :: N
    real(real64), intent(in) :: V(N)
    integer(int32) :: i

    print *, "Vector:"
    do i = 1, N
        print '(F8.3)', V(i)  ! Print each element with 3 decimal places
    end do

  end subroutine print_vector
  !> @brief Prints a vector with indices.
  !> 
  !> Prints each element of vector `V` with its corresponding index.
  !>
  !> @param[in] V The vector to print.
  !> @param[in] N The size of the vector.
  subroutine print_vector_with_indices(V, N)
    implicit none
    integer(int32), intent(in) :: N
    real(real64), intent(in) :: V(N)
    integer :: i
    print*,""
    ! Print header
    print *, "  Index    Value"
    print *, "----------------"
    do i = 1, N
        print '(I6, F10.4)', i, V(i)  ! Print index and value
    end do
    print*,""
  end subroutine print_vector_with_indices
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
  pure subroutine gen_observation_times(nobs, delta, times)
    implicit none
    integer(int32), intent(in) :: nobs
    integer(int32) :: i
    real(real64), intent(in) :: delta
    real(real64), intent(out) :: times(0:nobs)
    
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
  !> @param[out] lambda_numbers Computed eigenvalues.
  pure subroutine gen_lambdas(DIM, Nx, Ny, L1, L2, lambda_numbers)
    implicit none
    integer(int32), intent(in) :: DIM, Nx, Ny
    real(real64), intent(in) :: L1, L2
    real(real64), intent(out):: lambda_numbers(DIM)
    
    real(real64), parameter :: PI = 2.D0 * DASIN(1.D0)
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
              lambda_numbers(m) = lambda_ij
            endif
          enddo
        enddo
      enddo
    enddo
    return
  end subroutine gen_lambdas

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
  pure subroutine MB(DIM, lambdas, gamma, B)
    implicit none
    integer(int32), intent(in) :: DIM
    real(real64), intent(in) :: gamma
    real(real64), intent(in) :: lambdas(DIM)
    real(real64), intent(out) :: B(DIM, DIM)
    integer(int32) i

    B(:,:) = 0.0
    do i=1, DIM
      B(i,i) = lambdas(i) ** (-gamma)
    enddo
    return
  end subroutine MB

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
  pure subroutine gen_matrix_diag_B(DIM, lambdas, gamma, B)
    implicit none
    integer(int32), intent(in) :: DIM
    real(real64), intent(in) :: gamma
    real(real64), intent(in) :: lambdas(DIM)
    real(real64), intent(out) :: B(DIM)
    integer(int32) i

    B(:) = 0.0
    
    do i=1, DIM
      B(i) = lambdas(i) ** (-gamma)
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
  pure subroutine gen_lambda_matrix(DIM, lambdas, lambda_matrix)
    implicit none
    integer(int32), intent(in) :: DIM
    real(real64), intent(in) :: lambdas(DIM)
    real(real64), intent(out) :: lambda_matrix(DIM, DIM)
    integer(int32) i

    lambda_matrix(:,:) = 0.0
    do i=1,DIM
      lambda_matrix(i,i) = lambdas(i)
    enddo
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
  pure subroutine MA(DIM, Nx, Ny, AM, A)
    implicit none
    integer(int32), intent(in) :: DIM, Nx, Ny
    real(real64), intent(in) :: AM(DIM * DIM)
    real(real64), intent(out) :: A(DIM, DIM)

    integer(int32) :: i, j, k, l, m, n, tot
    tot=1

    do i = 1, Nx
      do j = 1, Ny
        m = i + (j - 1) * Ny
          do k = 1, Nx
            do l = 1, Ny
                n = k + (l - 1) * Ny
                A(m, n) = AM(tot)
                tot = tot + 1
            enddo
          enddo
      enddo
    enddo
    return
  end subroutine MA

end module mod_par_generators
