

!> @brief Module for data input/output operations.
!> 
!> This module provides subroutines for printing and writing matrices and vectors.
!> It includes functionality to print matrices and vectors with or without indices,
!> and to write matrices to files in a structured format.
!>
!> @details
!> The module includes the following subroutines:
!> - `print_matrix`: Prints a matrix to the console with formatted floating-point numbers.
!> - `print_matrix_with_indices`: Prints a matrix with row and column indices for better readability.
!> - `write_matrix_to_file`: Writes a matrix to a specified file in a structured format.
!> - `print_vector`: Prints a vector to the console with each element on a new line.
!> - `print_vector_with_indices`: Prints a vector with its corresponding indices for better readability.
!>
!> @note
!> - The matrix and vector elements are printed with specified formatting for floating-point numbers.
!> - The file writing operation overwrites the file if it already exists.
module mod_data_io
    use iso_fortran_env, only: int32, real64
    implicit none

!    public :: read_matrix_from_file, write_matrix_to_file
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
  !> @param[in]  name  string, label for matrix A.
  !> @param[in]  A     The matrix to print.
  !> @param[in]  rows  Number of rows in the matrix.
  !> @param[in]  cols  Number of columns in the matrix.
  subroutine print_matrix_with_indices(name, A, rows, cols)
    implicit none
    character(len=*), intent(in) :: name
    integer(int32), intent(in) :: rows, cols
    real(real64), intent(in) :: A(rows, cols)
    integer :: i, j

    print *, "(++++) ", name
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
  subroutine print_vector_with_indices(name, V, N)
    implicit none
    character(len=*), intent(in) :: name
    integer(int32), intent(in) :: N
    real(real64), intent(in) :: V(N)
    integer :: i
    print*,"(++++) vector ", name
    ! Print header
    print *, "  Index    Value"
    print *, "----------------"
    do i = 1, N
        print '(I6, F10.4)', i, V(i)  ! Print index and value
    end do
    print*,""
  end subroutine print_vector_with_indices
  
end module mod_data_io