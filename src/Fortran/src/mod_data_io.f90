

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
    integer, parameter :: col_width = 15  ! Match ES13.5 format
    character(len=6) :: row_label
    print *, "(++++) ", trim(name)
    print *, repeat("-", 4 + cols * col_width)
    
    ! Print column headers
    write(*, '(A)', advance="no") repeat(" ", 6)  ! Space for row index
    do j = 1, cols
      write(*, '(I15)', advance="no") j
    end do
    print *  ! Newline after headers
    
    ! Print matrix with row indices and aligned entries
    do i = 1, rows
      write(*, '(I6)', advance="no") i
      write(*, '( *(ES15.6) )') A(i, 1:cols)
    end do
    
    print *, repeat("=", 4 + cols * col_width)
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
        print '(I6, A1, ES15.7E2)', i, char(9), V(i)  ! Print index and value
    end do
    print*,""
  end subroutine print_vector_with_indices

!> @brief Adds a column of row indices (starting from 1) to the input array.
!! @param[in] input A 2D real(real64) array of shape (n, m)
!! @param[out] output A 2D real(real64) array of shape (n, m+1) with index in first column

  subroutine add_index_column(input, output)

    use iso_fortran_env, only: real64
    implicit none
  
    real(real64), intent(in)  :: input(:,:)
    real(real64), intent(out) :: output(:,:)
    integer :: n, m, i
    n = size(input, 1)
    m = size(input, 2)
    do i = 1, n
      output(i, 1) = real(i, kind=real64)       ! Row index (starting from 1)
      output(i, 2:m+1) = input(i, :)            ! Copy row from input
    end do
  end subroutine add_index_column
  
!> @brief Save a 2D real array to a CSV file with a header row.
!! 
!! This subroutine writes the contents of a 2D real array to a CSV file, 
!! including a header line with column names.
!!
!! @param[in] filename Name of the CSV output file.
!! @param[in] array 2D real array of size (rows, cols) to write to file.
!! @param[in] rows Number of rows in the array.
!! @param[in] cols Number of columns in the array.
!! @param[in] header Array of column names (length must equal cols).
!!
!! @note Each value is written with 4 decimal places (F10.4 format).
!! @warning Header size must match the number of columns.
  
  subroutine save_array_to_csv_with_header(filename, array, rows, cols, header)
    implicit none
  
    character(len=*), intent(in) :: filename
    real(real64), intent(in) :: array(:,:)
    integer(int32), intent(in) :: rows, cols
    character(len=*), dimension(:), intent(in) :: header
    integer :: i, j
    integer :: unit
  
    ! Open a unit number and the file for writing
    open(newunit=unit, file=filename, status='replace', action='write')
  
    ! Write the header
    do j = 1, cols
      write(unit, '(A)', advance='no') trim(header(j))
      if (j < cols) then
        write(unit, '(A)', advance='no') ','
      else
        write(unit, *) ''
      end if
    end do
  
    ! Write the data
    do i = 1, rows
      do j = 1, cols
        if(j < 2) then
          write(unit, '(I6)', advance='no') int(array(i, j), kind=Int32)
        else
          write(unit, '(F12.8)', advance='no') array(i, j)
        end if
        
        if (j < cols) then
          write(unit, '(A)', advance='no') ','
        else
          write(unit, *) ''
        end if
      end do
    end do
  
    close(unit)
  end subroutine save_array_to_csv_with_header
  
  !> @brief Writes a vector to a file.
  !> 
  !> This subroutine writes a vector to a specified file.
  !> Each element is written in a new line.
  !>
  !> @param[in]  filename Name of the file to save the vector.
  !> @param[in]  vector   The vector to write.
  subroutine write_data(filename, vector)
    implicit none
    character(len=*), intent(in) :: filename
    real(real64), intent(in) :: vector(:)
    integer :: i, unit_number
    
    ! Open a file for writing
    open(newunit=unit_number, file=filename, status='replace')
    
    ! Write each element of the vector on a separate line
    do i = 1, size(vector)
        write(unit_number, '(F12.6)') vector(i)
    end do
    
    ! Close the file
    close(unit_number)
  end subroutine write_data

  subroutine save_real64_2d_array_to_binary(filename, array)
    use iso_fortran_env, only: real64
    implicit none
    
    ! Arguments
    character(len=*), intent(in) :: filename
    real(real64), intent(in) :: array(:,:)
    integer(int32) :: nobs, DIM
    
    ! Local
    integer :: unit
    nobs = size(array, 1)
    DIM = size(array, 2)
    ! Open file in stream (raw binary) mode and write the 2D array
    open(newunit=unit, file=filename, access="stream", form="unformatted", &
            status="replace", action="write")
    write(unit) nobs, DIM ! write shape first
    write(unit) array       ! then write data
    close(unit)
  
  end subroutine save_real64_2d_array_to_binary
  !> Find all NaN entries in a 2D real(real64) array.
  !!
  !! This subroutine scans a two-dimensional array of real(real64) values
  !! and returns the row and column indices of all elements that are NaN
  !! (Not a Number), according to the IEEE standard.
  !!
  !! The output arrays `row_nan` and `col_nan` contain the row and column
  !! indices, respectively, for each detected NaN entry.
  !!
  !! @param[in]  A         The input 2D array of real(real64) values
  !! @param[in]  rows      Number of rows in the array A
  !! @param[in]  cols      Number of columns in the array A
  !! @param[out] row_nan   Integer array containing the row indices of NaNs
  !! @param[out] col_nan   Integer array containing the column indices of NaNs
  !! @param[out] count_nan The number of NaN entries found in the array
  !!
  !! Example:
  !!   A = reshape([1.0, 0.0/0.0, 3.0, 4.0], [2,2])
  !!   call find_nan_indices_2d(A, 2, 2, row_nan, col_nan, count_nan)
  !!   ! Now: count_nan = 1, row_nan(1) = 1, col_nan(1) = 2
  !!
  subroutine find_nan_indices_2d(A, rows, cols, row_nan, col_nan, count_nan)
    use iso_fortran_env, only: real64, int32
    use ieee_arithmetic, only: ieee_is_nan
    implicit none
    
    ! Arguments
    integer(int32), intent(in) :: rows, cols
    real(real64), intent(in) :: A(rows, cols)
    integer(int32), allocatable, intent(out) :: row_nan(:), col_nan(:)
    integer(int32), intent(out) :: count_nan
    
    ! Locals
    integer :: i, j, k
    
    ! First pass: count NaNs
    count_nan = count(ieee_is_nan(A))
    
    allocate(row_nan(count_nan), col_nan(count_nan))
    
    ! Second pass: collect indices
    k = 0
    do i = 1, rows
      do j = 1, cols
        if (ieee_is_nan(A(i,j))) then
          k = k + 1
          row_nan(k) = i
          col_nan(k) = j
        end if
      end do
    end do
  end subroutine find_nan_indices_2d
end module mod_data_io
