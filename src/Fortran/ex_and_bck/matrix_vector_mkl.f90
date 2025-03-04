program matrix_vector_mkl
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none

    ! Include the MKL module
    include 'mkl_blas.fi'

    ! Variable declarations
    integer(int32) :: m, n, i, j, info
    integer(int32) :: lda
    real(real64), allocatable :: A(:,:), x(:), y(:)
    real(real64) :: alpha, beta
    character(1) :: trans

    ! Matrix and vector dimensions
    m = 3  ! Number of rows in matrix A
    n = 3  ! Number of columns in matrix A
    lda = m ! Leading dimension of A

    ! Allocate memory for A, x, and y
    allocate(A(lda,n), x(n), y(m), stat=info)
    if (info /= 0) then
        write(*,*) "Error: Memory allocation failed!"
        stop
    end if

    ! Initialize matrix A and vector x with some values
    ! A is an m x n matrix
    A = reshape([ 1.0_real64, 4.0_real64, 7.0_real64, &
                  2.0_real64, 5.0_real64, 8.0_real64, &
                  3.0_real64, 6.0_real64, 9.0_real64 ], [m, n])

    ! x is a vector of length n
    x = [1.0_real64, 2.0_real64, 3.0_real64]

    ! Initialize y with zeros
    y = 0.0_real64

    ! Print input matrix and vector
    write(*,*) "Matrix A:"
    do i = 1, m
        write(*,*) (A(i,j), j = 1, n)
    end do

    write(*,*) "Vector x:"
    write(*,*) x

    ! Set parameters for DGEMV
    ! y = alpha*A*x + beta*y
    alpha = 1.0_real64
    beta = 0.0_real64
    trans = 'N'  ! No transpose

    ! Perform matrix-vector multiplication using MKL: y = alpha*A*x + beta*y
    call DGEMV(trans, m, n, alpha, A, lda, x, 1, beta, y, 1)

    ! Print the result
    write(*,*) "Result vector y = A*x:"
    write(*,*) y

    ! Expected result for this example:
    ! y(1) = 1*1 + 2*2 + 3*3 = 14
    ! y(2) = 4*1 + 5*2 + 6*3 = 32
    ! y(3) = 7*1 + 8*2 + 9*3 = 50
    write(*,*) "Expected result:"
    write(*,*) [14.0_real64, 32.0_real64, 50.0_real64]

    ! Clean up
    deallocate(A, x, y)

    write(*,*) "Matrix-vector multiplication completed successfully!"

end program matrix_vector_mkl

