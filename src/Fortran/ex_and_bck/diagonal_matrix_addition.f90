program diagonal_matrix_addition

    implicit none
    include 'mkl_blas.fi'  ! Use MKL BLAS interfaces
    integer, parameter :: N = 4
    real(8), dimension(N) :: D ! Diagonal elements
    real(8), dimension(N, N) :: A ! Full matrix
    integer :: i

    ! Initialize diagonal matrix as a vector
    D = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0 /)

    ! Initialize A
    A = reshape([1.0d0, 2.0d0, 3.0d0, 4.0d0, &
                 5.0d0, 6.0d0, 7.0d0, 8.0d0, &
                 9.0d0, 10.0d0, 11.0d0, 12.0d0, &
                 13.0d0, 14.0d0, 15.0d0, 16.0d0], [N, N])

    ! Efficiently add the diagonal elements to the matrix
    do i = 1, N
        A(i, i) = A(i, i) + D(i)
    end do

    ! Print the result
    print *, "Updated Matrix A:"
    do i = 1, N
        print *, A(i, :)
    end do

end program diagonal_matrix_addition
