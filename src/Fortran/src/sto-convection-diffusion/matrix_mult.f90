program matrix_multiplication
    implicit none
    integer, parameter :: n = 3  ! Matrix size (can be changed)
    double precision, dimension(n, n) :: A, B, C
    double precision :: alpha, beta
    integer :: lda, ldb, ldc

    ! Initialize matrices A and B
    A = reshape((/1.0d0, 2.0d0, 3.0d0, &
                  4.0d0, 5.0d0, 6.0d0, &
                  7.0d0, 8.0d0, 9.0d0/), shape(A))

    B = reshape((/9.0d0, 8.0d0, 7.0d0, &
                  6.0d0, 5.0d0, 4.0d0, &
                  3.0d0, 2.0d0, 1.0d0/), shape(B))

    ! Initialize result matrix C to zero
    C = 0.0d0

    ! Scaling factors
    alpha = 1.0d0
    beta = 0.0d0

    ! Leading dimensions (typically the row size)
    lda = n
    ldb = n
    ldc = n

    ! Call DGEMM to perform C = alpha*A*B + beta*C
    call dgemm('N', 'N', n, n, n, alpha, A, lda, B, ldb, beta, C, ldc)

    ! Print the result
    print *, "Resultant Matrix C:"
    print *, C

end program matrix_multiplication
