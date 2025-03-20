module mod_global_parameters_and_shared_data
    use iso_fortran_env, only: int32, real64
    use mod_alloc
    implicit none
    integer(int32), parameter :: Nx = 10
    integer(int32), parameter :: Ny = 10
    integer(int32), parameter :: DIM = Nx * Ny
    integer(int32), parameter :: SEED = 765431
    integer(int32), parameter :: nobs = 1000
    real(real64), parameter :: PI = 2.D0 * DASIN(1.D0)
    real(real64), parameter :: theta = 0.5_real64
    real(real64), parameter :: beta = 0.5_real64
    real(real64), parameter :: gamma = 1.0_real64
    real(real64), parameter :: sigma = 0.2_real64
    real(real64), parameter :: delta = 0.0001_real64
    real(real64), parameter :: L1 = 5.0_real64
    real(real64), parameter :: L2 = 5.0_real64
!-----------------------------------------    
    real(real64), allocatable :: times(:)
    real(real64), allocatable :: eigen_values(:), lambdas(:), B_(:)
    real(real64), allocatable :: hs(:), startx(:), Ls(:)
    real(real64), allocatable :: U(:), U_(:), vector_drift(:)
    real(real64), allocatable :: AM(:)  ! Ensure AM is declared as allocatable
    real(real64), allocatable :: vector_diffusion(:)
    real(real64), allocatable :: vectorial_winner_delta(:)
    real(real64), allocatable :: initial_vector_winner(:)
    real(real64), allocatable :: gaussian_sample(:)
    
    real(real64), allocatable :: lambda_matrix(:,:)
    real(real64), allocatable :: path(:,:)
    real(real64), allocatable :: B(:, :)
    real(real64), allocatable :: drift_mat(:,:), diffusion_mat(:,:)
               
    real(real64), allocatable :: brownian(:,:)
    real(real64) :: mean_a, std_a, winner_delta
    
contains

    subroutine allocate_dynamic_memory()
        ! Load matrix A entries
        allocate(AM(DIM * DIM))  ! Dynamically allocate memory for AM
        open(99, file="../MatrixA.dat")
            read(99, *) AM
        close(99)
        call alloc_vector(times, nobs)
        call alloc_vector(eigen_values, DIM)
        call alloc_vector(U, DIM)
        call alloc_vector(U_, DIM)
        call alloc_vector(vector_drift, DIM)
        call alloc_vector(vector_diffusion, DIM)
        call alloc_vector(vectorial_winner_delta, DIM)
        call alloc_vector(initial_vector_winner, DIM)
        allocate(path(0:nobs, DIM))

        call alloc_array(B, DIM, DIM)
        call alloc_array(drift_mat, DIM, DIM)
        call alloc_array(diffusion_mat, DIM, DIM)
        call alloc_array(lambda_matrix, DIM, DIM)
        call alloc_array(brownian, nobs, DIM)
    end subroutine allocate_dynamic_memory
end module mod_global_parameters_and_shared_data