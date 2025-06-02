module mod_global_parameters_and_shared_data
    use iso_fortran_env, only: int32, real64
    use mod_alloc
    implicit none
    
    ! Make all entities public by default
    public
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
    ! Vector variables (1D)
    real(real64), allocatable :: times(:)
    real(real64), allocatable :: eigenvalues(:)
    real(real64), allocatable :: lambdas(:)
    real(real64), allocatable :: b(:)      ! Vector of coefficients (no underscore for compiler compatibility)
    real(real64), allocatable :: hs(:)
    real(real64), allocatable :: startx(:)
    real(real64), allocatable :: Ls(:)
    real(real64), allocatable :: u(:)
    real(real64), allocatable :: AM(:)      ! Flattened matrix for file I/O
    
    ! Matrix variables (2D)
    real(real64), allocatable :: lambdamatrix(:,:)
    real(real64), allocatable :: path(:,:)
    real(real64), allocatable :: a(:,:)
    real(real64), allocatable :: bmat(:,:)    ! Matrix version
    real(real64), allocatable :: driftmat(:,:)
    real(real64), allocatable :: diffusionmat(:,:)
    real(real64), allocatable :: brownian(:,:)
               
    real(real64) :: mean_a, std_a, winner_delta
    
contains

    subroutine allocate_dynamic_memory()
        ! Load matrix A entries
        allocate(AM(DIM * DIM))  ! Dynamically allocate memory for AM
        open(99, file="../data/MatrixA.dat")
            read(99, *) AM
        close(99)
        
        ! Allocate and initialize 1D arrays (vectors)
        allocate(times(nobs))
        times = 0.0_real64
        
        allocate(eigenvalues(DIM))
        eigenvalues = 0.0_real64
        
        allocate(b(DIM))
        b = 0.0_real64
        
        allocate(u(DIM))
        u = 0.0_real64
        
        allocate(hs(DIM))
        hs = 0.0_real64
        
        allocate(startx(DIM))
        startx = 0.0_real64
        
        allocate(Ls(DIM))
        Ls = 0.0_real64
        
        allocate(lambdas(DIM))
        lambdas = 0.0_real64
        
        ! Allocate and initialize 2D arrays (matrices)
        allocate(path(0:nobs, DIM))
        path = 0.0_real64
        
        allocate(a(DIM, DIM))
        a = 0.0_real64
        
        ! Reshape AM into a(DIM,DIM)
        a = reshape(AM, [DIM, DIM])
        
        allocate(bmat(DIM, DIM))
        bmat = 0.0_real64
        
        allocate(driftmat(DIM, DIM))
        driftmat = 0.0_real64
        
        allocate(diffusionmat(DIM, DIM))
        diffusionmat = 0.0_real64
        
        allocate(lambdamatrix(DIM, DIM))
        lambdamatrix = 0.0_real64
        
        allocate(brownian(nobs, DIM))
        brownian = 0.0_real64
        
    end subroutine allocate_dynamic_memory
end module mod_global_parameters_and_shared_data