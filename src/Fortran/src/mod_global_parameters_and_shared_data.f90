module mod_global_parameters_and_shared_data
    !
    ! This module wraps the main Fortran functionality to be called from C
    !! It contains the main function and the setup for the SDE solver.
    !! It also includes the logging functionality to write the results to a file.
    !! The module uses the MKL library for random number generation and Gaussian sampling.
    !! The module also includes the setup for the SDE solver, including the generation of
    !! observation times, eigenvalues, and the drift and diffusion matrices.
    
    
    use iso_fortran_env, only: int32, real64
    use mod_alloc
    implicit none
    
    ! Make all entities public by default
    public
    integer(int32), parameter :: Nx = 50
    integer(int32), parameter :: Ny = 50
    integer(int32), parameter :: DIM = Nx * Ny
    integer(int32), parameter :: SEED = 765431
    integer(int32), parameter :: nobs = 50000
    integer(int32), parameter :: NUM_GAUSSIAN_SUB_STEPS = 10
    real(real64), parameter :: PI = acos(-1.0_real64)
    real(real64), parameter :: theta = 1.0_real64
    real(real64), parameter :: beta = 0.1_real64
    real(real64), parameter :: gamma = 2.0_real64
    real(real64), parameter :: sigma = 0.0001_real64
    real(real64), parameter :: delta = 1.0e-5_real64
    real(real64), parameter :: L1 = 5.0_real64
    real(real64), parameter :: L2 = 5.0_real64
    !-----------------------------------------
    ! Vector variables (1D)
    real(real64), allocatable :: times(:)
    real(real64), allocatable :: eigenvalues(:)
    real(real64), allocatable :: lambdas(:)
    real(real64), allocatable :: b(:)      ! Vector of coefficients (no underscore for compiler compatibility)
    real(real64), allocatable :: hs(:)
    real(real64), allocatable :: u_zero(:)
    real(real64), allocatable :: Ls(:)
    real(real64), allocatable :: u(:)
    real(real64), allocatable :: AM(:)
    real(real64), allocatable :: diffusion_derivative(:)
    ! Flattened matrix for file I/O
    
    
    ! Matrix variables (2D)
    real(real64), allocatable :: lambdamatrix(:,:)
    real(real64), allocatable :: path(:,:)
    real(real64), allocatable :: a(:, :)
    real(real64), allocatable :: bmat(:,:)    ! Matrix version
    real(real64), allocatable :: driftmat(:,:)
    real(real64), allocatable :: diffusionmat(:,:)
    real(real64), allocatable :: brownian(:,:)
    real(real64), allocatable :: gaussian_samples(:,:)
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
        
        allocate(u_zero(DIM))
        u_zero = 0.0_real64
        open(99, file="../data/u0_proj_row.dat")
        read(99, *) u_zero
        close(99)
        
        allocate(Ls(DIM))
        Ls = 0.0_real64
        
        allocate(lambdas(DIM))
        lambdas = 0.0_real64
        
        allocate(diffusion_derivative(DIM))
        diffusion_derivative = 0.0_real64
        
        ! Allocate and initialize 2D arrays (matrices)
        allocate(path(0:nobs, DIM))
        path = 0.0_real64
        
        allocate(a(DIM, DIM))
        a = 0.0_real64
        
        ! Reshape AM into a(DIM,DIM)
        call reshape_matix_A()
        
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
        
        allocate(gaussian_samples(nobs * NUM_GAUSSIAN_SUB_STEPS, DIM))
        gaussian_samples = 0.0_real64
    end subroutine allocate_dynamic_memory
    
    subroutine reshape_matix_A()
        implicit none
        integer i, j, m, k, l, n, tot
        tot=0
        do i = 1,Nx
            do j = 1,Ny
                m = i + (j - 1) * Nx
                do k = 1, Nx
                    do l = 1,Ny
                        n = k + (l - 1) * Nx
                        tot=tot+1
                        a(m, n) = AM(tot)
                    enddo
                enddo
            enddo
        enddo
    end subroutine reshape_matix_A
    subroutine display_parameters()
        use iso_fortran_env, only: real64
        implicit none
        character(len=*), parameter :: fmt_int = "(A30, A3, I8)"
        character(len=*), parameter :: fmt_real = "(A30, A3, ES12.6)"
        print *, "============================================================"
        print *, "              Global Simulation Parameters"
        print *, "------------------------------------------------------------"
        ! Integer Parameters Table
        print *, ">> Integer Parameters"
        print *, "Parameter                      |   Value"
        print *, "------------------------------+----------------------------"
        write(*, fmt_int) "Nx", " : ", Nx
        write(*, fmt_int) "Ny", " : ", Ny
        write(*, fmt_int) "DIM", " : ", DIM
        write(*, fmt_int) "nobs", " : ", nobs
        write(*, fmt_int) "NUM_GAUSSIAN_SUB_STEPS", " : ", NUM_GAUSSIAN_SUB_STEPS
        write(*, fmt_int) "SEED", " : ", SEED
        
        print *, ""
        
        ! Real Parameters Table
        print *, ">> Real Parameters"
        print *, "Parameter                      |   Value"
        print *, "------------------------------+----------------------------"
        write(*, fmt_real) "L1", " : ", L1
        write(*, fmt_real) "L2", " : ", L2
        write(*, fmt_real) "delta", " : ", delta
        write(*, fmt_real) "theta", " : ", theta
        write(*, fmt_real) "beta", " : ", beta
        write(*, fmt_real) "gamma", " : ", gamma
        write(*, fmt_real) "sigma", " : ", sigma
        write(*, fmt_real) "PI", " : ", PI
        
        print *, "============================================================"
    end subroutine display_parameters
    
    !!> @brief Deallocates all allocatable arrays defined in the module.
    !!>
    !!> This subroutine checks whether each allocatable array is currently allocated,
    !!> and if so, it deallocates the array to release memory resources.
    !!>
    !!> @note This subroutine should be called before program termination or when restarting
    !!>       simulations to prevent memory leaks and ensure clean reuse of data structures.
    subroutine deallocate_all_shared_data()
        use mod_alloc, only: free_vector
        implicit none
        if (allocated(AM)) deallocate(AM)
        if (allocated(Ls)) deallocate(Ls)
        if (allocated(a)) deallocate(a)
        if (allocated(b)) deallocate(b)
        if (allocated(bmat)) deallocate(bmat)
        if (allocated(brownian)) deallocate(brownian)
        if (allocated(diffusion_derivative)) deallocate(diffusion_derivative)
        if (allocated(diffusionmat)) deallocate(diffusionmat)
        if (allocated(driftmat)) deallocate(driftmat)
        if (allocated(eigenvalues)) deallocate(eigenvalues)
        if (allocated(gaussian_samples)) deallocate(gaussian_samples)
        if (allocated(hs)) deallocate(hs)
        if (allocated(lambdamatrix)) deallocate(lambdamatrix)
        if (allocated(lambdas)) deallocate(lambdas)
        if (allocated(path)) deallocate(path)
        if (allocated(times)) deallocate(times)
        if (allocated(u)) deallocate(u)
        if (allocated(u_zero)) deallocate(u_zero)
    end subroutine deallocate_all_shared_data

end module mod_global_parameters_and_shared_data