!> @ingroup modules
!> @author F. Baltazar-Larios, S. Diaz-Infante
!> @brief This module implements the numerical solution of a sde
!>  Accordingly with the Milstein scheme
!> @f$$
!>  X_{n+1} ^ {\Delta} = f(X_n) \Delta + g(X_n) \Delta W_{n+1}
!>  + (Milstein corr)
!! see Kloeden & Platten (1994)
 module mod_sde_solver
  use iso_fortran_env, only: int32, real64
  use ieee_arithmetic, only: ieee_is_nan
  use mod_random_number_generator
  use mod_par_generators
  use mod_sde_coefficients
  use mod_data_io
  implicit none
  contains

  !> @brief Computes the increment of a Wiener process over a given time interval.
  !!
  !! This subroutine calculates the increment of a Wiener process (also known as Brownian motion)
  !! over a specified time interval `delta`. The increment is computed using a Gaussian sampler.
  !!
  !! @param[in] delta The time interval over which the Wiener process increment is computed.
  !! @param[in] NUM_GAUSSIAN_SUB_STEPS The number of subintervals to divide `delta` into for the Gaussian sampling.
  !! @param[in] user_seed (optional) An optional seed for the random number generator.
  !! @param[in] winner_0 The initial value of the Wiener process.
  !! @param[out] winner_delta The computed increment of the Wiener process over the interval `delta`.
  !!
  !! The subroutine uses the MKL Gaussian sampler to generate random samples from a normal distribution
  !! with mean 0 and standard deviation 1. These samples are then scaled by the square root of the
  !! subinterval length and summed to produce the Wiener process increment.
  subroutine scalar_winner_increment(winner_0, winner_delta, user_seed)
    implicit none
    integer, intent(in), optional :: user_seed
    real(real64), intent(in) :: winner_0
    real(real64), intent(out) :: winner_delta
    
    real(real64) dd
    real(real64), allocatable :: ddW(:)

    dd = delta / NUM_GAUSSIAN_SUB_STEPS
    if (present(user_seed)) then
      call mkl_gaussian_sampler(NUM_GAUSSIAN_SUB_STEPS, 0.0_real64, 1.0_real64, ddW, user_seed)
    else
      call mkl_gaussian_sampler(NUM_GAUSSIAN_SUB_STEPS, 0.0_real64, 1.0_real64, ddW)
    end if
    ddW = sqrt(dd) * ddW
    winner_delta = winner_0 + sum(ddW(:))
    return
  end subroutine scalar_winner_increment
  
  
  !> @brief Computes the increment of a vectorial Wiener process.
  !>
  !> This subroutine calculates the increment of a vectorial Wiener process
  !> over a given time step `delta`. The Wiener process is represented by
  !> `winner_0` and the increment is stored in `winner_delta`.
  !>
  !> @param delta The time step over which the increment is computed.
  !> @param DIM The dimension of the Wiener process.
  !> @param NUM_GAUSSIAN_SUB_STEPS The number of sub-steps within the time step `delta`.
  !> @param user_seed (Optional) Seed for the random number generator.
  !> @param winner_0 The initial value of the Wiener process.
  !> @param winner_delta The computed increment of the Wiener process.
  !>
  !> The subroutine uses the MKL library to sample from a Gaussian distribution
  !> and scales the samples appropriately to compute the Wiener process increment.
  subroutine vectorial_winner_increment(&
          &winner_0, &
          &winner_delta, &
          &user_seed)
    implicit none
    integer, intent(in), optional :: user_seed
    real(real64), intent(in) :: winner_0(DIM)
    real(real64), intent(out) :: winner_delta(DIM)
    
    real(real64) dd
    real(real64), allocatable, dimension(:, :) :: ddW
    integer(int32) i

    dd = delta / NUM_GAUSSIAN_SUB_STEPS
    if (present(user_seed)) then
      call mkl_array_gaussian_sampler(DIM, NUM_GAUSSIAN_SUB_STEPS, 0.0_real64, 1.0_real64, ddW, user_seed)
    else
      call mkl_array_gaussian_sampler(DIM, NUM_GAUSSIAN_SUB_STEPS, 0.0_real64, 1.0_real64, ddW)
    end if
    ddW = sqrt(dd) * ddW
    do i=1, DIM
      winner_delta(i) = winner_0(i) + sum(ddW(i, :))
    end do
    
    ! Deallocate local array to prevent memory leaks
    if (allocated(ddW)) deallocate(ddW)
    
    return
  end subroutine vectorial_winner_increment
 
!---------------------------------------------------------------------------
  !> @brief Performs a single Brownian motion step in a stochastic differential
  !> equation solver.
  !>
  !> This subroutine updates the position vector `endx` by performing a single
  !> Brownian motion step starting from the position vector `startx`. The step
  !> size is determined by the parameter `delta`.
  !>
  !> @param[in] DIM The dimension of the position vectors.
  !> @param[in] delta The time step size for the Brownian motion.
  !> @param[in] NUM_GAUSSIAN_SUB_STEPS The number of random variables (not used in current implementation).
  !> @param[in] startx The starting position vector of dimension `DIM`.
  !> @param[inout] endx The updated position vector of dimension `DIM`.
!---------------------------------------------------------------------------
  subroutine BrownianStep(startx, endx)
    implicit none
    real(real64), intent(in) :: startx(DIM)
    real(real64), intent(inout) :: endx(DIM)
    !!
    real(real64) xi, dd, dW(DIM)
    integer(int32) i
    
    !! TODO: Call mkl rng to create a vectorized version of the below code.
    call vectorial_winner_increment(startx, dW)
    endx(:) = startx(:) + dW(:)
    return
  end subroutine BrownianStep
  
  
  
!---------------------------------------------------------------------------
  !> @brief Performs a single Milstein step for a stochastic differential equation (SDE).
  !>
  !> @param[in] DIM The dimension of the SDE system.
  !> @param[in] delta The time increment for the step.
  !> @param[in] startx The starting values of the SDE system.
  !> @param[in] alpha The drift coefficients of the SDE system.
  !> @param[in] sigma The diffusion coefficients of the SDE system.
  !> @param[inout] endx The resulting values of the SDE system after the step.
  !> @param[inout] brown The Brownian motion increments used in the step.
  !>
  !> This subroutine performs a single step of the Milstein method for solving
  !> stochastic differential equations. It updates the state of the system
  !> from `startx` to `endx` using the given drift (`alpha`) and diffusion
  !> (`sigma`) coefficients, and the Brownian motion increments (`brown`).

  subroutine milstein_step(&
    &current_u, &
    &brownian_increment, &
    &next_u &
  &)
    implicit none
    real(real64), intent(in) :: current_u(DIM)
    real(real64), intent(in) :: brownian_increment(DIM)
    real(real64), allocatable, intent(out) :: next_u(:)
    
    real(real64), allocatable :: u_drift (:), u_diffusion(:)
    real(real64), allocatable :: u_euler_maruyama(:), u_milstein_correction(:)
  
    call alloc_vector(u_drift, DIM)
    call alloc_vector(u_diffusion, DIM)
    call alloc_vector(u_euler_maruyama, DIM)
    call alloc_vector(next_u, DIM)
    call alloc_vector(u_milstein_correction, DIM)
    call eval_drift_at_u(current_u, u_drift)
    call eval_diagonal_diffusion_at_u(current_u, u_diffusion)
    
    u_euler_maruyama(:) = current_u(:) &
      + u_drift(:) * delta &
      + u_diffusion(:) * brownian_increment(:)
    call compute_milstein_correction(brownian_increment, u_milstein_correction)
    next_u(:) = u_euler_maruyama(:) + u_milstein_correction(:)
    
    ! Deallocate local arrays to prevent memory leaks
    call free_vector(u_drift)
    call free_vector(u_diffusion)
    call free_vector(u_euler_maruyama)
    call free_vector(u_milstein_correction)
    
    return
  end subroutine milstein_step
  
  subroutine compute_milstein_correction(brownian_increment, milstein_correction)
    use ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
    implicit none
    real(real64), intent(in) :: brownian_increment(DIM)
    real(real64), allocatable, intent(out) :: milstein_correction(:)
    real(real64), allocatable :: square_brownian_increment(:)
    real(real64), parameter :: max_correction = 1.0e6_real64
    integer :: i
   
    call alloc_vector(milstein_correction, DIM)
    call alloc_vector(square_brownian_increment, DIM)
    
    square_brownian_increment = brownian_increment ** 2
    milstein_correction = sigma * b * (square_brownian_increment - delta)
    
    ! Check for NaN or infinite values and apply bounds
    do i = 1, DIM
      if (.not. ieee_is_finite(milstein_correction(i))) then
        print *, "Warning: Non-finite Milstein correction at index", i
        print *, "Value:", milstein_correction(i), "b(i):", b(i)
        print *, "brownian_inc(i):", brownian_increment(i)
        print *, "square_brownian_increment(i):", square_brownian_increment(i)
        print *, "(square - delta):", square_brownian_increment(i) - delta
        print *, "sigma:", sigma
        milstein_correction(i) = 0.0_real64
      else if (abs(milstein_correction(i)) > max_correction) then
        ! Limit extremely large corrections
        milstein_correction(i) = sign(max_correction, milstein_correction(i))
      end if
    end do
    
    ! Deallocate local array
    call free_vector(square_brownian_increment)
  end subroutine compute_milstein_correction
   
   
   
   !> @brief Solves a stochastic differential equation using the Milstein scheme.
  !>
  !> This subroutine implements the numerical solution of a stochastic differential
  !> equation using the Milstein scheme. It initializes the necessary variables,
  !> advances the solution step-by-step, and returns a status code indicating
  !> success or failure.
  !>
  !> @param[out] status Status code: 0 for success, non-zero for failure.
  !>
  !> The subroutine uses the `milstein_step` subroutine to advance the solution
  !> over time, starting from an initial condition and progressing to a final time.
  subroutine solve_sde_with_milstein(status)
    implicit none
    logical, intent(out) :: status
    
    ! Local variables
    real(real64), allocatable :: u_current(:), u_next(:), brownian_inc(:)
    character(len=100) :: file_name
    character(len=20) :: header(DIM)
    integer :: i, n_steps, j
    
    ! Initialize status to success
    status = .FALSE.
    
    ! Allocate arrays
    call alloc_vector(u_current, DIM)
    call alloc_vector(brownian_inc, DIM)
    
    ! Set initial condition using global u vector
    u_current(:) = u_zero(:)
    path(0, :) = u_zero(:)
    ! Time stepping loop
    do i = 1, nobs
      ! Check for NaN or extremely large values in current solution
      if (any(ieee_is_nan(u_current))) then
        print *, "NanN found at iteration:"
        print *, "i: ", i
        stop
      end if
      if (any(abs(u_current) > 1.0e10_real64)) then
        print *, "Solution overflow detected at iteration:", i
        print *, "Max value:", maxval(abs(u_current))
        stop
      end if
      
      call vectorial_winner_increment(u_current, brownian_inc)
      
      ! Check brownian increment for overflow
      if (any(abs(brownian_inc) > 1.0e10_real64)) then
        print *, "Brownian increment overflow at iteration:", i
        print *, "Max brownian:", maxval(abs(brownian_inc))
        stop
      end if
      
      call milstein_step(u_current, brownian_inc, u_next)
      
      ! Apply bounds to the solution
      do j = 1, DIM
        if (abs(u_next(j)) > 1.0e6_real64) then
          u_next(j) = sign(1.0e6_real64, u_next(j))
        end if
      end do
      
      path(i, :) = u_next(:)! Update current solution
      u_current(:) = u_next(:)
    end do
    file_name="../data/path_sample.bin"
    ! Create simple header
    call save_real64_2d_array_to_binary(file_name, path)
    status = .TRUE.
    call print_matrix_with_indices('head(path)', path(0:5, 0:5), 5, 5)
  end subroutine solve_sde_with_milstein
   !> Reshapes a 1D array into a 2D array using column-major order.
   !!
   !! This subroutine maps a linear vector of length `Nx * Ny` into a 2D array
   !! of shape `(Nx, Ny)`, using column-major indexing:
   !!
   !! \f[
   !! u0\_proj\_np1(i, j) = u0\_proj\_row\_np1(i + (j - 1) \cdot Nx)
   !! \f]
   !!
   !! This is equivalent to MATLAB's reshaping logic:
   !! \code{.m}
   !! for i = 1:Nx
   !!     for j = 1:Ny
   !!         m = i + (j - 1) * Nx;
   !!         u(i,j) = u_vec(m);
   !!     end
   !! end
   !! \endcode
   !!
   !! @param[in]  u0_proj_row_np1  A 1D array of length Nx * Ny (modal vector)
   !! @param[in]  Nx               Number of grid points in x-direction
   !! @param[in]  Ny               Number of grid points in y-direction
   !! @param[out] u0_proj_array      2D reshaped array of size (Nx, Ny)
   subroutine reshape_to_2d(u0_proj_row, u0_proj_array)
     use iso_fortran_env, only: real64
     implicit none
     real(real64), intent(in) :: u0_proj_row(DIM)
     real(real64), intent(out) :: u0_proj_array(Nx, Ny)
     integer :: i, j, m
     do j = 1, Ny
       do i = 1, Nx
         m = i + (j - 1) * Nx
         u0_proj_array(i, j) = u0_proj_row(m)
       end do
     end do
   end subroutine reshape_to_2d
   
   !> Projects modal coefficients onto a uniform 2D spatial grid using a cosine basis.
   !!
   !! This subroutine evaluates a modal expansion of the form:
   !! \f[
   !! u(x_i, y_j) = \sum_{m=0}^{Nx-1} \sum_{n=0}^{Ny-1}
   !! u0\_proj(m,n) \cdot \phi_m(x_i) \cdot \phi_n(y_j)
   !! \f]
   !! where the basis functions \f$ \phi_k(z) \f$ are:
   !! \f[
   !! \phi_k(z) = \sqrt{1 + \mathrm{sign}(k)} \cdot \cos\left(\frac{\pi k z}{L} \right)
   !! \f]
   !!
   !! @param[in]  u0_proj   Modal coefficient matrix of size (0:Nx-1, 0:Ny-1)
   !! @param[in]  Nx        Number of spatial grid points in x-direction
   !! @param[in]  Ny        Number of spatial grid points in y-direction
   !! @param[in]  Lx        Length of the domain in x-direction
   !! @param[in]  L2        Length of the domain in y-direction
   !! @param[out] u_grid    Reconstructed solution array on the grid (Nx, Ny)
   !!
   !! @note Grid points are located at cell centers:
   !!       \f$ x_i = dx \cdot (i - 1/2),\quad y_j = dy \cdot (j - 1/2) \f$
   subroutine project_modal_to_grid(u_proj, u_grid)
     use iso_fortran_env, only: real64
     implicit none
     
     real(real64), intent(in) :: u_proj(0:Nx - 1, 0:Ny - 1)
          ! Output
     real(real64), intent(out) :: u_grid(Nx, Ny)
     real(real64) :: dx, dy, xi, yj, hi, hj
     integer :: ix, iy, i, j
     
     ! Compute uniform grid spacing
     dx = L1 / real(Nx, real64)
     dy = L2 / real(Ny, real64)
     
     ! Initialize grid solution to zero
     u_grid = 0.0_real64
     ! Evaluate modal expansion at each grid point
     do iy = 1, Ny
       yj = dy * (real(iy, real64) - 0.5d0)
       do ix = 1, Nx
         xi = dx * (real(ix, real64) - 0.5d0)
         do i = 0, Nx - 1
           hi = sqrt(1.0d0 + sign(1.0d0, real(i, real64))) &
                   &* cos(PI * real(i, real64) * xi / L1)
           do j = 0, Ny - 1
             hj = sqrt(1.0d0 + sign(1.0d0, real(j, real64))) &
                     & * cos(PI * real(j, real64) * yj / L2)
             u_grid(ix, iy) = u_grid(ix, iy) + u_proj(i, j) * hi * hj
           end do
         end do
       end do
     end do
   end subroutine project_modal_to_grid
 end module mod_sde_solver
