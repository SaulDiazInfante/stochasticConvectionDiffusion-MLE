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
  use mod_random_number_generator
  use mod_par_generators
  use mod_sde_coefficients
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
      & + u_drift(:) * delta &
      & + u_diffusion(:) * brownian_increment(:)
    call compute_milstein_correction(brownian_increment, u_milstein_correction)
    next_u(:) = u_euler_maruyama(:) + u_milstein_correction(:)
    return
  end subroutine milstein_step
  
  subroutine compute_milstein_correction(brownian_increment, milstein_correction)
    implicit none
    real(real64), intent(in) :: brownian_increment(DIM)
    real(real64), allocatable, intent(out) :: milstein_correction(:)
    real(real64), allocatable :: square_browinian_increment(:)
   
    call alloc_vector(milstein_correction, DIM)
    call alloc_vector(square_browinian_increment, DIM)
    
    square_browinian_increment = brownian_increment ** 2
    milstein_correction = sigma * b * (square_browinian_increment - delta)
    deallocate(square_browinian_increment)
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
  subroutine solve_sde(status)
    implicit none
    integer, intent(out) :: status
    
    ! Local variables
    real(real64), allocatable :: u_current(:), u_next(:), brownian_inc(:)
    integer :: i, n_steps
    
    ! Initialize status to success
    status = 0
    
    ! Initialize parameters
    n_steps = 100 ! Number of time steps, adjust as needed
    
    ! Allocate arrays
    call alloc_vector(u_current, DIM)
    call alloc_vector(brownian_inc, DIM)
    
    ! Set initial condition using global u vector
    u_current(:) = u(:)
    
    ! Time stepping loop
    do i = 1, n_steps
      ! Generate Brownian increment
      call vectorial_winner_increment(u_current, brownian_inc)
      
      ! Advance solution using Milstein scheme
      call milstein_step(u_current, brownian_inc, u_next)
      
      ! Update current solution
      u_current(:) = u_next(:)
      
      ! Deallocate u_next as it's reallocated in milstein_step
      if (allocated(u_next)) deallocate(u_next)
    end do
    
    ! Clean up
    ! Store the final state in the global u vector
    u(:) = u_current(:)
    
    ! Clean up
    if (allocated(u_current)) deallocate(u_current)
    if (allocated(brownian_inc)) deallocate(brownian_inc)
    if (allocated(u_next)) deallocate(u_next)
    
    return
  end subroutine solve_sde
end module mod_sde_solver
