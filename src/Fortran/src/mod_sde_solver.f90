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
  ! use MKL_VSL
  use mod_random_number_generator

  implicit none
  contains

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
  !> @param[in] n_omega The number of random variables (not used in current implementation).
  !> @param[in] startx The starting position vector of dimension `DIM`.
  !> @param[inout] endx The updated position vector of dimension `DIM`.
  !---------------------------------------------------------------------------
  pure subroutine BrownianStep(DIM, delta, n_omega, startx, endx)
    implicit none
    integer(int32), intent(in) :: DIM
    integer(int32), intent(in) :: n_omega
    real(real64), intent(in) :: delta
    real(real64), intent(in) :: startx(DIM)
    real(real64), intent(inout) :: endx(DIM)
    !!
    real(real64) xi, dd, ddW(n_omega), dW
    integer(int32) i
    
    !! TODO: Call mkl rng to create a vectorized version of the below code.
    do i =1, DIM
      !v call normalvar(var)
      endx(i) = startx(i) + sqrt(delta) !*var
    end do
    return
  end subroutine BrownianStep

  !> @brief Computes the increment of a Wiener process over a given time interval.
  !!
  !! This subroutine calculates the increment of a Wiener process (also known as Brownian motion)
  !! over a specified time interval `delta`. The increment is computed using a Gaussian sampler.
  !!
  !! @param[in] delta The time interval over which the Wiener process increment is computed.
  !! @param[in] n_omega The number of subintervals to divide `delta` into for the Gaussian sampling.
  !! @param[in] user_seed (optional) An optional seed for the random number generator.
  !! @param[in] winner_0 The initial value of the Wiener process.
  !! @param[out] winner_delta The computed increment of the Wiener process over the interval `delta`.
  !!
  !! The subroutine uses the MKL Gaussian sampler to generate random samples from a normal distribution
  !! with mean 0 and standard deviation 1. These samples are then scaled by the square root of the
  !! subinterval length and summed to produce the Wiener process increment.
  subroutine winner_increment(delta, n_omega, winner_0, winner_delta, user_seed)
    implicit none
    real(real64), intent(in) :: delta
    integer(int32), intent(in) :: n_omega
    integer, intent(in), optional :: user_seed
    real(real64), intent(in) :: winner_0
    real(real64), intent(out) :: winner_delta 
    
    real(real64) dd
    real(kind=8), allocatable :: ddW(:)

    dd = delta / n_omega
    if (present(user_seed)) then
      call mkl_gaussian_sampler(n_omega, 0.0_real64, 1.0_real64, ddW, user_seed)
    else
      call mkl_gaussian_sampler(n_omega, 0.0_real64, 1.0_real64, ddW)
    end if
    ddW = sqrt(dd) * ddW
    winner_delta = winner_0 + sum(ddW(:))
    return
  end subroutine winner_increment
  
  
  !> @brief Computes the increment of a vectorial Wiener process.
  !>
  !> This subroutine calculates the increment of a vectorial Wiener process
  !> over a given time step `delta`. The Wiener process is represented by
  !> `winner_0` and the increment is stored in `winner_delta`.
  !>
  !> @param delta The time step over which the increment is computed.
  !> @param DIM The dimension of the Wiener process.
  !> @param n_omega The number of sub-steps within the time step `delta`.
  !> @param user_seed (Optional) Seed for the random number generator.
  !> @param winner_0 The initial value of the Wiener process.
  !> @param winner_delta The computed increment of the Wiener process.
  !>
  !> The subroutine uses the MKL library to sample from a Gaussian distribution
  !> and scales the samples appropriately to compute the Wiener process increment.
  subroutine vectorial_winner_increment(&
      &delta, &
      &DIM, &
      &n_omega, & 
      &winner_0, & 
      &winner_delta, & 
      &user_seed)
    implicit none
    real(real64), intent(in) :: delta
    integer(int32), intent(in) :: DIM
    integer(int32), intent(in) :: n_omega
    integer, intent(in), optional :: user_seed
    real(real64), intent(in) :: winner_0(DIM) 
    real(real64), intent(out) :: winner_delta(DIM) 
    
    real(real64) dd
    real(real64), allocatable, dimension(:, :) :: ddW
    integer(int32) i

    dd = delta / n_omega
    if (present(user_seed)) then
      call mkl_array_gaussian_sampler(DIM, n_omega, 0.0_real64, 1.0_real64, ddW, user_seed)
    else
      call mkl_array_gaussian_sampler(DIM, n_omega, 0.0_real64, 1.0_real64, ddW)
    end if
    ddW = sqrt(dd) * ddW
    do i=1, DIM
      winner_delta(i) = winner_0(i) + sum(ddW(i, :))
    end do
    return
  end subroutine vectorial_winner_increment

  
  
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
  pure subroutine MilsteinStep(DIM, delta, startx, alpha, sigma, endx, brown)
    implicit none
    integer(int32), intent(in) :: DIM
    real(real64), intent(in) :: delta
    real(real64), intent(in) :: startx(DIM)
    real(real64), intent(inout) :: endx(DIM)
    real(real64), intent(in) :: alpha(DIM)
    real(real64), intent(in) :: sigma(DIM)
    real(real64), intent(inout) :: brown(DIM)
    real(real64) W(DIM), xx(DIM), sum_aux

    integer(int32) i, j, n_omega
    n_omega = 100
    
    do i=1,DIM
      xx(i)=0.0
    end do

    call BrownianStep(DIM, delta, n_omega, xx, W)
    do i=1,DIM
      sum_aux = 0.0
      sum_aux = sum_aux + sigma(i) * W(i)
      endx(i) = startx(i) + alpha(i) * delta + sum_aux
      brown(i)=W(i)
    end do
   return
  end subroutine MilsteinStep

  end module mod_sde_solver
