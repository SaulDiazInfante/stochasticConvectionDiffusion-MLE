module mod_sde_solver
  use iso_fortran_env, only: int32, real64
  implicit none
  contains
  pure subroutine MilsteinStep(DIM, delta, startx, alpha, sigma, endx, brown)
    implicit none
    !! Accordingly with the Milstein scheme
    !! X_{n+1} ^ {\Delta} = f(X_n) \Delta + g(X_n) \Delta W_{n+1}  + (Milstein corr)
    !! see Kloeden & Platten (1994)
    integer(int32), intent(in) :: DIM
    real(real64), intent(in) :: delta
    !! Step-size
    real(real64), intent(in) :: startx(DIM)
    !! $X_n$ in the Milstein recurrence
    real(real64), intent(inout) :: endx(DIM)
    !! $X_{n+1}$ in the above expression
    !TODO; Review the following parameters
    real(real64), intent(in) :: alpha(DIM)
    real(real64), intent(in) :: sigma(DIM)
    real(real64), intent(inout) :: brown(DIM)

    real(real64) W(DIM), xx(DIM), sum_aux

    integer(int32) i,j, n_omega
    n_omega = 100

    do i=1,DIM
      xx(i)=0.0
    enddo
    call BrownianStep(DIM, delta, n_omega, xx, W)
    do i=1,DIM
      sum_aux = 0.0
      sum_aux = sum_aux + sigma(i) * W(i)
      endx(i) = startx(i) + alpha(i) * delta + sum_aux
      brown(i)=W(i)
    enddo
   return
  end subroutine MilsteinStep

  pure subroutine BrownianStep(DIM, delta, n_omega, startx, endx)
    implicit none
    integer(int32), intent(in) :: DIM
    !! Dimension
    integer(int32), intent(in) :: n_omega
    !! Number of Gaussian random variables to generate the winner increment.
    real(real64), intent(in) :: delta
    !! Step-size
    real(real64), intent(in) :: startx(DIM)
    !! old value in the recursive method
    real(real64), intent(inout) :: endx(DIM)

    real(real64) xi
    integer(int32) i
    !! TODO: Call mkl rng to create a vectorized version of the below code.
    do i=1,DIM
      !v call normalvar(var)
      endx(i) = startx(i) + sqrt(delta) !*var
    enddo
    return
  end subroutine BrownianStep
end module mod_sde_solver
