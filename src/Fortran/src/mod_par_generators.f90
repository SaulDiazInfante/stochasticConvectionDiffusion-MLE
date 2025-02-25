module mod_par_generators
  use iso_fortran_env, only: int32, real64
  implicit none
contains
  pure subroutine gen_observation_times(nobs, delta, times)
    !! Returns the array times which is  a stencil of len n_obs
    !! with step-size delta
    implicit none
    ! Arguments
    integer(int32), intent(in) :: nobs
      !! The number of observation in the whole time interval
    integer(int32) :: i
    real(real64), intent(in) :: delta
      !! Step-size, such that nobs * delta = \(T\)
    real(real64), intent(out) :: times(0:nobs)
      !! Time stencil.
    do i=0, nobs
      times(i) = i * delta
    end do
    return
  end subroutine gen_observation_times

  pure subroutine gen_lambdas(DIM, L1, L2, lambda_numbers)
    implicit none
    integer(int32), intent(in) :: DIM
      !! Spatial dimension of the vector space
      !! \( \mathcal{O}\subset \mathbb{R}^{d}\).
    real(real64), intent(in) :: L1
    real(real64), intent(in) :: L2
          !!
    real(real64), intent(out):: lambda_numbers(DIM)
      !! Entries of matrix \( \Lambda \) from equation
      !!  $$
      !!    dU =
      !!      \big( -\beta \Lambda\,   -\theta A\,\big)  U   dt
      !!      + \sigma B\, U\, dW(t)
      !!  $$
    real(real64), parameter :: PI = 2.d0 * DACOS (1.d0)
    real(real64)  pi_square, L1_res, L2_res, pi_square_sum_l1_res_l2_res
    integer(int32) i
!
    pi_square = PI ** 2
    L1_res = (L1) ** (-1)
    L2_res = (L2) ** (-1)
    pi_square_sum_l1_res_l2_res = (PI ** 2) * (L1_res + L2_res)
    do i=1, DIM
      lambda_numbers(i)= pi_square_sum_l1_res_l2_res * (i**2)
    enddo
     return
  end subroutine gen_lambdas

  pure subroutine MB(DIM, lambdas, gamma, B)
    implicit none
    integer(int32), intent(in) :: DIM
      !! Spatial dimension of the vector space
      !! \( \mathcal{O}\subset \mathbb{R}^{d}\).
    real(real64), intent(in) :: gamma
    real(real64), intent(in) :: lambdas(DIM)
    real(real64), intent(out) :: B(DIM, DIM)

    integer(int32) i
!
    B(:,:)=0.0
    do i=1, DIM
      B(i,i) = lambdas(i) ** (-gamma)
    enddo
    B(1,1)=1.0
    return
  end subroutine MB
  !TODO: Rename this matrix
  pure subroutine gen_lambda_matrix(DIM, lambdas, lambda_matrix)
    implicit none
    integer(int32), intent(in) :: DIM
    real(real64), intent(in) :: lambdas(DIM)
    real(real64), intent(out) :: lambda_matrix(DIM, DIM)

    integer(int32) i
    lambda_matrix(:,:)=0.0

    do i=1,DIM
      lambda_matrix(i,i)=lambdas(i)
    enddo
    return
  end subroutine gen_lambda_matrix
  pure subroutine MA(DIM, hs, A)
    implicit none
    integer(int32), intent(in) :: DIM
    real(real64), intent(in) :: hs(DIM)
    real(real64), intent(out) :: A(DIM,DIM)
    integer(int32) :: i
    do i=1,DIM
      A(i,:)=hs(:)
    enddo
    return
  end subroutine MA
end module mod_par_generators
