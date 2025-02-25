module mod_generators
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

  pure subroutine gen_lambdas(DIM, aum, lambdas)
    implicit none
    integer(int32), intent(in) :: DIM
      !! Spatial dimension of the vector space 
      !! \( \mathcal{O}\subset \mathbb{R}^{d}\).
    real(real64), intent(in) :: aum
      !! 
    real(real64), intent(out):: lambdas(DIM)
      !! Entries of matrix \( \Lambda \) from equation
      !!  $$
      !!    dU =  
      !!      \big( -\beta \Lambda\,   -\theta A\,\big)  U   dt 
      !!      + \sigma B\, U\, dW(t)
      !!  $$
    integer(int32) i
!
    lambdas(1)=0.0
    do i=2,DIM
      lambdas(i)=aum+lambdas(i-1)
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

  pure subroutine MG(DIM, lambdas, G)
    implicit none
    integer(int32), intent(in) :: DIM
    real(real64), intent(in) :: lambdas(DIM)
    real(real64), intent(out) :: G(DIM, DIM)

    integer(int32) i
    G(:,:)=0.0

    do i=1,DIM
      G(i,i)=lambdas(i)
    enddo
    return
  end subroutine MG

  pure subroutine gen_hs(DIM, x, Ls, hs)
    implicit none
    integer(int32), intent(in) :: DIM

    real(real64), intent(in) :: Ls(DIM)
    real(real64), intent(in) :: x
    real(real64), intent(out) :: hs(DIM)

    integer(int32) i,j
    real(real64), parameter :: PI = 3.141592
    do i=1,DIM
      hs(i)=COS(PI * i * x / Ls(i))
    enddo

    return
  end subroutine gen_hs

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

  pure subroutine MDrift(DIM, theta, beta, G, A, Talpha)
    implicit none

    integer(int32), intent(in) :: DIM
    real(real64), intent(in) :: theta
    real(real64), intent(in) :: beta
    real(real64), intent(in) :: G(DIM, DIM)
    real(real64), intent(in) :: A(DIM, DIM)
    real(real64), intent(out) :: Talpha(DIM, DIM)

    integer(int32) i,j

    Talpha(:,:) = -beta * G(:,:) - theta * A(:,:)
    return
  end subroutine MDrift


end module mod_generators
