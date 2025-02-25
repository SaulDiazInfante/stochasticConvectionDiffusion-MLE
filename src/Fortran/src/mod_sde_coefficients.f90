module mod_sde_coefficients
  use iso_fortran_env, only: int32, real64
  implicit none
contains
  subroutine MDrift(DIM, theta, beta,lambda_mat, A, drift_mat)
    implicit none
    integer(int32), intent(in) ::  DIM
    !! dimension $(N_x * N_y)$
    real(real64), intent(in) :: theta
    !! meaning
    real(real64), intent(in) :: beta
    !! meaning
    real(real64), intent(in) :: A(DIM,DIM)
    !! meaning
    real(real64), intent(in) :: lambda_mat(DIM,DIM)
    !!
    real(real64), intent(out) ::lambda_mat(DIM,DIM)
    !!
    integer(int32) i, j
  end subrotine Mdrift

  subroutine MDiff(DIM,sigma,B,Tsigma)
    implicit none
        integer DIM
    real sigma,B(DIM,DIM),Tsigma(DIM,DIM)



         Tsigma(:,:)=sigma*B(:,:)


    return
    end

end mod_sde_coefficients
