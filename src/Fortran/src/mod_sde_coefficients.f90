module mod_sde_coefficients
  use iso_fortran_env, only: int32, real64
  implicit none
contains
  subroutine gen_drift_matrix(DIM, theta, beta, lambda_mat, A, drift_mat)
    implicit none
    integer(int32), intent(in) ::  DIM
    !! dimension $(N_x * N_y)$
    real(real64), intent(in) :: theta
    !! meaning
    real(real64), intent(in) :: beta
    !! meaning
    real(real64), intent(in) :: A(DIM, DIM)
    !! meaning
    real(real64), intent(in) :: lambda_mat(DIM, DIM)
    !!
    real(real64), intent(out) :: drift_mat(DIM,DIM)
    
    drift_mat = beta * lambda_mat(:,:) + theta * A(:, :)
    return
  end subroutine gen_drift_matrix

  subroutine gen_diffusion_matrix(DIM, sigma, B, diffusion_mat)
    implicit none
    integer(int32), intent(in) :: DIM
    !! >  dimension
    real(real64), intent(in) :: sigma
    !! Noise intensity
    real(real64), intent(in) :: B(DIM, DIM)

    real(real64), intent(out) :: diffusion_mat(DIM, DIM)

    diffusion_mat(:, :) = sigma * B(:, :)
    return
  end subroutine gen_diffusion_matrix
end module mod_sde_coefficients
