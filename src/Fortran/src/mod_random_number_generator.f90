!> @ingroup modules
!> @author F. Baltazar-Larios, S. Diaz-Infante
!> @brief This module implements the Box-Muller algorithm to generate
!> random variables with standard Gaussian distribution from a
!> uniform distributed random variable. This module enclose tree functions.
!> @see Kloeden & Platen 1992
! include 'mkl_vsl.f90'
! include 'mkl_vsl.f90'
module mod_random_number_generator
    use mkl_vsl_type
    use mkl_vsl
    use iso_fortran_env, only : int32, real32
    implicit none
contains

!> @brief Returns an array of size buffer_size of gaussian random variables with mean
!> mean_a amd standard deviation std_a
    subroutine mkl_gaussian_sampler(&
        &buffer_size, &
        &mean_a, &
        &std_a, &
        &gaussian_sample, &
        &user_seed, &
        &debug&
    &)
        integer, intent(in) :: buffer_size
        real(kind=8), intent(in) :: mean_a, std_a 
        real(kind=8), intent(out), allocatable :: gaussian_sample(:)
        integer, intent(in), optional :: user_seed
        logical, intent(in), optional :: debug
        TYPE (VSL_STREAM_STATE) :: stream
        
        integer(kind=4) errcode
        integer(kind=4) i, j
        real (kind=8) :: rand
        integer brng, seed, method, n
        integer, allocatable :: new (:), old(:)
        integer, parameter :: lower = 1, upper = 100
        if (present(user_seed)) then
            seed = user_seed
        else
            call random_seed()  ! Initialize the random number generator
            call random_number(rand)
            seed = lower + int(rand * (upper - lower + 1))
        end if
        gaussian_sample = [(0.0, i=1, buffer_size)]
    !       ***** Initializing *****
        brng = VSL_BRNG_MT19937 !! Mersenne Twister
        method = VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2 !! For Gaussian distribution
        errcode = vslNewStream(stream, brng, seed)
        errcode = vdRngGaussian(&
            method, stream, buffer_size, gaussian_sample, mean_a, std_a&
        )
        if (present(debug)) then
            if (debug) then
                print *,"Error structure = ", stream
                print *, "Gaussian stream: ", gaussian_sample(1:5)
            end if
        end if
    !       ***** Deinitialize *****
        errcode = vslDeleteStream(stream)
    end subroutine mkl_gaussian_sampler

    subroutine mkl_array_gaussian_sampler(&
        &n_row, &
        &n_col, &
        &mean_a, &
        &std_a, &
        &gaussian_sample, &
        &user_seed, &
        &debug&
    &)
        integer(int32), intent(in) :: n_row, n_col
        real(kind=8), intent(in) :: mean_a, std_a 
        real(kind=8), intent(out), dimension(:,:), allocatable :: gaussian_sample
        integer, intent(in), optional :: user_seed
        logical, intent(in), optional :: debug
        TYPE (VSL_STREAM_STATE) :: stream
        
        integer(kind=4) errcode
        integer(kind=4) i, j
        real (kind=8) :: rand
        integer brng, seed, method, n
        integer, allocatable :: new (:), old(:)
        integer, parameter :: lower = 1, upper = 100
        if (present(user_seed)) then
            seed = user_seed
        else
            call random_seed()  ! Initialize the random number generator
            call random_number(rand)
            seed = lower + int(rand * (upper - lower + 1))
        end if
        gaussian_sample = reshape([(0.0, i=1, n_row*n_col)], [n_row, n_col])

    !       ***** Initializing *****
        brng = VSL_BRNG_MT19937 !! Mersenne Twister
        method = VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2 !! For Gaussian distribution
        errcode = vslNewStream(stream, brng, seed)
        errcode = vdRngGaussian(&
            &method, &
            & stream, & 
            & n_row*n_col, &
            & gaussian_sample, & 
            &mean_a, &
            & std_a&
        )
        if (present(debug)) then
            if (debug) then
                print *,"Error structure = ", stream
                print *, "Gaussian stream: ", gaussian_sample(1:5, 1:5)
            end if
        end if
    !       ***** Deinitialize *****
        errcode = vslDeleteStream(stream)
    end subroutine mkl_array_gaussian_sampler
end module mod_random_number_generator

! ifort -c -i8 -I"${MKLROOT}/include" mod_random_number_generator.f90
