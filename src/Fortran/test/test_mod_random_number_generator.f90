program test_mkl_gaussian_sampler
    use mod_random_number_generator
    use mkl_vsl_type
    use mkl_vsl
    use iso_fortran_env, only : int32, real64
    implicit none
    
    integer, parameter :: n = 100000
    integer, parameter :: n_row = 2500
    integer, parameter :: n_col = 10000
    real(real64) :: mean_a, std_a, diff_mean, diff_std
    real(real64), allocatable :: samples(:)
    real(real64), allocatable :: gaussian_matrix_sample(:, :)
    real(real64) :: mean_val, std_val
    integer :: i
    real(real64), parameter :: tolerance = 5.0e-3
    integer :: user_seed
    logical :: test_passed
    logical :: debug
    
    ! Test parameters
    mean_a = 0.0_real64
    std_a = 1.0_real64
    user_seed = 123456
    debug = .false.
    
    call mkl_gaussian_sampler(n, mean_a, std_a, samples, user_seed)
    
    ! Compute mean
    mean_val = sum(samples) / n
    
    ! Compute standard deviation
    std_val = sqrt(sum((samples - mean_val)**2) / (n - 1))
    
    ! Check results
    test_passed = abs(mean_val - mean_a) < tolerance .and. abs(std_val - std_a) < tolerance
    
    if (test_passed) then
        print *, 'TEST PASSED'
    else
        print *, 'TEST FAILED'
        print *, '  Mean = ', mean_val
        print *, '  Std  = ', std_val
    end if
    
    
    call mkl_array_gaussian_sampler(&
        &n_row, &
        &n_col, &
        &mean_a, &
        &std_a, &
        &gaussian_matrix_sample, &
        &user_seed,&
        &debug&
    &)
    
    mean_val = sum(gaussian_matrix_sample) / real(n_row * n_col, kind=real64)
    std_val = sqrt(sum((gaussian_matrix_sample - mean_val)**2) / real(n_row * n_col - 1, kind=real64))
    
    ! Compare
    diff_mean = abs(mean_val - mean_a)
    diff_std  = abs(std_val - std_a)
    
    test_passed = (diff_mean < tolerance) .and. (diff_std < tolerance)
    
    if (test_passed) then
        print *, 'TEST PASSED'
    else
        print *, 'TEST FAILED'
        print *, '  Mean expected: ', mean_a, ' got: ', mean_val
        print *, '  Std  expected: ', std_a, ' got: ', std_val
        stop 1
    end if
end program test_mkl_gaussian_sampler
