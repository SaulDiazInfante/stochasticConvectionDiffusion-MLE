program test_mod_sde_solver
    use iso_fortran_env, only: int32, real64
    use mod_alloc
    use mod_global_parameters_and_shared_data
    use mod_sde_coefficients
    use mod_sde_solver
    implicit none
    
    ! Variables for test
    integer :: status, i
    
    ! Initialize status
    status = 0
    
    ! Initialize required memory and data structures
    call allocate_dynamic_memory()
    
    ! Set initial values for testing
    do i = 1, DIM
        lambdas(i) = real(i, real64)
        u(i) = 1.0_real64
        b(i) = 0.1_real64
    end do
    
    ! Initialize matrices A and B with test values
    do i = 1, DIM
        a(i, i) = 1.0_real64
        bmat(i, i) = 0.5_real64
    end do
    
    ! Generate drift and diffusion matrices
    call gen_drift_matrix()
    call gen_diffusion_matrix()
    
    ! Call the SDE solver with the correct name
    call solve_sde(status)
    
    if (status == 0) then
        print *, "PASS: solve_sde executed without error"
    else
        print *, "FAIL: solve_sde returned nonzero status"
        stop 1
    end if
    
    ! Clean up (deallocate memory if needed)
    ! Note: In a real application, you might need more cleanup
    
end program test_mod_sde_solver
