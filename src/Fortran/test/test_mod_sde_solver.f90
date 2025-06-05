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
    
    call gen_drift_matrix()
    call gen_diffusion_matrix()
    call solve_sde(status)
    
    if (status == 0) then
        print *, "PASS: solve_sde executed without error"
    else
        print *, "FAIL: solve_sde returned nonzero status"
        stop 1
    end if
end program test_mod_sde_solver
