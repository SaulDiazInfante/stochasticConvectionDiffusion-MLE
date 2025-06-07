program test_mod_sde_solver
    use iso_fortran_env, only: int32, real64
    use mod_data_io
    use mod_alloc
    use mod_global_parameters_and_shared_data
    use mod_par_generators
    use mod_sde_coefficients
    use mod_random_number_generator
    use mod_setup_and_logs
    use mod_sde_solver
    implicit none
    integer ::  i
    logical :: status
    character(len=20), dimension(2) :: header
    character(len=50) :: file_name
    real(real64), dimension(DIM, 2) :: indexed_times
    ! Removed unused variable declarations
    real(real64), allocatable :: vectorial_winner_delta(:), initial_vector_winner(:)
    real(real64), allocatable :: current_brownian_point_path(:), next_brownian_point_path(:)
    real(real64), allocatable :: current_u(:), next_u(:)
    real(real64), parameter :: eps = 1.0e-12_real64
    
    call build_sde()
    call display_parameters()
    call scalar_winner_increment(&
            &0.0_real64,&
            &winner_delta,&
            &SEED&
            &)
    status = .FALSE.
    if ( winner_delta /= 0.0_real64 ) then
        status = .TRUE.
        print *, 'Scalar Winner Increment TEST PASSED'
        print '(A10, ES12.5)',  'DeltaW: ', winner_delta
    else
        print *, 'Scalar Winner Increment  TEST FAILED'
    endif
    
    status = .FALSE.
    call alloc_vector(initial_vector_winner, DIM)
    call alloc_vector(vectorial_winner_delta, DIM)
    initial_vector_winner = 0.0_real64
    call vectorial_winner_increment(&
            &initial_vector_winner, &
            &vectorial_winner_delta, &
            &SEED)
    if (all(abs(vectorial_winner_delta) < eps)) then
        print *, "ERROR: all entries are zero"
        print *, "Vectporial Winner Increment  TEST FAILED"
        call print_vector_with_indices("Delta_W", vectorial_winner_delta(1:5), 5)
    else
        status = .TRUE.
        print *, 'Vecorial Winner Increment TEST PASSED'
        call print_vector_with_indices("Delta_W", vectorial_winner_delta(1:5), 5)
    end if
    
    allocate(next_brownian_point_path(DIM))
    allocate(current_brownian_point_path(DIM))
    
    next_brownian_point_path = 0.0_real64
    current_brownian_point_path = 0.0_real64
    call BrownianStep(current_brownian_point_path, next_brownian_point_path)
    status = .FALSE.
    if (all(abs(next_brownian_point_path) < eps)) then
        print *, "ERROR: Vectorial Browninan path stuck TEST FAILED"
        call print_vector_with_indices("W_{t} + \Dela W: ", next_brownian_point_path(1:5), 5)
    else
        print *, "BrownianStep TEST PASSED"
        call print_vector_with_indices("W_{t} + \Dela W: ", next_brownian_point_path(1:5), 5)
    end if
    
    call alloc_vector(current_u, DIM)
    call alloc_vector(next_u, DIM)
    call milstein_step(&
        &current_u, &
        &vectorial_winner_delta, &
        &next_u &
    &)
    if (all(abs(next_u) < eps)) then
        print *, "ERROR: Milstein step is stuck TEST FAILED"
        call print_vector_with_indices("U_{n+1} ", next_u(1:5), 5)
    else
        print *, "BrownianStep TEST PASSED"
        call print_vector_with_indices("U_{n+1} ", next_u(1:5), 5)
    end if
    
    deallocate(vectorial_winner_delta, initial_vector_winner)
    deallocate(next_brownian_point_path, current_brownian_point_path, current_u, next_u)
end program test_mod_sde_solver
