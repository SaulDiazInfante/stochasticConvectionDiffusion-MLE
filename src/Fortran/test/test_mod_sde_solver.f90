program test_mod_sde_solver
    use iso_fortran_env, only: int32, real64
    use ieee_arithmetic, only: ieee_is_nan
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
    logical :: status, has_nan
    character(len=20), dimension(2) :: header
    character(len=50) :: file_name
    real(real64), dimension(DIM, 2) :: indexed_times
    ! Removed unused variable declarations
    real(real64), allocatable :: vectorial_winner_delta(:), initial_vector_winner(:)
    real(real64), allocatable :: current_brownian_point_path(:), next_brownian_point_path(:)
    real(real64), allocatable :: current_u(:), next_u(:), u_proj(:), u_grid(:,:)
    real(real64), allocatable :: u_proj_2d(:, :)
    real(real64), parameter :: eps = 1.0e-12_real64
    integer(int32), allocatable :: row_nan(:), col_nan(:)
    integer(int32) :: rows, cols, k, count_nan
    
    call build_sde()
    call display_parameters()
    call display_domain_problem_arrays()
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
    ! next_u is allocated by milstein_step, so don't allocate here
    call milstein_step(&
        &current_u, &
        &vectorial_winner_delta, &
        &next_u &
    &)
    if (all(abs(next_u) < eps)) then
        print *, "ERROR: Milstein step is stuck TEST FAILED"
        call print_vector_with_indices("U_{n+1} ", next_u(1:5), 5)
    else
        print *, "Non zero milstein_step TEST PASSED"
        call print_vector_with_indices("U_{n+1} ", next_u(1:5), 5)
    end if
    has_nan = any(ieee_is_nan(next_u))
    if (has_nan) then
        print *, "milsten_step returns an array with at least a NaN. TEST FAILED"
        status = .FALSE.
    else
        print *, "milsten_step returns an array without NaN values. TEST PASSED"
        status = .TRUE.
    end if
    
    status = .FALSE.
    has_nan = .TRUE.
    call solve_sde_with_milstein(status)
    has_nan = any(ieee_is_nan(path))
    if (has_nan) then
        print *, "At least one observation of the sampled path has a NaN. TEST FAILED"
        call find_nan_indices_2d(path, 1001, 2500, row_nan, col_nan, count_nan)
        print *, "NaN found at:"
        do k = 1, count_nan
            print *, "  (", row_nan(k), ",", col_nan(k), ")"
        end do
        status = .FALSE.
        
    else
        print *, "The sampled path does not contain NaN. TEST PASSED"
        status = .TRUE.
    end if
    ! call save_real64_2d_array_to_binary("../data/sampled_path.bin", path)
    
    call alloc_vector(u_proj, DIM)
    call alloc_array(u_grid, Nx, Ny)
    call alloc_array(u_proj_2d, Nx, Ny)
    u_proj = path(nobs, :)
    print *, "--------------------------------------------------"
    call print_vector_with_indices("u_proj:", u_proj(1:5), 5)
    call reshape_to_2d(u_proj, u_proj_2d)
    call print_matrix_with_indices("reshape(u_proj)", u_proj_2d(1:5, 1:5) ,5, 5)
    
    call project_modal_to_grid(u_proj_2d, u_grid)
    call print_matrix_with_indices("U_{grid}: ", u_grid(1:5, 1:5), 5, 5)
    ! Deallocate local arrays
    call free_vector(initial_vector_winner)
    call free_vector(vectorial_winner_delta)
    call free_vector(current_u)
    call free_vector(next_u)
    if (allocated(next_brownian_point_path)) deallocate(next_brownian_point_path)
    if (allocated(current_brownian_point_path)) deallocate(current_brownian_point_path)
    if (allocated(row_nan)) deallocate(row_nan)
    if (allocated(col_nan)) deallocate(col_nan)
    
    ! call deallocate_all_shared_data()  ! Commented out to avoid segfault with TBB allocator
end program test_mod_sde_solver
