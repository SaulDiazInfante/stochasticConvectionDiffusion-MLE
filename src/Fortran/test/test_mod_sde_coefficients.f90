program test_mod_sde_coefficients
    use iso_fortran_env, only: real64, int32
    use mod_sde_coefficients
    use mod_global_parameters_and_shared_data
    use mod_data_io
    use mod_setup_and_logs
    implicit none
    
    integer(int32) :: i, j
    real(real64) :: eps
    logical :: status
    
    real(real64), allocatable :: vector_drift(:)
    real(real64), allocatable :: vector_diffusion(:)
    real(real64), allocatable :: expected_mat(:,:)
    real(real64), allocatable :: expected_vec(:), vec_diff(:), diffusion_derivative(:)
    
    eps = 1.0e-12_real64
    
    status = .FALSE.
    
    ! Allocate memory first
    call allocate_dynamic_memory()
    call gen_observation_times()
    call gen_eigen_values()
    call gen_matrix_diag_B()
    call gen_matrix_B()
    call gen_lambda_matrix()
    call assemble_matrix_A()
    call gen_drift_matrix()
    call gen_diffusion_matrix()
    
   ! Allocate local test arrays
    call alloc_array(expected_mat, DIM, DIM)
    call alloc_vector(expected_vec, DIM)
    call alloc_vector(vector_diffusion, DIM)
    call alloc_vector(vec_diff, DIM)
    call alloc_vector(diffusion_derivative, DIM)
    
    u(:) = 1.0_real64
    expected_mat = driftmat
    call display_parameters()
    call display_domain_problem_arrays()
    
    call eval_whole_drift(u, vector_drift)
    call eval_drift_at_u(u, expected_vec)
    vec_diff = -1.0 * vector_drift + expected_vec
    
    if (all(abs(vec_diff) < eps)) then
        status = .TRUE.
        print *, "eval_whole_drift and eval_dirft_at_u TEST PASSED"
        print '(A10, ES30.16)', '|error|:', sum(abs(vec_diff))
        call print_vector_with_indices("whole_vector_drift output:", vector_drift(1:10), 10)
        call print_vector_with_indices("vector_drift_at u output:", expected_vec(1:10), 10)
        call print_vector_with_indices("error", vec_diff(1:10), 10)
    else
        print *, "eval_whole_drift and eval_dirft_at_u TEST FAILED"
        print '(A10, ES30.16)', '|error|:', sum(abs(vec_diff))
        call print_vector_with_indices("whole_vector_drift output:", vector_drift(1:10), 10)
        call print_vector_with_indices("vector_drift_at u output:", expected_vec(1:10), 10)
        call print_vector_with_indices("error", vec_diff(1:10), 10)
    end if
    
    expected_vec = 0.0_real64
    vec_diff = 0.0_real64
    call eval_diffusion_at_u(u, vector_diffusion)
    call eval_diagonal_diffusion_at_u(u, expected_vec)
    vec_diff = vector_diffusion - expected_vec
    
    if (all(abs(vec_diff) < eps)) then
        status = .TRUE.
        print *, "eval_diffusion_at_u and eval_diagonal_at_u TEST PASSED"
        print '(A10, ES30.16)', '|error|:', sum(abs(vec_diff))
        call print_vector_with_indices("whole_vector_drift output:", vector_diffusion(1:10), 10)
        call print_vector_with_indices("vector_drift_at u output:", expected_vec(1:10), 10)
        call print_vector_with_indices("error", vec_diff(1:10), 10)
    else
        print *, "eval_diffusion_at_u and eval_diagonal_at_u TEST FAILED"
        print '(A10, ES30.16)', '|error|:', sum(abs(vec_diff))
        call print_vector_with_indices("eval_diffusion_at_u output:", vector_diffusion(1:10), 10)
        call print_vector_with_indices("eval_diagonal_diffusion_at_u output:", expected_vec(1:10), 10)
        call print_vector_with_indices("error", vec_diff(1:10), 10)
    end if
    
    diffusion_derivative = 0.0_real64
    call compute_diffusion_derivative(diffusion_derivative)
    expected_vec = sigma * b(:)
    vec_diff = diffusion_derivative - expected_vec
    
    if (all(abs(vec_diff) < eps)) then
        status = .TRUE.
        print *, "eval_diffusion_derivative_at_u and eval_diagonal_at_u TEST PASSED"
        print '(A10, ES30.16)', '|error|:', sum(abs(vec_diff))
        call print_vector_with_indices("eval_diffusion_derivative output:", diffusion_derivative(1:10), 10)
        call print_vector_with_indices("expected (sigma * b) output:", expected_vec(1:10), 10)
        call print_vector_with_indices("error", vec_diff(1:10), 10)
    else
        print *, "eval_diffusion_derivative at_u and eval_diagonal_at_u TEST FAILED"
        print '(A10, ES30.16)', '|error|:', sum(abs(vec_diff))
        call print_vector_with_indices("eval_diffusion_derivative output:", diffusion_derivative(1:10), 10)
        call print_vector_with_indices("eval_diagonal_diffusion_derivative:", expected_vec(1:10), 10)
        call print_vector_with_indices("error", vec_diff(1:10), 10)
    end if
    
end program test_mod_sde_coefficients
