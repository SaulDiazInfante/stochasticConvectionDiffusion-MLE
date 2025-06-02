program test_mod_alloc
    use mod_alloc
    implicit none
    
    real(8), allocatable :: vec(:)
    integer :: n
    
    n = 5
    
    call alloc_vector(vec, n)
    
    if (allocated(vec) .and. size(vec) == n) then
        print *, "PASS: alloc_vector"
    else
        print *, "FAIL: alloc_vector"
        stop 1
    end if
end program test_mod_alloc
