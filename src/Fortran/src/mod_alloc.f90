module mod_alloc
    use iso_fortran_env, only: int32, real64
    implicit none
  
    private
    public :: alloc_vector, free_vector, alloc_array
  
  contains
  
    subroutine alloc_vector(a, n)
      real(real64), allocatable, intent(in out) :: a(:)
      integer(int32), intent(in) :: n
      integer(int32) :: stat
      character(len=100) :: errmsg
      if (allocated(a)) call free_vector(a)
      allocate(a(n), stat=stat, errmsg=errmsg)
      if (stat > 0) error stop errmsg
      a = 0.0_real64
    end subroutine alloc_vector
  
    subroutine free_vector(a)
      real(real64), allocatable, intent(in out) :: a(:)
      integer(int32) :: stat
      character(len=100) :: errmsg
      if (.not. allocated(a)) return
      deallocate(a, stat=stat, errmsg=errmsg)
      if (stat > 0) error stop errmsg
    end subroutine free_vector

    subroutine alloc_array(B, rows, cols)
        real(real64), allocatable, intent(inout) :: B(:,:)
        integer(int32), intent(in) :: rows, cols

        allocate(B(rows, cols))
        B = 0.0_real64
    end subroutine alloc_array

    subroutine free_array(B)
        real(real64), allocatable, intent(inout) :: B(:,:)
        if (allocated(B)) then
            deallocate(B)
            print *, "Array deallocated in subroutine."
        end if
    end subroutine free_array
  end module mod_alloc