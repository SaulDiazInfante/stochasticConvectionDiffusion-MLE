program check_snrm
    use, intrinsic :: iso_fortran_env, only : real32
    implicit none
    real(real32), external :: snrm2
    print *, snrm2(1, [0._real32], 1)
end program