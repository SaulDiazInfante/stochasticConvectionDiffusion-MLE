program	CargarDatos
	use iso_fortran_env, only: int32, real64
	implicit none
	integer(int32) :: i,j,k,l,m,n, Nx, Ny, N2
	real(real64) , dimension(:, :), allocatable :: A

	!----------------PARAMETERS----------------------------------------------
	!Grid size:
   	i=0
	j=0
	k=0
	l=0
	m=0
	n=0
	Nx  = 10 
	Ny  = 10
   	N2 = Nx * Ny

	!-------Allocation----------------------

	allocate(A(N2,N2)); 
!---------------Load the variables----------
	open(unit=99, &
		file='../../data/MatrixA.dat', & 
		form="formatted", & 
		status="old", &
		action="read" &
	)
	print*, Nx
	
!#if 0
	do i = 1, Nx
		do j = 1, Ny
			m = i + (j-1) * Ny;
 			do k = 1,Nx
				do l = 1,Ny
					n = k+(l-1)*Ny;
					read(99,*) A(m,n);
				enddo
			enddo
		enddo
	enddo
!#endif
close(99)

!------------------------------------------------------------------
	print*, A
	deallocate(A)
end program CargarDatos

!ifort -fpp CargarDatos.f90 -o datos -lm -C 



