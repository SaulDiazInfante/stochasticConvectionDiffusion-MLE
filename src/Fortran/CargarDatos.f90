
!ifort CargarDatos.f90 -o datos -lm -C 

program	CargarDatos
use, intrinsic :: iso_c_binding
implicit none
	!Integer
	integer(C_INT) :: i,j,k,l,m,n
	!Allocatable
	!Real
	real(C_DOUBLE) ,dimension(:,:), allocatable :: A
!----------------PARAMETERS----------------------------------------------

	!Grid size:
   	Nx  = 10; Ny  = 10;
   	N2 = Nx*Ny;

!-------Allocation----------------------

	!Real
	allocate(A(N2,N2)); 
	
!---------------Load the variables----------
	
	open(unit=99,file='Arow.dat',form="formatted",status="old",action="read")
	do i = 1,Nx
		do j = 1,Ny
			m = i+(j-1)*Ny;
			do k = 1,Nx
				do l = 1,Ny
					n = k+(l-1)*Ny;
					read(99,*) A(m,n);
			end do
		end do
	end do
	close(99)

!------------------------------------------------------------------

	deallocate(A)

end program CargarDatos





