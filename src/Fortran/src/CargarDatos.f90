program	CargarDatos
!	use iso_fortran_env, only: int32, real64
	implicit none
	integer :: i,j,k,l,m,n
	integer, parameter :: Nx=10, Ny=10
	!----------------PARAMETERS----------------------------------------------
	!Grid size:
	integer, parameter :: DIM = Nx *Ny
	real A(DIM, DIM), L1, L2
	real lambda_numbers(DIM)


	!-------Allocation----------------------

!	allocate(A(DIM, DIM)); 
!---------------Load the variables----------

!	open(unit=99, &
!		file='../../data/MatrixA.dat', & 
!		form="formatted", & 
!		status="old", &
!		action="read" &
!	)
!print*, Nx
	
	call MA(DIM,Nx,Ny, A)
	print*, "matrix A", A(1,1:2) 

	call gen_lambdas(DIM,lambda_numbers)
	print*, "matrix lambda", lambda_numbers(1:2) 
end program CargarDatos

subroutine MA(DIM,Nx,Ny, A)
	implicit none
    integer DIM,i,j,m,k,l,n,tot,Nx,Ny
    real A(DIM,DIM)
    real, parameter :: PI_sq = (3.14159264**2)
	real, parameter :: L1= 5.0, L2=5.0
	tot=1

    do i = 1,Nx
        do j = 1,Ny
            m = i+(j-1)*Ny
            do k = 1,Nx
                do l = 1,Ny
                    n = k+(l-1)*Ny

                    A(m,n)= PI_sq * ((m/L1)**2 + (n/L2)**2)

                enddo
            enddo
        enddo
    enddo

     return
end

pure subroutine gen_lambdas(DIM,lambda_numbers)
implicit none
integer, intent(in) :: DIM
  !! Spatial dimension of the vector space
  !! \( \mathcal{O}\subset \mathbb{R}^{d}\).
!real, intent(in) :: L1
!real, intent(in) :: L2
	  !!
real, intent(out):: lambda_numbers(DIM)
  !! Entries of matrix \( \Lambda \) from equation
  !!  $$
  !!    dU =
  !!      \big( -\beta \Lambda\,   -\theta A\,\big)  U   dt
  !!      + \sigma B\, U\, dW(t)
  !!  $$
real, parameter :: PI = 2.D0*DASIN(1.D0), L1=5.0, L2=5.0
real  pi_square, L1_res, L2_res, pi_square_sum_l1_res_l2_res, lambda_ij
integer i, j, k
integer, parameter :: Nx=10, Ny=10
pi_square = PI ** 2
L1_res = (L1) ** (-1)
L2_res = (L2) ** (-1)
pi_square_sum_l1_res_l2_res = (L1_res + L2_res) * (PI ** 2)
k=0
do i=1, Nx
  do j =1, Ny
	lambda_ij =  pi_square  * ((i/L1)**2 + (j/L1)**2)
	k = k+1
	lambda_numbers(k) =   lambda_ij
  enddo
enddo
return
end subroutine gen_lambdas

!ifort -fpp CargarDatos.f90 -o datos -lm -C 



