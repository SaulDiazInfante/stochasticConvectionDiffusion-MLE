program multivariate
  implicit none
  integer seed,DIM,nobs
  real :: theta, beta, gamma, sigma, delta, aum, PI
  parameter(seed=765431,DIM=2,nobs=10,theta=0.5,beta=0.5,gamma=1.0,sigma=0.1,delta=0.01,aum=0.1,PI=3.1416)
  real x, G(DIM,DIM), A(DIM,DIM), path(0:nobs,DIM), lambdas(DIM),B(DIM,DIM)
  real hs(DIM),startx(DIM),Ls(DIM)
  real Talpha(DIM,DIM),Tsigma(DIM,DIM),brownian(nobs,DIM),HT,times(0:nobs)
  call srand(seed)
  ! generate times
  call gen_times(nobs,delta,times)
  ! print*,"times",times
  ! generate lambas
  call gen_lambdas(DIM,aum,lambdas)
  print*,"lambdas",lambdas
  call MB(DIM,lambdas,gamma,B)
  print*,"B",B
  call MG(DIM,lambdas,G)
  print*,"GAMMA",G
  call gen_hs(DIM,x,PI,Ls,hs)
  print*,"hs",hs
  call MA(DIM,hs,A)
  print*,"A",A
  call MDrift(DIM,theta,beta,G,A,Talpha)
  print*,"Talpha",Talpha
  call MDiff(DIM,sigma,B,Tsigma)
  print*,"Tsigma",Tsigma
  startx(1)=1.0
  startx(2)=1.0
  call diffusion(DIM,Talpha,Tsigma,delta,startx,nobs,path,brownian)
  print*,"path1",path(:,1)
  print*,"path2",path(:,2)
end
!!  Generic routine to generate a uniform random number. Change the interior of !!  the routine
!!  in order to adopt your own method.

      FUNCTION UNIF(IX)
       call random_number(x)
       UNIF=x
      RETURN
      END


!  On increment in a Brownian motion. In dimension DIM, from initial point startx
!  a point endx is generated a time delta ahead. The result endx of the routine is a
!  point delta time ahead in the Brownian motion.


    subroutine BrownianStep(DIM,delta,startx,endx)
	implicit none
        integer DIM,i
	real delta,startx(DIM),endx(DIM),var

        do i=1,DIM
        call normalvar(var)
        endx(i)=startx(i)+sqrt(delta)*var
        enddo

	return
	end

    ! generate times
     subroutine gen_times(nobs,delta,times)
    implicit none
        integer nobs,i
    real delta,times(0:nobs)


        do i=0,nobs
         times(i)=i*delta
        enddo

    return
    end
! generate lambas
     subroutine gen_lambdas(DIM,aum,lambdas)
    implicit none
        integer DIM,i
    real aum,lambdas(DIM)

        lambdas(1)=0.0
        do i=2,DIM
         lambdas(i)=aum+lambdas(i-1)
        enddo

    return
    end

  ! generate matrix B
  subroutine MB(DIM,lambdas,gamma,B)
    implicit none
        integer DIM,i
    real lambdas(DIM),B(DIM,DIM),gamma

        B(:,:)=0.0
        do i=1,DIM
         B(i,i)=lambdas(i)**(-gamma)
            enddo
     B(1,1)=1.0
    return
    end

    ! generate matrix GAMMA
     subroutine MG(DIM,lambdas,G)
    implicit none
        integer DIM,i
    real lambdas(DIM),G(DIM,DIM)

        G(:,:)=0.0
        do i=1,DIM
         G(i,i)=lambdas(i)
          enddo

    return
    end

    ! generate hs
  subroutine gen_hs(DIM,x,PI,Ls,hs)
    implicit none
        integer DIM,i,j
    real hs(DIM),PI,Ls(DIM),x



        do i=1,DIM
         hs(i)=COS(PI*i*x/Ls(i))
        enddo

    return
    end

    ! generate matrix A
     subroutine MA(DIM,hs,A)
    implicit none
        integer DIM,i
    real hs(DIM),A(DIM,DIM)


        do i=1,DIM
         A(i,:)=hs(:)
          enddo

    return
    end

     ! generate matrix drift
     subroutine MDrift(DIM,theta,beta,G,A,Talpha)
    implicit none
        integer DIM,i,j
    real theta,beta,A(DIM,DIM),G(DIM,DIM),Talpha(DIM,DIM)


    Talpha(:,:)=-beta*G(:,:)-theta*A(:,:)


    return
    end

        ! generate matrix diffusion
     subroutine MDiff(DIM,sigma,B,Tsigma)
    implicit none
        integer DIM
    real sigma,B(DIM,DIM),Tsigma(DIM,DIM)



         Tsigma(:,:)=sigma*B(:,:)


    return
    end

! Brownian motion routine in one dimension

        subroutine BM1(numsteps,delta,points)
        implicit none
        integer i,numsteps
	real delta,points(0:numsteps)
        points(0)=0.0
        do i=1,numsteps
           call BrownianStep(1,delta,points(i-1),points(i))
        enddo
        return
        end



! Generic routine to generate a standard normal variate. Change the interior to
! use a specific method. Box-Muller is used here.

        subroutine normalvar(x)
        real x,u
        integer seed
        common seed
        x=boxmuller(seed)
        return
        end

!  Box-Muller method for generating a normal random variate. Uses seed dummy
!  which is not of any importance unless a method for random number generation
!  in unif is used with a seed (which can be the case). In the default case, unif
!  has seed as an input but does not currently use it.

        function boxmuller(seed)
        integer seed
        real boxmuller,u
        integer iset
        real fac,gset,rsq,v1,v2,unif
        save iset,gset
        data iset/0/

 1      if (iset.eq.0) then
         v1=2.*unif(seed)-1.
         v2=2.*unif(seed)-1.
         rsq=v1**2+v2**2
         if (rsq.ge.1..or.rsq.eq.0)goto 1
         fac=sqrt(-2.*log(rsq)/rsq)
         gset=v1*fac
         boxmuller=v2*fac
         iset=1
        else
         boxmuller=gset
         iset=0
        endif
        return
        end

! Milstein calculates for stepsize delta the next point in a diffusion.
! startx is the current location of the diffusion, alpha=alpha(startx) is
! the drift evaluated at startx and sigma=sigma(startx) is the diffusion
! coefficient at startx. endx is the next point in the diffusion (output).
! alpha and sigma must be available by calling this routine.

      subroutine MilsteinStep(DIM,delta,startx,alpha,sigma,endx,brown)
      implicit none
      integer DIM,i,j
      real delta,startx(DIM),endx(DIM)
      real alpha(DIM),sigma(DIM,DIM),W(DIM),xx(DIM),sum,brown(DIM)

      do i=1,DIM
        xx(i)=0.0
      enddo

      call BrownianStep(DIM,delta,xx,W)

       do i=1,DIM
        sum=0.0
       do j=1,DIM
         sum=sum+sigma(i,j)*W(j)
       enddo
	endx(i)=startx(i)+alpha(i)*delta+sum
        brown(i)=W(i)
      enddo

      return
      end




! Using the Milstein scheme we simulate a diffusion from startx, n steps ahead at
! stepsizes delta. Diffusion depends on external functions alpha and sigma. brownian
! contains the original brownian increments used to construct the diffusion. brownian1
! contains the brownian increments of the inverse process, which is calculated from the


      subroutine diffusion(DIM,Talpha,Tsigma,delta,startx,nsteps,points,brownian)
      implicit none
      integer nsteps,i,j,d,DIM,m,k,errorflag
      real delta,startx(DIM),points(0:nsteps,DIM),y1(DIM),y2(DIM,DIM),x1(DIM),x2(DIM)
      real brownian(nsteps,DIM),brownian1(nsteps,DIM),brown(DIM)
      real Talpha(DIM,DIM),Tsigma(DIM,DIM),invmat(DIM,DIM),v1(DIM),v2(DIM)
! nsteps replaces 10000000
      do i=1,DIM
	points(0,i)=startx(i)
        x1(i)=startx(i)
      enddo

      do i=1,nsteps
       call DriftParameter(DIM,Talpha,x1,y1)
       call DiffusionParameter(DIM,Tsigma,x1,y2)
       call MilsteinStep(DIM,delta,x1,y1,y2,x2,brown)

       do j=1,DIM
        points(i,j)=x2(j)
        x1(j)=x2(j)
        brownian(i,j)=brown(j)
       enddo
      enddo



	return
	end


! Evaluates y=mu(x;theta) where theta is a parameter vector of length numpara

      subroutine DriftParameter(DIM,Talpha,x,y)
      implicit none
      integer i,j,DIM
      real Talpha(DIM,DIM),x(DIM),y(DIM),sum
       do i=1,DIM
        sum=0.0
        do j=1,DIM
          sum=sum+Talpha(i,j)*x(j)
        enddo
        y(i)=-sum
       enddo
      return
      end

	subroutine DiffusionParameter(DIM,Tsigma,x,y)
	implicit none
	integer DIM,i,j
	real Tsigma(DIM,DIM),x(DIM),y(DIM,DIM)
         do i=1,DIM
          do j=1,DIM
            y(i,j)=Tsigma(i,j)
          enddo
         enddo
	return
	end

     SUBROUTINE FINDInv(matrix, inverse, n, errorflag)
        IMPLICIT NONE
        !Declarations
        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
        REAL, INTENT(IN), DIMENSION(n,n) :: matrix  !Input matrix
        REAL, INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix

        LOGICAL :: FLAG = .TRUE.
        INTEGER :: i, j, k, l
        REAL :: m
        REAL, DIMENSION(n,2*n) :: augmatrix !augmented matrix

        !Augment input matrix with an identity matrix
        DO i = 1, n
                DO j = 1, 2*n
                        IF (j <= n ) THEN
                                augmatrix(i,j) = matrix(i,j)
                        ELSE IF ((i+n) == j) THEN
                                augmatrix(i,j) = 1
                        Else
                                augmatrix(i,j) = 0
                        ENDIF
                END DO
        END DO

        !Reduce augmented matrix to upper traingular form
        DO k =1, n-1
                IF (augmatrix(k,k) == 0) THEN
                        FLAG = .FALSE.
                        DO i = k+1, n
                                IF (augmatrix(i,k) /= 0) THEN
                                        DO j = 1,2*n
                                                augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
                                        END DO
                                        FLAG = .TRUE.
                                        EXIT
                                ENDIF
                                IF (FLAG .EQV. .FALSE.) THEN
                                        PRINT*, "Matrix is non - invertible"
                                        inverse = 0
                                        errorflag = -1
                                        return
                                ENDIF
                        END DO
                ENDIF
                DO j = k+1, n
                        m = augmatrix(j,k)/augmatrix(k,k)
                        DO i = k, 2*n
                                augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
                        END DO
                END DO
        END DO

        !Test for invertibility
        DO i = 1, n
                IF (augmatrix(i,i) == 0) THEN
                        PRINT*, "Matrix is non - invertible"
                        inverse = 0
                        errorflag = -1
                        return
                ENDIF
        END DO

        !Make diagonal elements as 1
        DO i = 1 , n
                m = augmatrix(i,i)
                DO j = i , (2 * n)
                           augmatrix(i,j) = (augmatrix(i,j) / m)
                END DO
        END DO

        !Reduced right side half of augmented matrix to identity matrix
        DO k = n-1, 1, -1
                DO i =1, k
                m = augmatrix(i,k+1)
                        DO j = k, (2*n)
                                augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
                        END DO
                END DO
        END DO

        !store answer
        DO i =1, n
                DO j = 1, n
                        inverse(i,j) = augmatrix(i,j+n)
                END DO
        END DO
        errorflag = 0
        END SUBROUTINE FINDinv


