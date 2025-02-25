program multivariate
  implicit none
  integer seed,DIM,nobs,Nx,Ny,i,unit
  real theta,beta,gamma,sigma,delta,aum,PI
  parameter (seed=765431,Nx=10,Ny=10,nobs=1000,theta=0.5,beta=0.5,gamma=1.0,sigma=0.5,delta=0.001,aum=0.1,PI=3.1416)
  parameter(DIM=Nx*Ny)
  real x,G(DIM,DIM),A(DIM,DIM),path(0:nobs,DIM),lambdas(DIM),B(DIM,DIM),hs(DIM),startx(DIM),Ls(DIM)
  real Talpha(DIM,DIM),Tsigma(DIM,DIM),brownian(nobs,DIM),HT,times(0:nobs),AM(DIM*DIM)
  call srand(seed)

  open(99,file="MatrixA.dat")
        read(99,*) AM
    !    print*,'Dimension=',AM
  close(99)

! generate times
  call gen_times(nobs,delta,times)
! print*,"times",times

! generate lambas
  call gen_lambdas(DIM,aum,lambdas)
! print*,"lambdas",lambdas

  call MB(DIM,lambdas,gamma,B)
 !print*,"B",B

  call ML(DIM,lambdas,G)
! print*,"GAMMA",G

  call MA(DIM,Nx,Ny,AM,A)
! print*,"A",A

  call MDrift(DIM,theta,beta,G,A,Talpha)
! print*,"Talpha",Talpha


  call MDiff(DIM,sigma,B,Tsigma)
!print*,"Tsigma",Tsigma

  startx(:)=0.1

  call diffusion(DIM,Talpha,Tsigma,delta,startx,nobs,path,brownian)

  print*,"path1",path(:,1)
! Open a file for writing
  unit = 10  ! File unit number
  open(unit, file="vector_output.txt", status="replace", action="write")
! Write the vector to the file
  do i = 0, nobs
    write(unit, *) path(i,1)
  end do
! Close the file
  close(unit)
!print*,"path2",path(:,2)
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Module 1 para generar los parámetros del modelo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

    ! generate matrix LAMBDA
     subroutine ML(DIM,lambdas,G)
    implicit none
        integer DIM,i
    real lambdas(DIM),G(DIM,DIM)

        G(:,:)=0.0
        do i=1,DIM
         G(i,i)=lambdas(i)
          enddo

    return
    end


    ! generate matrix A
     subroutine MA(DIM,Nx,Ny,AM,A)
    implicit none
        integer DIM,i,j,m,k,l,n,tot,Nx,Ny
    real A(DIM,DIM),AM(DIM*DIM)


        tot=1
    do i = 1,Nx
        do j = 1,Ny
            m = i+(j-1)*Ny
            do k = 1,Nx

                do l = 1,Ny
                    n = k+(l-1)*Ny
                    A(m,n)=AM(tot)
                    tot=tot+1
                enddo
            enddo
        enddo
    enddo

    return
    end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Module 2 coeficientes de la Ecuación 10
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! generate matrix drift
     subroutine MDrift(DIM,theta,beta,G,A,Talpha)
    implicit none
        integer DIM,i,j
    real theta,beta,A(DIM,DIM),G(DIM,DIM),Talpha(DIM,DIM)


    Talpha(:,:)=beta*G(:,:)+theta*A(:,:)


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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module 3....Milstein for generate paths
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 !Using the Milstein scheme we simulate a diffusion from startx, n steps ahead at
! stepsizes delta. Diffusion depends on external functions alpha and sigma. brownian
! contains the original brownian increments used to construct the diffusion. brownian1
! contains the brownian increments of the inverse process, which is calculated from the


      subroutine diffusion(DIM,Talpha,Tsigma,delta,startx,nsteps,points,brownian)
      implicit none
      integer nsteps,i,j,d,DIM,m,k,errorflag
      real delta,startx(DIM),points(0:nsteps,DIM),y1(DIM),y2(DIM),x1(DIM),x2(DIM)
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


! Drift parameter

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
!!!! Diffusion parameter
subroutine DiffusionParameter(DIM,Tsigma,x,y)
    implicit none
    integer DIM,i,j
    real Tsigma(DIM,DIM),x(DIM),y(DIM)
         do i=1,DIM

            y(i)=Tsigma(i,i)*x(i)

         enddo
    return
    end


!  Generic routine to generate a uniform random number. Change the interior of the routine
!  in order to adopt your own method.

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



!


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
      real alpha(DIM),sigma(DIM),W(DIM),xx(DIM),sum,brown(DIM)

      do i=1,DIM
        xx(i)=0.0
      enddo

      call BrownianStep(DIM,delta,xx,W)

       do i=1,DIM
        sum=0.0

         sum=sum+sigma(i)*W(i)

    endx(i)=startx(i)+alpha(i)*delta+sum
        brown(i)=W(i)
      enddo

      return
      end




!
