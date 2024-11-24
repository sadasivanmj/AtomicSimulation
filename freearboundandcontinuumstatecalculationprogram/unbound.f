      subroutine unbound(n,u,up,r,l,Z,q,energy)

      integer n,nnorm
      parameter (nnorm = 10)!10)
      real*8 u(n),up(n),r(n)
      real*8 l,Z,q,energy
      real*8 upara(0:3),uppara(0:3)
      real*8 uin(2),uout(2)
      real*8 A,Ap,App
      real*8 x(2),C(nnorm)
      real*8 alpha,norm,pihalf,relerr
      parameter (pihalf = 1.570796327d0)
      integer i,index

c     determination of expansion coefficients of wave function u(r) about r=0

      upara(0) = 1.d0
      upara(1) = -Z/(l+1.d0)
      upara(2) = Z*Z - (l+1.d0)*energy
      upara(2) = upara(2)/(l+1.d0)
      upara(2) = upara(2)/(2.d0*l+3.d0)
      upara(3) = Z*Z - (3.d0*l+4.d0)*energy
      upara(3) = -upara(3)*Z/3.d0
      upara(3) = upara(3)/(l+1.d0)
      upara(3) = upara(3)/(l+2.d0)
      upara(3) = upara(3)/(2.d0*l+3.d0)

c     determination of expansion coefficients of u'(r) about r=0

      uppara(0) = l+1.d0
      uppara(1) = upara(1)*(l+2.d0)
      uppara(2) = upara(2)*(l+3.d0)
      uppara(3) = upara(3)*(l+4.d0)

c     evaluation of unnormalized wave function at the first two grid points

      u(1) = 0.d0

      if (l.eq.0) then
         up(1) = 1.d0
      else
         up(1) = 0.d0
      end if

      u(2) = upara(3)
      u(2) = upara(2) + u(2)*r(2)
      u(2) = upara(1) + u(2)*r(2)
      u(2) = upara(0) + u(2)*r(2)
      u(2) = u(2)*r(2)**(l+1.d0)
      
      up(2) = uppara(3)
      up(2) = uppara(2) + up(2)*r(2)
      up(2) = uppara(1) + up(2)*r(2)
      up(2) = uppara(0) + up(2)*r(2)
      up(2) = up(2)*r(2)**l

c     numerical integration of radial Schroedinger equation
c        write(6,*)  "inside unbound program",l,Z,q,energy,n 
      do i = 2, n-1
	 uin(1) = u(i)
         uin(2) = up(i)
c        write(6,*)  "before rk4",i,l,Z,q,energy,n 
         call rk4(uin,r(i),uout,r(i+1),l,energy)
c        write(6,*)  "After rk4",i,l,Z,q,energy,n 
         u(i+1) = uout(1)
         up(i+1) = uout(2)
      end do

c     energy normalization

      upara(0) = 2.d0*energy
      upara(1) = l*(l+1.d0)
      upara(2) = 2.d0*q
      upara(3) = upara(1)*3.d0

      i = n - nnorm
      uppara(0) = upara(1)/r(i)
      uppara(1) = uppara(0) - q
      uppara(2) = uppara(1) - q
      uppara(3) = 3.d0*uppara(0) - 2.d0*q
      A = upara(0) - uppara(2)/r(i)
      Ap = 0.5d0*uppara(1)/r(i)**2
      Ap = Ap/A
      Ap = 5.d0*Ap**2
      App = 0.5d0*uppara(3)/r(i)**3
      App = App/A
      x(1) = dabs(A + Ap + App)
      x(1) = dsqrt(x(1))
      uin(1) = dsqrt(x(1))*u(i)

      norm = 0.d0
      do i = n - nnorm + 1, n
         uppara(0) = upara(1)/r(i)       
         uppara(1) = uppara(0) - q       
         uppara(2) = uppara(1) - q       
         uppara(3) = 3.d0*uppara(0) - 2.d0*q
         A = upara(0) - uppara(2)/r(i)
         Ap = 0.5d0*uppara(1)/r(i)**2
         Ap = Ap/A
         Ap = 5.d0*Ap**2
         App = 0.5d0*uppara(3)/r(i)**3
         App = App/A
         x(2) = dabs(A + Ap + App)
         x(2) = dsqrt(x(2))
         uin(2) = dsqrt(x(2))*u(i)
         alpha = 0.5d0*(x(1) + x(2))*(r(i) - r(i-1))
         index = n + 1 - i
         C(index) = uin(1)**2 + uin(2)**2 
         C(index) = C(index) - 2.d0*uin(1)*uin(2)*dcos(alpha)
         C(index) = C(index)*pihalf
         C(index) = dabs(dsqrt(C(index))/dsin(alpha))
c         write(98,*) index,r(i),C(index)
         norm = norm + C(index)
         x(1) = x(2) 
         uin(1) = uin(2)
      end do

      norm = norm/dble(nnorm)
      relerr = dabs(norm - C(1))/norm
      if (relerr .lt. 1.d-4) then
         norm = C(1)
      else
         write(*,*) 'normalization failed in unbound.f'
         stop
      end if
      norm = 1.d0/norm

      do i = 1, n
         u(i) = u(i)*norm
         up(i) = up(i)*norm
      end do
	
      return 

      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

