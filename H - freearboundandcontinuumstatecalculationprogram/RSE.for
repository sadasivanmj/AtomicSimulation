        subroutine RSE(r,u,udot,l,energy)
      implicit none
      integer neq
      parameter(neq = 2) 
      real*8 r
      real*8 u(neq),udot(neq)
      real*8 l,energy
      real*8 pot
      real*8 Z,a1,a2,a3,a4,a5,a6
      Z=1
        a1=16.0390
        a2=2.007
        a3=-25.543
        a4=4.525
        a5=0.961
        a6=0.443
       pot = -(Z+a1*exp(-a2*r)+a3*r*exp(-a4*r)+a5*exp(-a6*r))
      udot(1) = u(2) 

      udot(2) = l*(l+1.d0)/r + 2.d0*pot
C      if (r.ge.5.3 .AND. r.le.8.1) then
C      ndot(2)=ndot(2)-(0.6*r)
C      end if
      udot(2) = udot(2)/r - 2.d0*energy
      udot(2) = udot(2)*u(1)
c	write(6,*) "value of",r,pot,udot(2)	
      return
      end
