        subroutine rk4(u0,r0,uf,rf,l,energy)

      implicit none

      integer neq
      parameter (neq = 2)
      real*8 u0(neq),uf(neq) 
      real*8 r0,rc,rf,dr
      real*8 l,energy
      real*8 work(neq,5)
      integer i

      rc = 0.5d0*(r0 + rf)
      dr = rf - r0
c        write(6,*)  "in rk4",neq,i,l,energy

      call RSE(r0,u0,work(1,1),l,energy)

      do i = 1, neq
         work(i,5) = u0(i) + 0.5d0*dr*work(i,1)
      end do

      call RSE(rc,work(1,5),work(1,2),l,energy)

      do i = 1, neq
         work(i,5) = u0(i) + 0.5d0*dr*work(i,2)
      end do
 
      call RSE(rc,work(1,5),work(1,3),l,energy)

      do i = 1, neq
         work(i,5) = u0(i) + dr*work(i,3)
      end do

      call RSE(rf,work(1,5),work(1,4),l,energy)

      do i = 1, neq
         work(i,5) = work(i,1) + 2.d0*work(i,2) 
     &        + 2.d0*work(i,3) + work(i,4)
         uf(i) = u0(i) + work(i,5)*dr/6.d0
c	 write(6,*)  "in rk4",neq,i,l,uf(i),energy
      end do

      return 
      end
