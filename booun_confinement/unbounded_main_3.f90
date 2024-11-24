      Program onephcoulmb
      !calculation of one-photon xsection using v-form. Final state wavefn is Coulmb waves
      Implicit none
      REAL*8:: rmax,dr,alpha,energy,pi,xsec,egyeV,k,elegy      
      REAL*8:: Z,q
      REAL*8:: M=0
      INTEGER:: l,lm
      REAL*8:: ui(3001),u(3001),r(3001),uf(3001),rd(3001) 
      REAL*8:: ufp(3001),p(3001),X(3001),jn(0:20),jnp(0:20)
      REAL*8:: b1s(3001)
      INTEGER:: i,j,ngrd
      ngrd=3001
      alpha=1/137.036d0
      rmax=26.d0
      dr = rmax/dble(ngrd-1)
      pi=3.14159265d0

      Z = 1.d0
      q = 1.d0
      l = 2.d0

      r(1) = 0.d0
      X(1) = 0.d0

C      open(unit=25,file='onecoulmb.dat')
      do j=1,1350

      energy=15.75d0+(j-1)*0.1!energy=j*0.5
      energy=energy/27.2113834d0
      elegy=energy-0.5d0

      !! Calculation of 1s  radial wavefunction of Hydrogen
      do i = 2, ngrd
      r(i) = r(i-1) + dr
      end do

C      open(unit=9,file='1s_interpolated.dat',form='formatted',
C     &   status='unknown')
C      DO i=1, ngrd
C      READ(9,*) rd(i),b1s(i)
C      END DO     

      !! To obtain continuum state function
      call unbound(ngrd,uf,ufp,r,dble(l),Z,q,elegy)
      do i=1,ngrd
           p(i)=uf(i)*rd(i)*b1s(i)
           write(28,*)rd(i),p(i)
           if(j.eq.660) then
           write(30,*)energy,r(i),uf(i)
           end if
      enddo
           write(29,*)energy,r(1617),uf(1617),r(1618),uf(1618)
       !!     M = trap(p,ngrd,dr)
C      do i=1,ngrd-1
C      M=M+(dr/2)*(p(i)+p(i+1))
C      enddo
C      xsec = (16.d0*pi**2*alpha*M**2)/(3.d0*energy)
C      xsec=xsec*0.52918D-8**2
C      egyeV=energy*27.2113834d0
C      close(9)     
C      write(25,*)egyeV,xsec
      enddo
     
C      close(25)
      end program onephcoulmb 
