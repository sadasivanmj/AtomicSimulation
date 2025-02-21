	program asymp_continuum_wfn
c   Simple integration program
	implicit none
        real*8:: un(30000),Rn(30000)
	real*8:: rd(30000)
	real*8:: uf(30000),ufp(30000)
	real*8:: Rf(30000),E(30000)
	real*8:: h,d,Z,q,energy
	real*8:: rd1,rd2,uf1,uf2,Rf1,Rf2
	INTEGER:: i,j,n1,n3,l
	 	
c	RNorm is normalized 3p and Rm is normalized continuum d
c	Rn & un are  unnonrmalized 3p,uf & Rf normalized continuum
  n1=0
	n3=8340
	h=0.003
  d=0.001				
	E(1)=0.01
	Z=18.0
	q=1.0
        Open   (10, FILE = 'unormbwfn_total.dat')

	Open   (31, FIlE = 'asymp_scatwfn.dat')
  Open(18, FILE = "normwfn.dat")

c       normalized wavefunction
	
	do  i=n1,n3
	Read   (10,*) energy,rd(i),un(i),Rn(i)
	write(11,*) i,rd(i)
 	end do
	
	do l=0,2,2
	
	do j=2,4990

 	E(j)=E(j-1)+d
	write(17,*)  j,E(j)

	call unbound(n3,uf,ufp,rd,dble(l),Z,q,E(j))
c	From unbound we can calculate unbounwfn

        do  i=n1,n3
	Rf(i)=uf(i)/rd(i)
	write(18,*) l,E(1990),rd(i),uf(i),Rf(i)
	end do

	rd1=rd(4704)
	rd2=rd(4740)
	uf1=uf(4704)
	uf2=uf(4740)
	Rf1=Rf(4704)
	Rf2=Rf(4740)

	write(31,*) l,E(j),rd1,uf1,Rf1,rd2,uf2,Rf2

       	end do
	end do
        end program asymp_continuum_wfn
