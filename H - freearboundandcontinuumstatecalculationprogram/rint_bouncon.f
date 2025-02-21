 	program integration
c   Simple integration program
	implicit none
        real*8:: fn(30000),un(30000),Rn(30000)
	real*8:: rd(30000),RNorm(30000),fm(30000)
	real*8:: matelement(30000),uf(30000),ufp(30000)
	real*8:: Rf(30000),E(30000),UNorm(30000)
	real*8:: fo(30000),fq(30000)
	real*8:: h,d,Z,q,energy,result,ele
	real*8:: rd1,rd2,uf1,uf2,Rf1,Rf2
	INTEGER:: i,j,n1,n2,n3,l
	 	
c	RNorm is normalized 3p and Rm is normalized continuum d
c	Rn & un are  unnonrmalized 3p,uf & Rf normalized continuum
        n1=0
        n2=4357
	h=0.003
        d=0.001				
	E(1)=0.01
	Z=1.0
	q=1.0
        Open   (10, FILE = 'unormbwfn_total.dat')
        Open   (34, FIlE = 'bwfn_unorm_3p.dat')
	Open   (16, FIlE = 'bwfn_norm_R(3p)_result.dat')
	Open   (29, FIlE = 'bouncontinuum.dat')
	Open   (30, FIlE = 'r_continuum.dat')

c       normalized wavefunction
	
	do  i=n1,n2
	Read   (10,*)energy,rd(i),un(i),Rn(i)
	fn(i)=(rd(i)*rd(i)*Rn(i)*Rn(i))
	write(12,*) i,rd(i),energy,un(i),Rn(i),fn(i)
 	end do
        call rint(fn,n1,n2,7,h,result)
	write(34,*) result
	do  i=n1,n2
	RNorm(i)=Rn(i)/sqrt(result)
	UNorm(i)=rd(i)*RNorm(i)
	write(16,*) energy,rd(i),RNorm(i),UNorm(i)
	end do

c	calculate <R31/r/RE2> AND <R31/r/RE0>, E cotinuum energy
	
	do l=0,2,2
	
	do j=2,4990

 	E(j)=E(j-1)+d
	write(17,*)  j,E(j)

	call unbound(n2,uf,ufp,rd,dble(l),Z,q,E(j))
c	From unbound we can calculate unbounwfn

        do  i=n1,n2
	Rf(i)=uf(i)/rd(i)
	fm(i)=(rd(i)*rd(i)*rd(i)*RNorm(i)*Rf(i))
	fo(i)=(UNorm(i)*uf(i))
	fq(i)=(rd(i)*UNorm(i)*uf(i))
	if(j==910) then
	write(30,*) l,E(910),rd(i),UNorm(i),uf(i),fo(i),fq(i)
	end if
	end do
	call rint(fm,n1,n2,7,h,ele)
        matelement(j)=ele
	write(29,*) l,E(j),matelement(j)

	end do
	end do

	
        end program integration
