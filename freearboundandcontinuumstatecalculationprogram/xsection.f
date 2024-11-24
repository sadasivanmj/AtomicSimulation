	program Total_xsection_calculation

	implicit none
        real*8:: delta0(30000),delta2(30000)
	real*8:: sigma0(30000),sigma2(30000),fm(30000)
	real*8:: matele0(30000),matele2(30000),E(30000)
	real*8:: Re(30000),Im(30000),xsection(30000)
	real*8:: k(30000)
	INTEGER:: i,n1,n2,l0,l2
		
        n1=1
        n2=4989
	
        Open   (36, FILE = '3p_d.dat')
        Open   (39, FIlE = '3p_s.dat')
	Open   (26, FIlE = 'FrAr_del2_3p_d')
	Open   (19, FIlE = 'FrAr_del0_3p_s')
	Open   (16, FIlE = 'Zigms_3p_s_d.dat')
	Open   (12, FIlE = 'xsection_fr_Ar.dat')

	
	do  i=n1,n2
	Read   (36,*) l2,E(i),matele2(i)
	Read   (39,*) l0,E(i),matele0(i)
	Read   (26,*) l2, E(i),delta2(i)
	Read   (19,*) l0, E(i),delta0(i)
	Read   (16,*) E(i),l0,sigma0(i),l2,sigma2(i)
	write(42,*) E(i),l0,matele0(i),delta0(i),sigma0(i)
	write(43,*) E(i),l2,matele2(i),delta2(i),sigma2(i)
	k(i)=sqrt(2*E(i))

	xsection(i)=(((ABS(matele0(i)))**2/4.0)+(ABS(matele2(i)))**2)/k(i)
	
	
	write(12,*) E(i),xsection(i)
 	end do
        end program Total_xsection_calculation
