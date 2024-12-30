 	program Total_PShift_calculation

	implicit none
        real*8:: delta0(30000),delta2(30000)
	real*8:: sigma0(30000),sigma2(30000),fm(30000)
	real*8:: matele0(30000),matele2(30000),E(30000)
	real*8:: Re(30000),Im(30000),Thita(30000)
	INTEGER:: i,n1,n2,l0,l2
	 	
        n1=1
        n2=4989
	
        Open   (36, FILE = '3p_d.dat')
        Open   (39, FIlE = '3p_s.dat')
	Open   (26, FIlE = 'FrAr_del2_3p_d')
	Open   (19, FIlE = 'FrAr_del0_3p_s')
	Open   (16, FIlE = 'Zigms_3p_s_d.dat')
	Open   (12, FIlE = 'PS_fr_Ar.dat')

	
	do  i=n1,n2
	Read   (36,*) l2,E(i),matele2(i)
	Read   (39,*) l0,E(i),matele0(i)
	Read   (26,*) l2, E(i),delta2(i)
	Read   (19,*) l0, E(i),delta0(i)
	Read   (16,*) E(i),l0,sigma0(i),l2,sigma2(i)
	write(42,*) E(i),l0,matele0(i),delta0(i),sigma0(i)
	write(43,*) E(i),l2,matele2(i),delta2(i),sigma2(i)
	Re(i)=(cos(delta0(i)+sigma0(i))*(matele0(i)/2.0))
     1           -(cos(delta2(i)+sigma2(i))*matele2(i))
	Im(i)=(sin(delta0(i)+sigma0(i))*(matele0(i)/2.0))
     1           -(sin(delta2(i)+sigma2(i))*matele2(i))
	Thita(i)=Atan(Im(i)/Re(i))
	write(12,*) E(i),Thita(i)
 	end do
        end program Total_PShift_calculation
