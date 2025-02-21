      program integration
c   Simple integration program
	implicit doubleprecision(a-h,o-z)
        dimension fn(20000),rd(20000),u(20000),R(20000),RNorm(20000)
c	dimension result(20000)
        n1=0
        n2=13044
	h=0.001
        Open   (1, FILE = 'bwfn_total.dat')
        Open   (34, FIlE = 'bwfn_norm.dat')
	Open   (16, FIlE = 'bwfn_norm_R(3p)_result.dat')
        do  i=n1,n2
        Read   (1,*)   E,rd(i),u(i),R(i)
        fn(i)=(rd(i)*rd(i)*R(i)*R(i))
	write(4,*) i,fn(i)
 	end do
        call rint(fn,n1,n2,7,h,result)
        write(34,*) result
	do  i=n1,n2
	RNorm(i)=R(i)/sqrt(result)
	write(16,*) RNorm(i)
	end do
        end program integration

