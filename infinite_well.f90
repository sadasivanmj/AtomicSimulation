       program main
        implicit real*8(a-h,o-z)
        parameter (km=10000)
        dimension z1(km),z2(km),x(km)    
        dimension a1(km),a2(km),a3(km),a4(km)
        dimension b1(km),b2(km),b3(km),b4(km)
        open(456,file="state_energy.txt",status="unknown")

        psi = 0.0d0
        dpsi = 1d-3
        rlambda = 10.0d0
        xi   = -rlambda
        xf   = +rlambda
        npts = 5000
        nstate = 0


        x(1)    = xi
        x(npts) = xf
        z1(1)   = psi
        z2(1)   = dpsi
        h = dabs(x(npts)-x(1))/dfloat(npts)
        
        do neps = 0,20000,1 
        eps = neps/1000.0d0

        do n  = 1,npts,1
        a1(n) = h*f1(x(n),z1(n),z2(n),eps)
        b1(n) = h*f2(x(n),z1(n),z2(n),eps)

        a2(n)=h*f1(x(n)+h/2.0d0,z1(n)+a1(n)/2.0d0,z2(n)+b1(n)/2.0d0,eps)
        b2(n)=h*f2(x(n)+h/2.0d0,z1(n)+a1(n)/2.0d0,z2(n)+b1(n)/2.0d0,eps)

        a3(n)=h*f1(x(n)+h/2.0d0,z1(n)+a2(n)/2.0d0,z2(n)+b2(n)/2.0d0,eps)
        b3(n)=h*f2(x(n)+h/2.0d0,z1(n)+a2(n)/2.0d0,z2(n)+b2(n)/2.0d0,eps)

        a4(n)=h*f1(x(n)+h,z1(n)+a3(n),z2(n)+b3(n),eps)
        b4(n)=h*f2(x(n)+h,z1(n)+a3(n),z2(n)+b3(n),eps)

        x(n+1)  =   x(n) + h
        z1(n+1) =   z1(n) + (a1(n)+ 2.0d0*(a2(n)+a3(n)) + a4(n))/6.0d0
        z2(n+1) =   z2(n) + (b1(n)+ 2.0d0*(b2(n)+b3(n)) + b4(n))/6.0d0
        end do
        if ((dabs(z1(n))).lt.1d-7) then
        nstate = nstate+1 
        write(456,*)nstate,eps
        end if
        end do
        end program main



        function f1(x,z1,z2,eps)
        implicit real*8(a-h,o-z)
        f1 = z2
        end function

        function f2(x,z1,z2,eps)
        implicit real*8(a-h,o-z)
        rmass = 938.272
        hbarc = 197.330
        f2 = -2.0*rmass*eps*z1/hbarc**2
        end function
