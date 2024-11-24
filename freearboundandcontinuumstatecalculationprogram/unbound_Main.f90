PROGRAM MAIN
  IMPLICIT NONE
  INTEGER, PARAMETER :: n = 100
  REAL*8 :: u(n), up(n), r(n), l, Z, q, energy
  INTEGER :: i
  
  l = 1.0d0
  Z = 1.0d0
  q = 0.0d0
  energy = 1.0d0

  DO i = 1,n
    r(i) = i*10.0d0/n
  END DO

  CALL unbound(n,u,up,r,l,Z,q,energy)

  DO i =1,n
    print*, r(i), u(i), up(i)
  END DO

END PROGRAM MAIN
