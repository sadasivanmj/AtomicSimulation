PROGRAM RUNGEKUTTA
  IMPLICIT NONE
  
  INTEGER :: I,N
  REAL :: X(0:10)=0,Y(0:10)=0,H,K1,K2,F,K,X_I
  PRINT*, "ENTER X0, Y0, X, H"
  READ*, X(0),Y(0),X_I,H
  N=(X_I-X(0))/H

  DO I=1,N
  X(I) = X(I-1)
  Y(I) = Y(I-1)
  K1 = H*F(X(I),Y(I))
  K2 = (H*F(X(I)+H,Y(I)+K1))
  K = (K1+K2)/2.0
  Y(I) = Y(I)+K
  X(I) = X(I)+H
  PRINT*, "Y", I, "=", Y(I),"WHEN X = ",I,"=",X(I)
  END DO

END PROGRAM

REAL FUNCTION F(X,Y)
  IMPLICIT NONE
  REAL INTENT(IN) :: X,Y
  F=(X-Y)/2.0
END FUNCTION F

