C---- COULOMB DRIVER    (SIMPLE VERSION)        FILE = COUL90_S.FOR
C---- A.R.BARNETT       MANCHESTER              EASTER 1995 (REVISED)
C---- COMPUTES COULOMBS, SPHERICAL BESSELS, AND CYLINDRICAL BESSELS
C----          USES COUL90, SBESJY AND RICBES, AND COMPARES THEM
      PROGRAM  COUL90_S
      IMPLICIT         NONE
      DOUBLE PRECISION FC1(0:1200), GC1(0:1200), FCP1(0:1200)
      DOUBLE PRECISION GCP1(0:1200)
      DOUBLE PRECISION FC2(0:1200), GC2(0:1200), FCP2(0:1200)
      DOUBLE PRECISION GCP2(0:1200)
      DOUBLE PRECISION XJ(0:1200),  Y(0:1200), XJP(0:1200),  YP(0:1200)
      DOUBLE PRECISION PSI (0:1200), CHI (0:1200)
      DOUBLE PRECISION PSIP(0:1200), CHIP(0:1200)
      DOUBLE PRECISION X1,X2, ETA, XM, ZERO,Z,r1,r2,U1,U2
      DOUBLE PRECISION  K,E,DELTA,TDELTA,R1,R2
      INTEGER          LRANGE, KFN, IFAIL, M1, L
      CHARACTER        TEXT*72
      LOGICAL          SPHRIC
      ZERO = 0.D0
      OPEN  (1, FILE = 'COULSMPL1.IN')
      OPEN  (2, FILE = 'COULSMPL2.OUT')
      OPEN  (3, FILE = 'final_wave_r_14au.in') 
      OPEN  (4, FILE = 'PS_Sim_Ar_NUM.out') 
      WRITE (2, 9000)
      READ  (1, '(A)') TEXT
      READ  (1, 9010) XM, LRANGE, KFN
 100  CONTINUE
      READ  (3,*, END = 300) E,r1,U1,R1,r2,U2,R2
C       E is photo electron energy(already)
        K=sqrt(2*E)
        ETA=-(1/K)
	X1=K*r1
        X2=K*r2
	Do l=0,2,2
      CALL COUL90(X1, ETA, XM, LRANGE, FC1, GC1, FCP1, GCP1, KFN, IFAIL)
      CALL COUL90(X2, ETA, XM, LRANGE, FC2, GC2, FCP2, GCP2, KFN, IFAIL)
      IF (IFAIL .NE. 0) WRITE (2,10000) ' COUL90 ERROR! IFAIL=',IFAIL
           SPHRIC = (KFN.NE.2 .AND. ETA.EQ.ZERO)
C      COMPUTES DELTA(PHASE SHIFT)
        Z=((U2)/(U1))
        TDELTA=(Z*FC1(L)-FC2(L))/(GC2(L)-Z*GC1(L))
C        Z=((U2*r1)/(U1*r2))
C        TDELTA=(Z*FC1(2)-FC2(2))/(-GC2(2)+Z*GC1(2))
        DELTA=ATAN(TDELTA)
         WRITE(4,*) L,E, DELTA
	End Do
C      IF ( SPHRIC ) THEN
C         CALL SBESJY(X, LRANGE,  XJ,   Y,  XJP,   YP, IFAIL)
C      IF (IFAIL .NE. 0) WRITE (2,10000) ' SBESJY ERROR! IFAIL=',IFAIL       
C         CALL RICBES(X, LRANGE, PSI, CHI, PSIP, CHIP, IFAIL)
C      IF (IFAIL .NE. 0) WRITE (2,10000) ' RICBES ERROR! IFAIL=',IFAIL     
C      END IF
      IF      ( KFN.EQ.0 ) THEN
         WRITE (2, 9020) ETA, X1, XM, KFN, ' (Coulomb)'
         WRITE (2, 9020) ETA, X2, XM, KFN, ' (Coulomb)'
      ELSE IF ( KFN.EQ.1 ) THEN
C         WRITE (2, 9030)      X, XM, KFN, ' (SphBess)'
C      ELSE IF ( KFN.EQ.2 ) THEN
C         WRITE (2, 9040)      X, XM, KFN, ' (CylBess)'
      END IF
      IF (IFAIL .NE. 0)         GOTO 100          ! seek more data
      TEXT = 
     + ' The 3 sets of results are COUL90(KFN), SBESJY & (1/X) RICBES'
      IF ( SPHRIC ) WRITE (2, 9050) TEXT
      M1 = IDINT(XM)
      DO 200 L = M1, M1 + LRANGE, MAX(LRANGE,1)
C           WRITE  (2, 9060) L, FC(L), GC(L), FCP(L), GCP(L)
C           IF ( SPHRIC ) THEN
C            WRITE (2, 9060) L, XJ(L),  Y(L), XJP(L),  YP(L)
C            WRITE (2, 9060) L, PSI(L)/X, CHI(L)/X, PSIP(L)/X, CHIP(L)/X
C          END IF
  200    CONTINUE
                                GOTO 100          ! seek more data
C---- 
  300 CONTINUE
      CLOSE (1)
      CLOSE (2)
C---- 
 9000 FORMAT (8X, 
     +        ' TEST OF THE CONTINUED-FRACTION COULOMB & BESSEL', 
     +        ' PROGRAM - COUL90'//12X,
     +        ' WHEN LAMBDA IS AN INTEGER (L-VALUE) '//8X,
     +        '  F IS REGULAR AT THE ORIGIN ( X = 0 ) WHILE',/8X, 
     +        '  G IS IRREGULAR ( => -INFINITY AT X = 0 )'//4X,
     +        'L', 6X, ' F(ETA,X,L)', 6X, ' G(ETA,X,L)', 4X, 
     +        ' D/DX (F)', 9X, ' D/DX (G)'/)
 9010 FORMAT (1F10.3, 2I5)
 9020 FORMAT (1X, 'ETA =', F9.3, 4X, ' X = ', F10.3, 4X, 'XLMIN = ', 
     +        F6.2, 4X, 'KFN = ', I2, A10)
 9030 FORMAT (5X, 'j,y Bessels', 3X, ' X = ', F10.3, 4X, 'XLMIN = ', 
     +        F6.2, 4X, 'KFN = ', I2, A10)
 9040 FORMAT (5X, 'J,Y Bessels', 3X, ' X = ', F10.3, 4X, 'XLMIN = ', 
     +        F6.2, 4X, 'KFN = ', I2, A10)
 9050 FORMAT (7X, A70)
 9060 FORMAT (1X, I4, 1P4D17.7)
10000 FORMAT (1X, A20, I5)
      END
