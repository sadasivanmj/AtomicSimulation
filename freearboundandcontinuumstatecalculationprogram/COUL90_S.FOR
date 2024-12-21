C---- COULOMB DRIVER    (SIMPLE VERSION)        FILE = COUL90_S.FOR
C---- A.R.BARNETT       MANCHESTER              EASTER 1995 (REVISED)
C---- COMPUTES COULOMBS, SPHERICAL BESSELS, AND CYLINDRICAL BESSELS
C----          USES COUL90, SBESJY AND RICBES, AND COMPARES THEM
      PROGRAM  COUL90_S
      IMPLICIT         NONE
      DOUBLE PRECISION FC(0:1200), GC(0:1200), FCP(0:1200), GCP(0:1200)
      DOUBLE PRECISION XJ(0:1200),  Y(0:1200), XJP(0:1200),  YP(0:1200)
      DOUBLE PRECISION PSI (0:1200), CHI (0:1200)
      DOUBLE PRECISION PSIP(0:1200), CHIP(0:1200)
      DOUBLE PRECISION X, ETA, XM, ZERO
      INTEGER          LRANGE, KFN, IFAIL, M1, L
      CHARACTER        TEXT*72
      LOGICAL          SPHRIC
      ZERO = 0.D0
      OPEN  (1, FILE = 'COULSMPL.IN')
      OPEN  (2, FILE = 'COULSMPL.OUT')
      WRITE (2, 9000)
      READ  (1, '(A)') TEXT
 100  CONTINUE
      READ  (1, 9010, END = 300) ETA, X, XM, LRANGE, KFN
      CALL COUL90(X, ETA, XM, LRANGE, FC, GC, FCP, GCP, KFN, IFAIL)
      IF (IFAIL .NE. 0) WRITE (2,10000) ' COUL90 ERROR! IFAIL=',IFAIL
           SPHRIC = (KFN.NE.2 .AND. ETA.EQ.ZERO)
      IF ( SPHRIC ) THEN
         CALL SBESJY(X, LRANGE,  XJ,   Y,  XJP,   YP, IFAIL)
      IF (IFAIL .NE. 0) WRITE (2,10000) ' SBESJY ERROR! IFAIL=',IFAIL       
         CALL RICBES(X, LRANGE, PSI, CHI, PSIP, CHIP, IFAIL)
      IF (IFAIL .NE. 0) WRITE (2,10000) ' RICBES ERROR! IFAIL=',IFAIL     
      END IF
      IF      ( KFN.EQ.0 ) THEN
         WRITE (2, 9020) ETA, X, XM, KFN, ' (Coulomb)'
      ELSE IF ( KFN.EQ.1 ) THEN
         WRITE (2, 9030)      X, XM, KFN, ' (SphBess)'
      ELSE IF ( KFN.EQ.2 ) THEN
         WRITE (2, 9040)      X, XM, KFN, ' (CylBess)'
      END IF
      IF (IFAIL .NE. 0)         GOTO 100          ! seek more data
      TEXT = 
     + ' The 3 sets of results are COUL90(KFN), SBESJY & (1/X) RICBES'
      IF ( SPHRIC ) WRITE (2, 9050) TEXT
      M1 = IDINT(XM)
      DO 200 L = M1, M1 + LRANGE, MAX(LRANGE,1)
           WRITE  (2, 9060) L, FC(L), GC(L), FCP(L), GCP(L)
           IF ( SPHRIC ) THEN
            WRITE (2, 9060) L, XJ(L),  Y(L), XJP(L),  YP(L)
            WRITE (2, 9060) L, PSI(L)/X, CHI(L)/X, PSIP(L)/X, CHIP(L)/X
           END IF
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
 9010 FORMAT (3F10.3, 2I5)
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
$INCLUDE: 'COUL90.FOR'
$INCLUDE: 'RICBES.FOR'
$INCLUDE: 'SBESJY.FOR'
