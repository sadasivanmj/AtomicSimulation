MODULE CSEED
        INTEGER :: ISEED
END MODULE CSEED

PROGRAM MAIN
        USE CSEED
        IMPLICIT NONE
        INTEGER :: I
        REAL :: RANDOM_VALUE

        DO I = 1,10
                RANDOM_VALUE=RANF()
                PRINT *, "RANDOM VALUE", I, RANDOM_VALUE
        END DO
END PROGRAM MAIN

FUNCTION RANF() 
        RESULT(R)
        USE CSEED
        IMPLICIT NONE
        INTEGER :: IH, IL, IT, IA, IC, IQ, IR
        DATA IA/16807/, IC/2147483647/, IQ/127773/, IR/2836/
        REAL :: R

        IH = ISEED/IQ
        IL = MOD(ISEED, IQ)
        IT = IA*IL-IR*IH
        IF(IT.GT.0) THEN
                ISEED = IT
        ELSE
                ISEED = IC + IT
        END IF
        R = ISEED/FLOAT(IC)
END FUNCTION RANF
