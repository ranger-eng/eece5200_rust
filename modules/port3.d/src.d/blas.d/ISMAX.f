        INTEGER FUNCTION ISMAX(N,X,INCX)
        INTEGER I,N,INCX
        REAL X(INCX,1),SMAX
C
C THIS FUNCTION RETURNS THE INDEX OF THE LARGEST (ALGEBRAIC)
C COMPONENT OF X.
C ONLY EVERY INCXTH COMPONENT OF X IS CONSIDERED.
C
        ISMAX=0
        IF(N.EQ.0) RETURN
C/6S
        IF(N.LT.0) CALL SETERR(12HISMAX-N.LT.0,12,1,2)
        IF(INCX.LT.1)CALL SETERR(15HISMAX-INCX.LT.1,15,2,2)
C/7S
C       IF(N.LT.0) CALL SETERR('ISMAX-N.LT.0',12,1,2)
C       IF(INCX.LT.1)CALL SETERR('ISMAX-INCX.LT.1',15,2,2)
C/
        SMAX=X(1,1)
        ISMAX=1
        DO 10 I=1,N
        IF(SMAX.GE.X(1,I)) GO TO 10
           SMAX=X(1,I)
           ISMAX=I
 10     CONTINUE
        RETURN
        END
