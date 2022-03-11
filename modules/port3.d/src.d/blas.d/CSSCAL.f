       SUBROUTINE CSSCAL(N,A,X,INCX)
C THIS SUBROUTINE SCALES EVERY INCXTH COMPONENT OF X
C BY A

      COMPLEX X(INCX,1)
       REAL A
      IF(N.EQ.0) RETURN

C/6S
      IF(N.LT.0) CALL SETERR(13HCSSCAL-N.LT.0,13,1,2)
      IF(INCX.LT.1)CALL SETERR(16HCSSCAL-INCX.LT.1,16,2,2)
C/7S
C     IF(N.LT.0) CALL SETERR('CSSCAL-N.LT.0',13,1,2)
C     IF(INCX.LT.1)CALL SETERR('CSSCAL-INCX.LT.1',16,2,2)
C/
c========CT MOD
      DO 10 I=1,N
         X(1,I)=X(1,I)*CMPLX(A,0.0)
 10   CONTINUE
c	iarg = 1-incx
c	do i=1,N
c	x(iarg+incx*i) = x(iarg+incx*i)*cmplx(A,0.0)
c	enddo
c========
      RETURN
      END
