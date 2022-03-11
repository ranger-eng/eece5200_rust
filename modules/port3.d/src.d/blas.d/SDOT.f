      REAL FUNCTION SDOT(N, SX, INCX, SY, INCY)
      INTEGER N, INCX, INCY
      REAL SX(INCX, 1), SY(INCY, 1)
      INTEGER I
C        RETURNS DOT PRODUCT OF SX AND SY
      SDOT = 0.
      IF (N .EQ. 0) RETURN
C/6S
      IF (N .LT. 0) CALL SETERR(17H  SDOT - N .LT. 0, 17, 1, 2)
      IF (INCX .LE. 0) CALL SETERR(20H  SDOT - INCX .LT. 0, 20, 1, 2)
      IF (INCY .LE. 0) CALL SETERR(20H  SDOT - INCY .LT. 0, 20, 1, 2)
C/7S
C     IF (N .LT. 0) CALL SETERR('  SDOT - N .LT. 0', 17, 1, 2)
C     IF (INCX .LE. 0) CALL SETERR('  SDOT - INCX .LT. 0', 20, 1, 2)
C     IF (INCY .LE. 0) CALL SETERR('  SDOT - INCY .LT. 0', 20, 1, 2)
C/
      DO  1 I = 1, N
         SDOT = SDOT+SX(1, I)*SY(1, I)
   1     CONTINUE
      RETURN
      END
