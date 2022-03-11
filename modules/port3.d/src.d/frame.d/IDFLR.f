      INTEGER FUNCTION IDFLR(X)
C
C  IDFLR RETURNS FLR(X)
C
      DOUBLE PRECISION X
C
      IDFLR = IDINT(X)
      IF (X .GE. 0.0D0) RETURN
      IF (DBLE(FLOAT(IDFLR)) .NE. X) IDFLR = IDFLR - 1
C
      RETURN
      END
