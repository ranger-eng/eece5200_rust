      SUBROUTINE S88FMT( N, W, IFMT )
C
C  S88FMT  REPLACES IFMT(1), ... , IFMT(N) WITH
C  THE CHARACTERS CORRESPONDING TO THE N LEAST SIGNIFICANT
C  DIGITS OF W.
C
      INTEGER N,W
C/6S
      INTEGER IFMT(N)
C/7S
C     CHARACTER*1 IFMT(N)
C/
C
      INTEGER NT,WT
C
C/6S
      INTEGER DIGITS(10)
      DATA DIGITS( 1) / 1H0 /
      DATA DIGITS( 2) / 1H1 /
      DATA DIGITS( 3) / 1H2 /
      DATA DIGITS( 4) / 1H3 /
      DATA DIGITS( 5) / 1H4 /
      DATA DIGITS( 6) / 1H5 /
      DATA DIGITS( 7) / 1H6 /
      DATA DIGITS( 8) / 1H7 /
      DATA DIGITS( 9) / 1H8 /
      DATA DIGITS(10) / 1H9 /
C/7S
C     CHARACTER*1 DIGITS(10)
C     DATA DIGITS( 1) / '0' /
C     DATA DIGITS( 2) / '1' /
C     DATA DIGITS( 3) / '2' /
C     DATA DIGITS( 4) / '3' /
C     DATA DIGITS( 5) / '4' /
C     DATA DIGITS( 6) / '5' /
C     DATA DIGITS( 7) / '6' /
C     DATA DIGITS( 8) / '7' /
C     DATA DIGITS( 9) / '8' /
C     DATA DIGITS(10) / '9' /
C/
C
      NT = N
      WT = W
C
 10   IF (NT .LE. 0) RETURN
        IDIGIT = MOD( WT, 10 )
        IFMT(NT) = DIGITS(IDIGIT+1)
        WT = WT/10
        NT = NT - 1
        GO TO 10
C
      END
