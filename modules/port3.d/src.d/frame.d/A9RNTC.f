      SUBROUTINE A9RNTC(A, NITEMS, IOUT, MCOL, W, D)
C
C  THIS IS THE DOCUMENTED ROUTINE APRNTC, BUT WITHOUT THE CALLS TO
C  SETERR- BECAUSE IT IS CALLED BY SETERR.
C
C  THIS SUBROUTINE PRINTS OUT NITEMS FROM THE COMPLEX ARRAY, A, ON
C  OUTPUT UNIT IOUT, USING A MAXIMUM OF MCOL PRINT SPACES.
C  THE OUTPUT FORMAT IS 2(1PEW.D).
C  THE PROGRAM PUTS AS MANY VALUES ON A LINE AS POSSIBLE.
C  W SHOULD BE INPUT AS THE ACTUAL WIDTH +1 FOR A SPACE BETWEEN VALUES.
C
C  DUPLICATE LINES ARE NOT ALL PRINTED, BUT ARE INDICATED BY ASTERISKS.
C
C  WRITTEN BY DAN WARNER, REVISED BY PHYL FOX, OCTOBER 21, 1982.
C
C  THE LINE WIDTH IS COMPUTED AS THE MINIMUM OF THE INPUT MCOL AND 160.
C  IF THE LINE WIDTH IS TO BE INCREASED ABOVE 160, THE BUFFERS LINE()
C  AND LAST(), WHICH THE VALUES TO BE PRINTED ON ONE LINE, MUST
C  BE DIMENSIONED ACCORDINGLY.
C
C  INPUT PARAMETERS -
C
C    A        - THE START OF THE COMPLEX ARRAY TO BE PRINTED
C
C    NITEMS   - THE NUMBER OF ITEMS TO BE PRINTED
C
C    IOUT     - THE OUTPUT UNIT FOR PRINTING
C
C    MCOL     - THE NUMBER OF SPACES ACROSS THE LINE
C
C    W        - THE WIDTH OF THE PRINTED VALUE (1PEW.D)
C
C    D        - THE NUMBER OF DIGITS AFTER THE DECIMAL POINT (1PEW.D)
C
C
C  ERROR STATES - NONE.  LOWER LEVEL ROUTINE CALLED BY
C  SETERR, SO IT CANNOT CALL SETERR.
C
      INTEGER  NITEMS, IOUT, MCOL, W, D
C/R
C     REAL A(2,NITEMS)
C/C
      COMPLEX  A(NITEMS)
C/
C
      INTEGER  MAX0, MIN0, WW, DD, EMIN, EMAX,
     1         EXPENT, I1MACH, ICEIL, IABS, I10WID
C/6S
      INTEGER  IFMT1(20), IFMT2(18), BLANK, STAR
      INTEGER IFMT1C(20), IFMT2C(18)
      EQUIVALENCE (IFMT1(1),IFMT1C(1)), (IFMT2(1),IFMT2C(1))
C/7S
C     CHARACTER*1  IFMT1(20), IFMT2(18), BLANK, STAR
C     CHARACTER*20 IFMT1C
C     CHARACTER*18 IFMT2C
C     EQUIVALENCE (IFMT1(1),IFMT1C), (IFMT2(1),IFMT2C)
C/
      INTEGER  INDW, NCOL, COUNT, I, J, K, ILINE, ILAST
      LOGICAL  DUP
C/R
C     REAL LINE(2,18), LAST(2,18)
C/C
      COMPLEX  LINE(18), LAST(18)
C/
      REAL  LOGETA
C
C/6S
      DATA BLANK/1H /, STAR/1H*/, INDW/7/, EXPENT/0/
C/7S
C     DATA BLANK/' '/, STAR/'*'/, INDW/7/, EXPENT/0/
C/
C
C  IFMT1 IS FOR THE ASTERISK LINES- IFMT2 FOR THE DATA LINES
C
C/6S
      DATA IFMT1( 1) /1H(/,  IFMT2( 1) /1H(/
      DATA IFMT1( 2) /1H1/,  IFMT2( 2) /1H1/
      DATA IFMT1( 3) /1HA/,  IFMT2( 3) /1HA/
      DATA IFMT1( 4) /1H1/,  IFMT2( 4) /1H1/
      DATA IFMT1( 5) /1H,/,  IFMT2( 5) /1H,/
      DATA IFMT1( 6) /1H5/,  IFMT2( 6) /1HI/
      DATA IFMT1( 7) /1HX/,  IFMT2( 7) /1H7/
      DATA IFMT1( 8) /1H,/,  IFMT2( 8) /1H,/
      DATA IFMT1( 9) /1H2/,  IFMT2( 9) /1H1/
      DATA IFMT1(10) /1HA/,  IFMT2(10) /1HP/
      DATA IFMT1(11) /1H1/,  IFMT2(11) /1H /
      DATA IFMT1(12) /1H,/,  IFMT2(12) /1HE/
      DATA IFMT1(13) /1H /,  IFMT2(13) /1H /
      DATA IFMT1(14) /1H /,  IFMT2(14) /1H /
      DATA IFMT1(15) /1HX/,  IFMT2(15) /1H./
      DATA IFMT1(16) /1H,/,  IFMT2(16) /1H /
      DATA IFMT1(17) /1H2/,  IFMT2(17) /1H /
      DATA IFMT1(18) /1HA/,  IFMT2(18) /1H)/
      DATA IFMT1(19) /1H1/
      DATA IFMT1(20) /1H)/
C/7S
C     DATA IFMT1( 1) /'('/,  IFMT2( 1) /'('/
C     DATA IFMT1( 2) /'1'/,  IFMT2( 2) /'1'/
C     DATA IFMT1( 3) /'A'/,  IFMT2( 3) /'A'/
C     DATA IFMT1( 4) /'1'/,  IFMT2( 4) /'1'/
C     DATA IFMT1( 5) /','/,  IFMT2( 5) /','/
C     DATA IFMT1( 6) /'5'/,  IFMT2( 6) /'I'/
C     DATA IFMT1( 7) /'X'/,  IFMT2( 7) /'7'/
C     DATA IFMT1( 8) /','/,  IFMT2( 8) /','/
C     DATA IFMT1( 9) /'2'/,  IFMT2( 9) /'1'/
C     DATA IFMT1(10) /'A'/,  IFMT2(10) /'P'/
C     DATA IFMT1(11) /'1'/,  IFMT2(11) /' '/
C     DATA IFMT1(12) /','/,  IFMT2(12) /'E'/
C     DATA IFMT1(13) /' '/,  IFMT2(13) /' '/
C     DATA IFMT1(14) /' '/,  IFMT2(14) /' '/
C     DATA IFMT1(15) /'X'/,  IFMT2(15) /'.'/
C     DATA IFMT1(16) /','/,  IFMT2(16) /' '/
C     DATA IFMT1(17) /'2'/,  IFMT2(17) /' '/
C     DATA IFMT1(18) /'A'/,  IFMT2(18) /')'/
C     DATA IFMT1(19) /'1'/
C     DATA IFMT1(20) /')'/
C/
C
C     EXPENT IS USED AS A FIRST-TIME SWITCH TO SIGNAL IF THE
C     MACHINE-VALUE CONSTANTS HAVE BEEN COMPUTED.
C
      IF (EXPENT .GT. 0) GO TO 10
         LOGETA = ALOG10(FLOAT(I1MACH(10)))
         EMIN   = ICEIL(LOGETA*FLOAT(IABS(I1MACH(12)-1)))
         EMAX   = ICEIL(LOGETA*FLOAT(I1MACH(13)))
         EXPENT = I10WID(MAX0(EMIN, EMAX))
C
C     COMPUTE THE FORMATS.
C
   10 WW = MIN0(99, MAX0(W, 5+EXPENT))
      CALL S88FMT(2, WW, IFMT2(13))
      DD = MIN0(D, (WW-(5+EXPENT)))
      CALL S88FMT(2, DD, IFMT2(16))
C
C  NCOL IS THE NUMBER OF VALUES TO BE PRINTED ACROSS THE LINE.
C
      NCOL = MAX0(1, MIN0(9, (MIN0(MCOL,160)-INDW)/(2*WW)))
      CALL S88FMT(1, (2*NCOL), IFMT2(11))
      WW = WW-2
C
C  THE ASTERISKS ARE POSITIONED RIGHT-ADJUSTED IN THE W-WIDTH SPACE.
      CALL S88FMT(2, WW, IFMT1(13))
C
C  I COUNTS THE NUMBER OF ITEMS TO BE PRINTED,
C  J COUNTS THE NUMBER ON A GIVEN LINE,
C  COUNT COUNTS THE NUMBER OF DUPLICATE LINES.
C
      I = 1
      J = 0
      COUNT = 0
C
C  THE LOGICAL OF THE FOLLOWING IS ROUGHLY THIS -
C  IF THERE ARE STILL MORE ITEMS TO BE PRINTED, A LINE-
C  FULL IS PUT INTO THE ARRAY, LINE.
C  WHENEVER A LINE IS PRINTED OUT, IT IS ALSO STUFFED INTO
C  THE ARRAY, LAST, TO COMPARE WITH THE NEXT ONE COMING IN
C  TO CHECK FOR REPEAT OR DUPLICATED LINES.
C  ALSO WHENEVER A LINE IS WRITTEN OUT, THE DUPLICATION
C  COUNTER, COUNT, IS SET TO ONE.
C  THE ONLY MILDLY TRICKY PART IS TO NOTE THAT COUNT HAS TO
C  GO TO 3 BEFORE A LINE OF ASTERISKS IS PRINTED BECAUSE
C  OF COURSE NO SUCH LINE IS PRINTED FOR JUST A PAIR OF
C  DUPLICATE LINES.
C
C  ILINE IS PRINTED AS THE INDEX OF THE FIRST ARRAY ELEMENT
C  IN A LINE.
C
C
   20 IF (I .GT. NITEMS)  GO TO 90
        J = J+1
C/R
C       LINE(1,J) = A(1,I)
C       LINE(2,J) = A(2,I)
C/C
        LINE(J) = A(I)
C/
        IF (J .EQ. 1) ILINE = I
        IF (J .LT. NCOL .AND. I .LT. NITEMS) GO TO 80
          IF (COUNT .EQ. 0) GO TO 50
            DUP = .TRUE.
            DO 30 K=1,NCOL
C/R
C             IF (LAST(1,K) .NE. LINE(1,K)  .OR.
C    1            LAST(2,K) .NE. LINE(2,K))
C    2            DUP = .FALSE.
C/C
              IF (REAL(LAST(K)) .NE. REAL(LINE(K))  .OR.
     1            AIMAG(LAST(K)) .NE. AIMAG(LINE(K)))
     2            DUP = .FALSE.
C/
   30       CONTINUE
            IF (I .EQ. NITEMS  .AND.  J .LT. NCOL) DUP = .FALSE.
            IF (.NOT. DUP .AND. COUNT .EQ. 1) GO TO 50
              IF (.NOT. DUP) GO TO 40
                COUNT = COUNT+1
                IF (COUNT .EQ. 3) WRITE(IOUT, IFMT1C) BLANK,
     1                                 STAR, STAR, STAR, STAR
                IF (I .EQ. NITEMS)  GO TO 50
                  GO TO 70
C/R
C  40         WRITE(IOUT, IFMT2C) BLANK, ILAST, (LAST(1,K),
C    1              LAST(2,K), K=1,NCOL)
C  50     WRITE(IOUT, IFMT2C) BLANK, ILINE, (LINE(1,K),
C    1              LINE(2,K), K=1,J)
C/C
   40         WRITE(IOUT, IFMT2C) BLANK, ILAST, (LAST(K), K=1,NCOL)
C*******COMMENTED OUT BY CHARLES THOMPSON**************
C reinstated for alliant fx8
  50     WRITE(IOUT, IFMT2C) BLANK, ILINE, (LINE(K), K=1,J)
C   50     WRITE(IOUT,*) 'BUG IN 50 A9RNTC'
C/
          COUNT = 1
          DO 60 K=1,NCOL
C/R
C           LAST(1,K) = LINE(1,K)
C  60       LAST(2,K) = LINE(2,K)
C/C
   60       LAST(K) = LINE(K)
C/
   70     ILAST = ILINE
          J = 0
   80   I = I+1
        GO TO 20
   90 RETURN
      END
