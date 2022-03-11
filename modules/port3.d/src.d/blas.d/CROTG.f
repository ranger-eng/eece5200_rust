      SUBROUTINE CROTG(CA,CB,C,S)
C
C     ADAPTED TO PORT 3 BY PHYL FOX, 11/8/83
C
      COMPLEX CA,CB,S
      REAL C
      REAL NORM,SCALE
      COMPLEX ALPHA
C
      IF (CABS(CA) .NE. 0.) GO TO 10
      IF (CABS(CB) .EQ. 0.) GO TO  5
         C = 0.
         S = (1.,0.)
         CA = CB
         GO TO 20
    5    C = 1.0
         S = (0.0,0.0)
         GO TO 20
   10 CONTINUE
         SCALE = CABS(CA) + CABS(CB)
         NORM = SCALE * SQRT((CABS(CA/SCALE))**2 + (CABS(CB/SCALE))**2)
         ALPHA = CA /CABS(CA)
         C = CABS(CA) / NORM
         S = ALPHA * CONJG(CB) / NORM
         CA = ALPHA * NORM
   20 CONTINUE
         CB = (0.0,0.0)
      RETURN
      END
