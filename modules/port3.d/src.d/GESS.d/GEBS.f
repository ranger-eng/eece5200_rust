       SUBROUTINE GEBS(N,A,IA,B,IB,NB)
CMNEMONIC-GENERAL BACK SOLVE
C THIS SUBROUTINE SOLVES AX=B WHERE A IS
C AN UPPER TRIANGULAR MATRIX
C INPUT PARAMETERS
C    N    ORDER OF THE PROBLEM
C    A    AN IA X N ARRAY CONTAINING THE MATRIX A
C    IA   ROW DIMENSION OF THE A MATRIX
C    B    AN IB X NB MATRIX CONTAINING THE RIGHT HAND SIDES
C         DESTROYED ON OUTPUT
C    IB   ROW DIMENSION OF THE B ARRAY, MUST BE AT LEAST N
C    NB   NUMBER OF RIGHT HAND SIDES
C OUTPUT PARAMTERS
C   B    THE SOLUTION X
C THIS SUBROUTINE USES SDOT
C ERROR CONDITIONS
C  1      N LESS THAN 1     FATAL
C  2      IA .LT.N          FATAL
C  3      IB.LT.N           FATAL
C  4      NB.LT.1           FATAL
C  10+K   SINGULAR MATRIX OF RANK K      RECOVERABLE
C EXTRA STORAGE ALLOCATED-NONE
       INTEGER N,IA,IB,NB
       REAL A(IA,N),B(IB,NB),T
C/6S
       IF (N.LT.1) CALL SETERR(12H GEBS-N.LT.1,12,1,2)
       IF (IA.LT.N) CALL SETERR(13H GEBS-IA.LT.N,13,2,2)
       IF (IB.LT.N) CALL SETERR(13H GEBS-IB.LT.N,13,3,2)
       IF (NB.LT.1) CALL SETERR(13H GEBS-NB.LT.1,13,4,2)
C/7S
C      IF (N.LT.1) CALL SETERR(' GEBS-N.LT.1',12,1,2)
C      IF (IA.LT.N) CALL SETERR(' GEBS-IA.LT.N',13,2,2)
C      IF (IB.LT.N) CALL SETERR(' GEBS-IB.LT.N',13,3,2)
C      IF (NB.LT.1) CALL SETERR(' GEBS-NB.LT.1',13,4,2)
C/
        CALL ENTER(1)
C DO THE BACK SOLVE
 50    NP1=N+1
       DO 70 KB=1,N
          K=NP1-KB
          DO 65 I=1,NB
          T=B(K,I)
          IF(K.NE.N)T=T-SDOT(KB-1,A(K,K+1),IA,B(K+1,I),1)
          IF(A(K,K).NE.0.0) GO TO 60
C/6S
          CALL SETERR(38H GEBS-DIVISION BY ZERO-SINGULAR MATRIX,
     1    38,9+K,1)
C/7S
C         CALL SETERR(' GEBS-DIVISION BY ZERO-SINGULAR MATRIX',
C    1    38,9+K,1)
C/
          GO TO 75
 60       B(K,I)=T/A(K,K)
 65       CONTINUE
 70     CONTINUE
 75     CALL LEAVE
        RETURN
        END
