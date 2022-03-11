      SUBROUTINE DG4ECE(N,A,IA,ANORM,COND,IPVT,Z)
      INTEGER IA,N,IPVT(N)
      DOUBLE PRECISION A(IA, N), ANORM, COND, Z(N)
      INTEGER KB,  KP1,  J, K
      INTEGER L
      DOUBLE PRECISION EK, GREAT, SM, WK, DDOT
      DOUBLE PRECISION S, T, WKM, DASUM, YNORM
      DOUBLE PRECISION D1MACH,BIG,DSQRT
C THIS SUBROUTINE DETERMINES A LOWER BOUND ON THE CONDITION NUMBER
C OF THE DECOMPOSED MATRIX A VIA THE ALGORITHM USED IN LINPACK
C
C
C SOLVE A(TRANSPOSE)W = E
C WHERE E IS CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH
C IN THE COMPONENTS OF W
      IF (N.NE.1) GO TO 1
          COND=1.D0
           Z(1)=1.D0
          RETURN
 1    EK = 1.D0
      BIG=DSQRT(D1MACH(2))/DBLE(FLOAT(N))
      DO  2 J = 1, N
         Z(J) = 0.D0
  2   CONTINUE
      DO 30 K=1,N
        IF (DABS(Z(K)) .NE. 0.D0) EK=DSIGN(EK,-Z(K))
        IF (DABS(EK-Z(K)) .LE. DABS(A(K,K))) GO TO 20
           S=DABS(A(K,K))/DABS(EK-Z(K))
           CALL DSCAL(N,S,Z,1)
           EK=S*EK
 20     CONTINUE
        WK=EK - Z(K)
        WKM=-EK-Z(K)
        S = DABS(WK)
        SM = DABS(WKM)
        IF (A(K,K).EQ.0.D0) GO TO 22
        WK = WK/A(K,K)
        WKM = WKM/A(K,K)
        GO TO 23
 22     WK=1.D0
        WKM=1.D0
 23     CONTINUE
        KP1 = K +1
        IF (KP1. GT. N) GO TO 28
        DO 24 J=KP1,N
          SM = SM + DABS(Z(J) + WKM*A(K,J))
          Z(J) = Z(J) + WK * A(K,J)
          S = S + DABS(Z(J))
 24    CONTINUE
       IF ( S .GE. SM) GO TO 28
          T= WKM - WK
          WK = WKM
          DO 26 J=KP1,N
             Z(J) = Z(J) + T*A(K,J)
 26       CONTINUE
 28    CONTINUE
       Z(K) = WK
 30    CONTINUE
       S= 1.D0/DASUM(N,Z,1)
       CALL DSCAL(N,S,Z,1)
C
C FORM Y = L(TRANSPOSE) * W
C
      DO  40 KB = 2, N
         K = N+1-KB
         IF (K .LT. N) Z(K) = Z(K)+DDOT(N-K, A(K+1, K), 1, Z(K+1), 1)
         IF (DABS(Z(K)).LT.BIG) GO TO 35
             S=1.D0/DABS(Z(K))
             CALL DSCAL(N,S,Z,1)
 35      CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  40  CONTINUE
      S=1.D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = 1.D0
C
C   FORM W = L * Y
C
       IF (N.EQ.1) GO TO46
       NM1=N-1
       DO 45 K=1,NM1
          L=IPVT(K)
          T=Z(L)
          Z(L)=Z(K)
          Z(K)=T
          CALL DAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)
          IF(DABS(Z(K)).LE.BIG) GO TO 45
            S=1.D0/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM=S*YNORM
 45     CONTINUE
       S = 1.D0/DASUM(N,Z,1)
       IF (S.GT.1.D0) GO TO 46
       CALL DSCAL(N,S,Z,1)
       YNORM = YNORM*S
C
C   SOLVE U * Z = W
 46    DO  50 KB = 1, N
         K = N+1-KB
         IF (DABS(Z(K)) .LE. DABS(A(K,K))) GO TO 48
            S=DABS(A(K,K))/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = YNORM * S
 48      CONTINUE
         IF(A(K,K).NE.0.D0)Z(K) = Z(K)/A(K,K)
         IF (A(K,K).EQ.0.D0)Z(K)=1.D0
         T = -Z(K)
         IF  (K .NE. 1) CALL DAXPY(K-1,T,A(1,K),1,Z(1),1)
 50      CONTINUE
C    MAKE ZNORM = 1.D0
       S= 1.D0/DASUM(N,Z,1)
       CALL DSCAL(N,S,Z,1)
       YNORM = YNORM*S
C
C   SET COND = ESTIMATE OF THE CONDITION NUMBER OF A
C
       GREAT=D1MACH(2)
       IF (YNORM.GT.1.D0) GO TO 60
       IF (ANORM.LE.YNORM*GREAT) GO TO 60
       COND=GREAT
       RETURN
 60    COND=ANORM/YNORM
       RETURN
       END
