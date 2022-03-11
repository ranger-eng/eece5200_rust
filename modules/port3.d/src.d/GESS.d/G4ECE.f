      SUBROUTINE G4ECE(N,A,IA,ANORM,COND,IPVT,Z)
      INTEGER IA,N,IPVT(N)
      REAL A(IA, N), ANORM, COND, Z(N)
      INTEGER KB,  KP1,  J, K
      INTEGER L
      REAL EK, GREAT, SM, WK, SIGN, SDOT, ABS
      REAL S, T, WKM, SASUM, YNORM
      REAL R1MACH,BIG
C THIS SUBROUTINE DETERMINES A LOWER BOUND ON THE CONDITION NUMBER
C OF THE DECOMPOSED MATRIX A VIA THE ALGORITHM USED IN LINPACK
C
C
C SOLVE A(TRANSPOSE)W = E
C WHERE E IS CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH
C IN THE COMPONENTS OF W
      IF (N.NE.1) GO TO 1
          COND=1.0
           Z(1)=1.0
          RETURN
 1    EK = 1.0
      BIG=SQRT(R1MACH(2))/FLOAT(N)
      DO  2 J = 1, N
         Z(J) = 0.0
  2   CONTINUE
      DO 30 K=1,N
        IF (ABS(Z(K)) .NE. 0.0) EK=SIGN(EK,-Z(K))
        IF (ABS(EK-Z(K)) .LE. ABS(A(K,K))) GO TO 20
           S=ABS(A(K,K))/ABS(EK-Z(K))
           CALL SSCAL(N,S,Z,1)
           EK=S*EK
 20     CONTINUE
        WK=EK - Z(K)
        WKM=-EK-Z(K)
        S = ABS(WK)
        SM = ABS(WKM)
        IF (A(K,K).EQ.0.0) GO TO 22
        WK = WK/A(K,K)
        WKM = WKM/A(K,K)
        GO TO 23
 22     WK=1.0
        WKM=1.0
 23     CONTINUE
        KP1 = K +1
        IF (KP1. GT. N) GO TO 28
        DO 24 J=KP1,N
          SM = SM + ABS(Z(J) + WKM*A(K,J))
          Z(J) = Z(J) + WK * A(K,J)
          S = S + ABS(Z(J))
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
       S= 1.0/SASUM(N,Z,1)
       CALL SSCAL(N,S,Z,1)
C
C FORM Y = L(TRANSPOSE) * W
C
      DO  40 KB = 2, N
         K = N+1-KB
         IF (K .LT. N) Z(K) = Z(K)+SDOT(N-K, A(K+1, K), 1, Z(K+1), 1)
         IF (ABS(Z(K)).LT.BIG) GO TO 35
             S=1.0/ABS(Z(K))
             CALL SSCAL(N,S,Z,1)
 35      CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  40  CONTINUE
      S=1.0/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)
      YNORM = 1.0
C
C   FORM W = L * Y
C
       IF (N.EQ.1)GO TO 46
       NM1=N-1
       DO 45 K=1,NM1
          L=IPVT(K)
          T=Z(L)
          Z(L)=Z(K)
          Z(K)=T
          CALL SAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)
          IF(ABS(Z(K)).LE.BIG) GO TO 45
            S=1.0/ABS(Z(K))
            CALL SSCAL(N,S,Z,1)
            YNORM=S*YNORM
 45     CONTINUE
       S = 1.0/SASUM(N,Z,1)
       IF (S.GT.1.0) GO TO 46
       CALL SSCAL(N,S,Z,1)
       YNORM = YNORM*S
C
C   SOLVE U * Z = W
 46    DO  50 KB = 1, N
         K = N+1-KB
         IF (ABS(Z(K)) .LE. ABS(A(K,K))) GO TO 48
            S=ABS(A(K,K))/ABS(Z(K))
            CALL SSCAL(N,S,Z,1)
            YNORM = YNORM * S
 48      CONTINUE
         IF(A(K,K).NE.0.0)Z(K) = Z(K)/A(K,K)
         IF (A(K,K).EQ.0.0)Z(K)=1.0
         T = -Z(K)
         IF  (K .NE. 1) CALL SAXPY(K-1,T,A(1,K),1,Z(1),1)
 50      CONTINUE
C    MAKE ZNORM = 1.0
       S= 1.0/SASUM(N,Z,1)
       CALL SSCAL(N,S,Z,1)
       YNORM = YNORM*S
C
C   SET COND = ESTIMATE OF THE CONDITION NUMBER OF A
C
       GREAT=R1MACH(2)
       IF (YNORM.GT.1.0) GO TO 60
       IF (ANORM.LE.YNORM*GREAT) GO TO 60
       COND=GREAT
       RETURN
 60    COND=ANORM/YNORM
       RETURN
       END
