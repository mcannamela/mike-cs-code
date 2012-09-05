*DECK PFIND
      SUBROUTINE PFIND
C===============================================================PFIND
C      THIS SUBROUTINE DETERMINES THE CELL (ICP(IP)) THAT CONTAINS
C      THE PARCEL WHOSE INDEX IS IP.
C----------
C     CALLED BY LAVA
C===============================================================PFIND
      INCLUDE 'COML.h'
C----------
      DOUBLE PRECISION RRP
      INTEGER ICPSAV
C===============================================================PFIND
C     FOR 2-D AXISYMMETRIC WITH 3-D PARTCILE (PSEUDO 3-D)
C---------------------------------------------------------------PFIND
      IF(I3DP2D.EQ.1) THEN
        DO 100 IP=1,NP
        ICPSAV=ICP(IP)
        J=ICPSAV/NXT+1
        I=ICPSAV-(J-1)*NXT
        RRP=SQRT(XP(IP)**2+ZP(IP)**2)
C----------
  110   IF(RRP.GE.X(I-1) .AND. RRP.LT.X(I)) GO TO 120
          IF(RRP.LE.X(I-1)) I=I-1
          IF(RRP.GT.X(I)) I=I+1
          GO TO 110
C----------
  120   IF(YP(IP).GE.Y(J-1) .AND. YP(IP).LT.Y(J)) GO TO 130
          IF(YP(IP).LE.Y(J-1)) J=J-1
          IF(YP(IP).GT.Y(J)) J=J+1
          GO TO 120
C----------
  130   ICP(IP)=(J-1)*NXT+I
  100   CONTINUE
C---------------------------------------------------------------PFIND
C    FOR PLANAR 2-D, 2-D AXISYMMETRIC WITH 2-D PARTCILE, AND FULL 3-D
C---------------------------------------------------------------PFIND
      ELSE
        DO 200 IP=1,NP
        ICPSAV=ICP(IP)
        K=ICPSAV/NXYT+1
        J=MOD(ICPSAV,NXYT)/NXT+1
        I=ICPSAV-(J-1)*NXT-(K-1)*NXYT
C----------
  210   IF(XP(IP).GE.X(I-1) .AND. XP(IP).LT.X(I)) GO TO 220
          IF(XP(IP).LE.X(I-1)) I=I-1
          IF(XP(IP).GT.X(I)) I=I+1
          GO TO 210
C----------
  220   IF(YP(IP).GE.Y(J-1) .AND. YP(IP).LT.Y(J)) GO TO 230
          IF(YP(IP).LE.Y(J-1)) J=J-1
          IF(YP(IP).GT.Y(J)) J=J+1
          GO TO 220
C----------
  230   IF(ZP(IP).GE.Z(K-1) .AND. ZP(IP).LT.Z(K)) GO TO 240
          IF(ZP(IP).LE.Z(K-1)) K=K-1
          IF(ZP(IP).GT.Z(K)) K=K+1
          GO TO 230
C----------
  240   ICP(IP)=(K-1)*NXYT+(J-1)*NXT+I
  200   CONTINUE
      END IF
C===============================================================PFIND
      RETURN
      END
