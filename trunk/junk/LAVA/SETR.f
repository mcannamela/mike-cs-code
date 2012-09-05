*DECK SETR
      SUBROUTINE SETR
C================================================================SETR
C     THIS ROUTINE PROVIDES "R" FACTOR
C     R = (VOID FRACTION)*(RADIAL DIRECTION IN 2-D CYLINDRICAL COORD.)
C----------
C     CALLED BY LAVA
C----------
C     INPUT: NP, PARTICLE STUFF, XC, DX
C     OUTPUT: RC, RR, RF, RT, RRC
C================================================================SETR
      INCLUDE 'COML.h'
C----------
      DOUBLE PRECISION AXC,AXR,VOIDC,VOIDR,VOIDF,VOIDT
C
C================================================================SETR
      DO 100 K=1,NZT
      ICK=(K-1)*NXYT
      DO 100 J=1,NYT
      ICJ=ICK+(J-1)*NXT
      DO 100 I=1,NXT
      IC=I+ICJ
C----------------------------------------------------------------SETR
      AXC=ONE
      AXR=ONE
      IF(ICYL.EQ.1) THEN
        AXC=ABS(XC(IC))
        AXR=ABS(XC(IC)+HALF*DX(IC))
      END IF
C----------------------------------------------------------------SETR
      VOIDC=ONE
      VOIDR=ONE
      VOIDF=ONE
      VOIDT=ONE
C----------------------------------------------------------------SETR
      RC(IC)=VOIDC*AXC
      RR(IC)=VOIDR*AXR
      RF(IC)=VOIDF*AXC
      RT(IC)=VOIDT*AXC
      RRC(IC)=ONE
      IF(RC(IC).GT.SMALL) RRC(IC)=ONE/RC(IC)
  100 CONTINUE
C================================================================SETR
      RETURN
      END
