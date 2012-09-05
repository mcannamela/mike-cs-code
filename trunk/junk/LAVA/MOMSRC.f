*DECK MOMSRC
      SUBROUTINE MOMSRC
C==============================================================MOMSRC
C
C
C
C----------
C     CALLED BY LAVA
C==============================================================MOMSRC
      INCLUDE 'COML.h'
C----------
      DOUBLE PRECISION SWRL
C==============================================================MOMSRC
      DO 10 IC=IC1,IC2
      U(IC)=U(IC)+DT(IC)*GX
      V(IC)=V(IC)+DT(IC)*GY
      W(IC)=W(IC)+DT(IC)*GZ
   10 CONTINUE
C==============================================================MOMSRC
      IF(ISWIRL.EQ.1) THEN
        DO 20 IC=IC1U,IC2
        SWRL=(DX(IC)*WN(IC+1)+DX(IC+1)*WN(IC))/(DX(IC)+DX(IC+1))
        U(IC)=U(IC)+DT(IC)*SWRL*SWRL/MAX((XC(IC)+HALF*DX(IC)),SMALL)
   20   CONTINUE
      END IF
C==============================================================MOMSRC
      RETURN
      END
