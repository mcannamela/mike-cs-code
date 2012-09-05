*DECK CARRY
      SUBROUTINE CARRY
C===============================================================CARRY
C     THIS ROUTINE PROVIDES MASS, MOMENTUM, AND ENERGY SOURCES
C     DUE TO GASES INJECTED THROUGH NOZZLES.
C     THOSE SOURCES ARE PROVIDED BY THE SUBROUTINE RINPUT.
C----------
C     CALLED BY LAVA
C===============================================================CARRY
      INCLUDE 'COML.h'
C---------------------------------------------------------------CARRY
      INTEGER NOZ
C============================================= MASS SOURCE =====CARRY
      DO 200 NOZ=1,NUMNOZ
      IC=ICNOZ(NOZ)
      WT(IC)=DT(IC)*RC(IC)/VOL(IC)
        DO 100 ISP=1,NSP
        SPDR(IC,ISP)=SPDR(IC,ISP)+CARRYG(NOZ,ISP)*WT(IC)
  100   CONTINUE
  200 CONTINUE
C========================================= MOMENTUM SOURCE =====CARRY
C     MOMENTUM SOURCES ARE EQUALLY DIVIDED TO EACH OF MOMENTUM CELLS
C---------------------------------------------------------------CARRY
      DO 300 NOZ=1,NUMNOZ
      IC=ICNOZ(NOZ)
      RORU(IC)=RORU(IC)+DT(IC)*RR(IC)*CRMOMX(NOZ)/(P(IC)*VOL(IC))
      RORU(IC-1)=RORU(IC-1)
     1           +DT(IC-1)*RR(IC-1)*CRMOMX(NOZ)/(P(IC)*VOL(IC-1))
      RORV(IC)=RORV(IC)+DT(IC)*RF(IC)*CRMOMY(NOZ)/(P(IC)*VOL(IC))
      RORV(IC-NXT)=RORV(IC-NXT)
     1         +DT(IC-NXT)*RF(IC-NXT)*CRMOMY(NOZ)/(P(IC)*VOL(IC-NXT))
  300 CONTINUE
      IF(NZ.GT.1) THEN
      DO 400 NOZ=1,NUMNOZ
      IC=ICNOZ(NOZ)
      RORW(IC)=RORW(IC)+DT(IC)*RT(IC)*CRMOMZ(NOZ)/(P(IC)*VOL(IC))
      RORW(IC-NXYT)=RORW(IC-NXYT)
     1      +DT(IC-NXYT)*RT(IC-NXYT)*CRMOMZ(NOZ)/(P(IC)*VOL(IC-NXYT))
  400 CONTINUE
      END IF
C============================================ SWIRL SOURCE =====CARRY
C     SWIRL SOURCE IS AT THE CELL CENTER, AND IS TORQUE/RADIUS
C     FOR RORW
C---------------------------------------------------------------CARRY
      IF(ISWIRL.EQ.1) THEN
      DO 500 NOZ=1,NUMNOZ
      IC=ICNOZ(NOZ)
      RORW(IC)=RORW(IC)+DT(IC)*RC(IC)*TWO*CRMOMZ(NOZ)/VOL(IC)
  500 CONTINUE
      END IF
C=========================================== ENERGY SOURCE =====CARRY
      DO 600 NOZ=1,NUMNOZ
      IC=ICNOZ(NOZ)
      ROER(IC)=ROER(IC)+DT(IC)*RC(IC)*CRENGY(NOZ)/VOL(IC)
  600 CONTINUE
C===============================================================CARRY
      RETURN
      END
