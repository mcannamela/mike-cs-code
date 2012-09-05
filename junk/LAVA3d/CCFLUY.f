*DECK CCFLUY
      SUBROUTINE CCFLUY
C==============================================================CCFLUY
C     THIS ROUTINE CALCULATES CONVECTIVE FLUXING OF CELL-CENTERED
C     QUANTITIES IN Y-DIRECTION
C----------
C     CALLED BY LAVA
C----------
C     INPUT=SPDN,SPDR,RON,ROER,RORW,XC,DY,RDY,RF,DT,UN,VN,WN
C     OUTPUT=SPDR,ROER,RORW,
C            ROVU,ROVV,ROVW -- AT FRONT SIDE OF THE CELL
C     PARAMETER=ADC,BDC
C==============================================================CCFLUY
      INCLUDE 'COML.h'
C----------
      INTEGER IF
      DOUBLE PRECISION RDYCF,COUR,SS,DELF,DENC,DENF
C
C==============================================================CCFLUY
      IF(INCOMP.EQ.1) GO TO 290
C--------------------------------------------------------------CCFLUY
      DO 10 IC=IC1V,IC2
      RDYCF=TWO*HRDYCF(IC)
      DELF=DY(IC+NXT)*RDYCF
      COUR=DT(IC)*VN(IC)*RDYCF
      SS=SIGN(ONE,COUR)
C---------- INTERPOLATED DONOR-CELL
      WT(IC)=HALF*(DELF+ADC*(ONE-DELF+SS)+BDC*COUR)
   10 CONTINUE
C=============================== MASS AND SPECIES DENSITY =====CCFLUY
      DO 100 ISP=1,NSP
        DO 20 IC=IC1V,IC2
        DENC=SPDN(IC,ISP)
        DENF=SPDN(IC+NXT,ISP)
        FL(IC)=(DENC*WT(IC)+DENF*(ONE-WT(IC)))*RF(IC)*VN(IC)
   20   CONTINUE
        DO 30 IC=IC1,IC2
        SPDR(IC,ISP)=SPDR(IC,ISP)-DT(IC)*(FL(IC)-FL(IC-NXT))*RDY(IC)
   30   CONTINUE
  100 CONTINUE
C========================================== UPDATE ENERGY =====CCFLUY
      DO 210 IC=IC1V,IC2
      DENC=ROEN(IC)
      DENF=ROEN(IC+NXT)
      FL(IC)=(DENC*WT(IC)+DENF*(ONE-WT(IC)))*RF(IC)*VN(IC)
  210 CONTINUE
      DO 220 IC=IC1,IC2
      ROER(IC)=ROER(IC)-DT(IC)*(FL(IC)-FL(IC-NXT))*RDY(IC)
  220 CONTINUE
C================================= UPDATE ELECTRON ENERGY =====CCFLUY
      IF(NLTE.NE.0) THEN
        DO 260 IC=IC1V,IC2
        DENC=ROEEN(IC)
        DENF=ROEEN(IC+NXT)
        FL(IC)=(DENC*WT(IC)+DENF*(ONE-WT(IC)))*RF(IC)*VN(IC)
  260   CONTINUE
        DO 270 IC=IC1,IC2
        ROEER(IC)=ROEER(IC)-DT(IC)*(FL(IC)-FL(IC-NXT))*RDY(IC)
  270   CONTINUE
      END IF
  290 CONTINUE
C======================== UPDATE SWIRL VELOCITY (2D ONLY) =====CCFLUY
      IF(ISWIRL.EQ.1) THEN
        DO 310 IC=IC1V,IC2
        IF=IC+NXT
        DENC=RON(IC)*WN(IC)
        DENF=RON(IF)*WN(IF)
        FL(IC)=(DENC*WT(IC)+DENF*(ONE-WT(IC)))*RF(IC)*VN(IC)
  310   CONTINUE
        DO 320 IC=IC1,IC2
        RORW(IC)=RORW(IC)-DT(IC)*(FL(IC)-FL(IC-NXT))*RDY(IC)
  320   CONTINUE
      END IF
C========================== FOR TURBULENCT KINETIC ENERGY =====CCFLUY
      IF(ITURB.NE.0) THEN
        DO 410 IC=IC1V,IC2
        DENC=RON(IC)*TKEN(IC)
        DENF=RON(IC+NXT)*TKEN(IC+NXT)
        FL(IC)=(DENC*WT(IC)+DENF*(ONE-WT(IC)))*RF(IC)*VN(IC)
  410   CONTINUE
        DO 420 IC=IC1,IC2
        ROTKER(IC)=ROTKER(IC)-DT(IC)*(FL(IC)-FL(IC-NXT))*RDY(IC)
  420   CONTINUE
      END IF
C----------- FOR DISSIPATION OF TURBULENCT KINETIC ENERGY -----CCFLUY
      IF(ITURB.EQ.2) THEN
        DO 430 IC=IC1V,IC2
        DENC=RON(IC)*EPSN(IC)
        DENF=RON(IC+NXT)*EPSN(IC+NXT)
        FL(IC)=(DENC*WT(IC)+DENF*(ONE-WT(IC)))*RF(IC)*VN(IC)
  430   CONTINUE
        DO 440 IC=IC1,IC2
        ROEPSR(IC)=ROEPSR(IC)-DT(IC)*(FL(IC)-FL(IC-NXT))*RDY(IC)
  440   CONTINUE
      END IF
C==============================================================CCFLUY
      RETURN
      END
