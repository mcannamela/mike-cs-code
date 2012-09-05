*DECK CCFLUX
      SUBROUTINE CCFLUX
C==============================================================CCFLUX
C     THIS ROUTINE CALCULATES CONVECTIVE FLUXING OF CELL-CENTERED
C     QUANTITIES IN X-DIRECTION
C----------
C     CALLED BY LAVA
C----------
C     INPUT=SPDN,SPDR,RON,ROER,RORW,XC,DX,RDX,RR,DT,UN,VN,WN
C     OUTPUT=SPDR,ROER,RORW,
C            ROUU,ROUV,ROUW -- AT RIGHT SIDE OF THE CELL
C     PARAMETER=ADC,BDC
C==============================================================CCFLUX
      INCLUDE 'COML.h'
C----------
      INTEGER IR
      DOUBLE PRECISION RDXCR,COUR,SS,DELR,DENC,DENR
C
C==============================================================CCFLUX
      IF(INCOMP.EQ.1) GO TO 290
C--------------------------------------------------------------CCFLUX
      DO 10 IC=IC1U,IC2
      RDXCR=TWO*HRDXCR(IC)
      DELR=DX(IC+1)*RDXCR
      COUR=DT(IC)*UN(IC)*RDXCR
      SS=SIGN(ONE,COUR)
C---------- INTERPOLATED DONOR-CELL
      WT(IC)=HALF*(DELR+ADC*(ONE-DELR+SS)+BDC*COUR)
   10 CONTINUE
C=============================== MASS AND SPECIES DENSITY =====CCFLUX
      DO 100 ISP=1,NSP
        DO 20 IC=IC1U,IC2
        DENC=SPDN(IC,ISP)
        DENR=SPDN(IC+1,ISP)
        FL(IC)=(DENC*WT(IC)+DENR*(ONE-WT(IC)))*RR(IC)*UN(IC)
   20   CONTINUE
        DO 30 IC=IC1,IC2
        SPDR(IC,ISP)=SPDR(IC,ISP)-DT(IC)*(FL(IC)-FL(IC-1))*RDX(IC)
   30   CONTINUE
  100 CONTINUE
C========================================== UPDATE ENERGY =====CCFLUX
      DO 210 IC=IC1U,IC2
      DENC=ROEN(IC)
      DENR=ROEN(IC+1)
      FL(IC)=(DENC*WT(IC)+DENR*(ONE-WT(IC)))*RR(IC)*UN(IC)
  210 CONTINUE
      DO 220 IC=IC1,IC2
      ROER(IC)=ROER(IC)-DT(IC)*(FL(IC)-FL(IC-1))*RDX(IC)
  220 CONTINUE
C================================= UPDATE ELECTRON ENERGY =====CCFLUX
      IF(NLTE.NE.0) THEN
        DO 260 IC=IC1U,IC2
        DENC=ROEEN(IC)
        DENR=ROEEN(IC+1)
        FL(IC)=(DENC*WT(IC)+DENR*(ONE-WT(IC)))*RR(IC)*UN(IC)
  260   CONTINUE
        DO 270 IC=IC1,IC2
        ROEER(IC)=ROEER(IC)-DT(IC)*(FL(IC)-FL(IC-1))*RDX(IC)
  270   CONTINUE
      END IF
  290 CONTINUE
C======================== UPDATE SWIRL VELOCITY (2D ONLY) =====CCFLUX
      IF(ISWIRL.EQ.1) THEN
        DO 310 IC=IC1U,IC2
        IR=IC+1
        DENC=RON(IC)*XC(IC)*WN(IC)
        DENR=RON(IR)*XC(IR)*WN(IR)
        FL(IC)=(DENC*WT(IC)+DENR*(ONE-WT(IC)))*RR(IC)*UN(IC)
  310   CONTINUE
        DO 320 IC=IC1,IC2
        RORW(IC)=RORW(IC)-DT(IC)*(FL(IC)-FL(IC-1))*RDX(IC)/XC(IC)
  320   CONTINUE
      END IF
C========================== FOR TURBULENCT KINETIC ENERGY =====CCFLUX
      IF(ITURB.NE.0) THEN
        DO 410 IC=IC1U,IC2
        DENC=RON(IC)*TKEN(IC)
        DENR=RON(IC+1)*TKEN(IC+1)
        FL(IC)=(DENC*WT(IC)+DENR*(ONE-WT(IC)))*RR(IC)*UN(IC)
  410   CONTINUE
        DO 420 IC=IC1,IC2
        ROTKER(IC)=ROTKER(IC)-DT(IC)*(FL(IC)-FL(IC-1))*RDX(IC)
  420   CONTINUE
      END IF
C----------- FOR DISSIPATION OF TURBULENCT KINETIC ENERGY -----CCFLUX
      IF(ITURB.EQ.2) THEN
        DO 430 IC=IC1U,IC2
        DENC=RON(IC)*EPSN(IC)
        DENR=RON(IC+1)*EPSN(IC+1)
        FL(IC)=(DENC*WT(IC)+DENR*(ONE-WT(IC)))*RR(IC)*UN(IC)
  430   CONTINUE
        DO 440 IC=IC1,IC2
        ROEPSR(IC)=ROEPSR(IC)-DT(IC)*(FL(IC)-FL(IC-1))*RDX(IC)
  440   CONTINUE
      END IF
C==============================================================CCFLUX
      RETURN
      END
