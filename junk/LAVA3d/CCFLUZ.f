*DECK CCFLUZ
      SUBROUTINE CCFLUZ
C==============================================================CCFLUZ
C     THIS ROUTINE CALCULATES CONVECTIVE FLUXING OF CELL-CENTERED
C     QUANTITIES IN Z-DIRECTION
C----------
C     CALLED BY LAVA
C----------
C     INPUT=SPDN,SPDR,RON,ROER,RORW,DZ,RDZ,RT,DT,UN,VN,WN
C     OUTPUT=SPDR,ROER,RORW,
C            ROWU,ROWV,ROWW -- AT TOP SIDE OF THE CELL
C     PARAMETER=ADC,BDC
C==============================================================CCFLUZ
      INCLUDE 'COML.h'
C----------
      DOUBLE PRECISION RDZCT,COUR,SS,DELT,DENC,DENT
C
C==============================================================CCFLUZ
      IF(INCOMP.EQ.1) GO TO 290
C--------------------------------------------------------------CCFLUZ
      DO 10 IC=IC1W,IC2
      RDZCT=TWO*HRDZCT(IC)
      DELT=DZ(IC+NXYT)*RDZCT
      COUR=DT(IC)*WN(IC)*RDZCT
      SS=SIGN(ONE,COUR)
C---------- INTERPOLATED DONOR-CELL
      WT(IC)=HALF*(DELT+ADC*(ONE-DELT+SS)+BDC*COUR)
   10 CONTINUE
C=============================== MASS AND SPECIES DENSITY =====CCFLUZ
      DO 100 ISP=1,NSP
        DO 20 IC=IC1W,IC2
        DENC=SPDN(IC,ISP)
        DENT=SPDN(IC+NXYT,ISP)
        FL(IC)=(DENC*WT(IC)+DENT*(ONE-WT(IC)))*RT(IC)*WN(IC)
   20   CONTINUE
        DO 30 IC=IC1,IC2
        SPDR(IC,ISP)=SPDR(IC,ISP)-DT(IC)*(FL(IC)-FL(IC-NXYT))*RDZ(IC)
   30   CONTINUE
  100 CONTINUE
C========================================== UPDATE ENERGY =====CCFLUZ
      DO 210 IC=IC1W,IC2
      DENC=ROEN(IC)
      DENT=ROEN(IC+NXYT)
      FL(IC)=(DENC*WT(IC)+DENT*(ONE-WT(IC)))*RT(IC)*WN(IC)
  210 CONTINUE
      DO 220 IC=IC1,IC2
      ROER(IC)=ROER(IC)-DT(IC)*(FL(IC)-FL(IC-NXYT))*RDZ(IC)
  220 CONTINUE
C================================= UPDATE ELECTRON ENERGY =====CCFLUZ
      IF(NLTE.NE.0) THEN
        DO 260 IC=IC1W,IC2
        DENC=ROEEN(IC)
        DENT=ROEEN(IC+NXYT)
        FL(IC)=(DENC*WT(IC)+DENT*(ONE-WT(IC)))*RT(IC)*WN(IC)
  260   CONTINUE
        DO 270 IC=IC1,IC2
        ROEER(IC)=ROEER(IC)-DT(IC)*(FL(IC)-FL(IC-NXYT))*RDZ(IC)
  270   CONTINUE
      END IF
  290 CONTINUE
C========================== FOR TURBULENCT KINETIC ENERGY =====CCFLUZ
      IF(ITURB.NE.0) THEN
        DO 410 IC=IC1W,IC2
        DENC=RON(IC)*TKEN(IC)
        DENT=RON(IC+NXYT)*TKEN(IC+NXYT)
        FL(IC)=(DENC*WT(IC)+DENT*(ONE-WT(IC)))*RT(IC)*WN(IC)
  410   CONTINUE
        DO 420 IC=IC1,IC2
        ROTKER(IC)=ROTKER(IC)-DT(IC)*(FL(IC)-FL(IC-NXYT))*RDZ(IC)
  420   CONTINUE
      END IF
C----------- FOR DISSIPATION OF TURBULENCT KINETIC ENERGY -----CCFLUZ
      IF(ITURB.EQ.2) THEN
        DO 430 IC=IC1W,IC2
        DENC=RON(IC)*EPSN(IC)
        DENT=RON(IC+NXYT)*EPSN(IC+NXYT)
        FL(IC)=(DENC*WT(IC)+DENT*(ONE-WT(IC)))*RT(IC)*WN(IC)
  430   CONTINUE
        DO 440 IC=IC1,IC2
        ROEPSR(IC)=ROEPSR(IC)-DT(IC)*(FL(IC)-FL(IC-NXYT))*RDZ(IC)
  440   CONTINUE
      END IF
C==============================================================CCFLUZ
      RETURN
      END
