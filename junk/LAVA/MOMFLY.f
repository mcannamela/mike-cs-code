*DECK MOMFLY
      SUBROUTINE MOMFLY
C==============================================================MOMFLY
C     THIS ROUTINE UPDATES RORU, RORV, AND RORW
C     DUE TO CONVECTIVE FLUXING OF MOMENTUM IN THE Y DIRECTION
C----------
C     CALLED BY LAVA
C----------
C     INPUT: ROVU, ROVV, ROVW, DX, DY, DZ, RC, RF
C     OUTPUT: RORU, RORV, RORW
C==============================================================MOMFLY
      INCLUDE 'COML.h'
C----------
      DOUBLE PRECISION
     1     RDYCF,VC,SS,COUR,WTC,VNFR,VNFT,DELF,RONC,ROND,RONF,
     2     ROUC,ROUF,ROVF,ROVD,ROWF,ROWC,RFT,RFR,SSD,SSF,DRD,DRF
C
C==============================================================MOMFLY
      IF(NX.EQ.1) GO TO 200
      DO 110 IC=IC1V-1,IC2
      VNFR=(DX(IC)*VN(IC+1)+DX(IC+1)*VN(IC))*HRDXCR(IC)
      SS=SIGN(ONE,VNFR)
      RDYCF=TWO*HRDYCF(IC)
      DELF=DY(IC+NXT)*RDYCF
      COUR=VNFR*DT(IC)*RDYCF
      WTC=HALF*(DELF+ADC*(ONE-DELF+SS)+BDC*COUR)
      RONF=(DX(IC)*RON(IC+NXT+1)+DX(IC+1)*RON(IC+NXT))*HRDXCR(IC)
      RONC=(DX(IC)*RON(IC+1)+DX(IC+1)*RON(IC))*HRDXCR(IC)
      ROUF=RONF*UN(IC+NXT)
      ROUC=RONC*UN(IC)
      RFR=(DX(IC)*RF(IC+1)+DX(IC+1)*RF(IC))*HRDXCR(IC)
      RORVU(IC)=(ROUC*WTC+ROUF*(ONE-WTC))*RFR*VNFR
  110 CONTINUE
      DO 140 IC=IC1U,IC2
      RORU(IC)=RORU(IC)-DT(IC)*(RORVU(IC)-RORVU(IC-NXT))*RDY(IC)
  140 CONTINUE
  200 CONTINUE
C==============================================================MOMFLY
      DO 10 IC=IC1,IC2
      VC=HALF*(VN(IC)+VN(IC-NXT))
      SS=SIGN(ONE,VC)
      COUR=VC*DT(IC)*RDY(IC)
      WTC=HALF*(ONE+ADC*SS+BDC*COUR)
      ROND=(DY(IC)*RON(IC-NXT)+DY(IC-NXT)*RON(IC))*HRDYCF(IC-NXT)
      RONF=(DY(IC)*RON(IC+NXT)+DY(IC+NXT)*RON(IC))*HRDYCF(IC)
      ROVD=ROND*VN(IC-NXT)
      ROVF=RONF*VN(IC)
      RORVV(IC)=(ROVD*WTC+ROVF*(ONE-WTC))*RC(IC)*VC
   10 CONTINUE
C---------- EXTRAPOLATE OUTFLOW BOUNDARY ----------------------MOMFLY
      DO 30 K=1,NZT
      ICK=(K-1)*NXYT
           DO 20 I=1,NXT
           ICD=I+ICK
           ICF=ICD+NXYT-NXT
           SSD=SIGN(ONE,V(ICD))
           SSF=SIGN(ONE,-V(ICF-NXT))
           DRD=(DY(ICD)+DY(ICD+NXT))/(DY(ICD+NXT)+DY(ICD+NXT+NXT))
           DRF=(DY(ICF)+DY(ICF-NXT))/(DY(ICF-NXT)+DY(ICF-NXT-NXT))
           RORVV(ICD)=RORVV(ICD+NXT)+HALF*(ONE-SSD)*DRD
     1                *(RORVV(ICD+NXT)-RORVV(ICD+NXT+NXT))
           RORVV(ICF)=RORVV(ICF-NXT)+HALF*(ONE-SSF)*DRF
     1                *(RORVV(ICF-NXT)-RORVV(ICF-NXT-NXT))
   20 CONTINUE
   30 CONTINUE
C----------
      DO 40 IC=IC1V,IC2
      RORV(IC)=RORV(IC)-DT(IC)*(RORVV(IC+NXT)-RORVV(IC))*TWO*HRDYCF(IC)
   40 CONTINUE
C==============================================================MOMFLY
      IF(NZ.EQ.1) GO TO 300
      DO 210 IC=IC1V-NXYT,IC2
      VNFT=(DZ(IC)*VN(IC+NXYT)+DZ(IC+NXYT)*VN(IC))*HRDZCT(IC)
      SS=SIGN(ONE,VNFT)
      RDYCF=TWO*HRDYCF(IC)
      DELF=DY(IC+NXT)*RDYCF
      COUR=VNFT*DT(IC)*RDYCF
      WTC=HALF*(DELF+ADC*(ONE-DELF+SS)+BDC*COUR)
      RONF=(DZ(IC)*RON(IC+NXYT+NXT)+DZ(IC+NXYT)*RON(IC+NXT))*HRDZCT(IC)
      RONC=(DZ(IC)*RON(IC+NXYT)+DZ(IC+NXYT)*RON(IC))*HRDZCT(IC)
      ROWF=RONF*WN(IC+NXT)
      ROWC=RONC*WN(IC)
      RFT=(DZ(IC)*RF(IC+NXYT)+DZ(IC+NXYT)*RF(IC))*HRDZCT(IC)
      RORVW(IC)=(ROWC*WTC+ROWF*(ONE-WTC))*RFT*VNFT
  210 CONTINUE
      DO 240 IC=IC1W,IC2
      RORW(IC)=RORW(IC)-DT(IC)*(RORVW(IC)-RORVW(IC-NXT))*RDY(IC)
  240 CONTINUE
  300 CONTINUE
C==============================================================MOMFLY
      RETURN
      END
