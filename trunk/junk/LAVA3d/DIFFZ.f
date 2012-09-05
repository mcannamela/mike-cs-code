*DECK DIFFZ
      SUBROUTINE DIFFZ
C===============================================================DIFFZ
C     THIS ROUTINE CALCULATES DIFFUSION OF HEAT AND CHEMICAL SPECIES
C     IN THE Z-DIRECTION
C----------
C     CALLED BY LAVA
C===============================================================DIFFZ
      INCLUDE 'COML.h'
C----------
      DOUBLE PRECISION DTDZ,DTEDZ,FLEE,
     1     TB,FR,AKAPP1,AKAPP2,AKAPPA,AMU,AMU1,AMU2,
     2     DCOEFT,ROCC,ROCT,DROCDZ,
     3     ENTHC,ENTHT,ENTH,HSUBK,FRMC,FRMT,FRM,
     4     DKDZ,DEDZ,
     5     RECC,RECT,DRECDZ,CRQMET,SPDET,SPDNT
C
C===============================================================DIFFZ
      IF(INCOMP.EQ.1) GO TO 500
C========================================= HEAT CONDUCTION =====DIFFZ
      DO 10 IC=IC1W,IC2
      AKAPP1=(DZ(IC+NXYT)+DZ(IC))*COND(IC)*COND(IC+NXYT)
      AKAPP2=DZ(IC+NXYT)*COND(IC)+DZ(IC)*COND(IC+NXYT)
      AKAPPA=AKAPP1/AKAPP2
      DTDZ=(TEMP(IC+NXYT)-TEMP(IC))*TWO*HRDZCT(IC)
      FL(IC)=-AKAPPA*DTDZ*RT(IC)
      FLS(IC)=ZERO
      FLQ(IC)=ZERO
   10 CONTINUE
      DO 20 IC=IC1,IC2
      ROER(IC)=ROER(IC)-DT(IC)*(FL(IC)-FL(IC-NXYT))*RDZ(IC)
   20 CONTINUE
C--------------------------------------------- FOR NON-LTE -----DIFFZ
      IF(NLTE.NE.0) THEN
        DO 50 IC=IC1W,IC2
        AKAPP1=(DZ(IC+NXYT)+DZ(IC))*CONDE(IC)*CONDE(IC+NXYT)
        AKAPP2=DZ(IC+NXYT)*CONDE(IC)+DZ(IC)*CONDE(IC+NXYT)
        AKAPPA=AKAPP1/AKAPP2
        DTEDZ=(TE(IC+NXYT)-TE(IC))*TWO*HRDZCT(IC)
        FL(IC)=-AKAPPA*DTEDZ*RT(IC)
   50   CONTINUE
        DO 60 IC=IC1,IC2
        FLEE=-DT(IC)*(FL(IC)-FL(IC-NXYT))*RDZ(IC)
        ROER(IC)=ROER(IC)+FLEE
        ROEER(IC)=ROEER(IC)+FLEE
   60   CONTINUE
      END IF
C======================================= SPECIES DIFFUSION =====DIFFZ
C     NOTE THAT ISP DOES NOT GO OVER IELC
C     SO WE DON'T HAVE TO USE DCOEF OF ELECTRON.
C---------------------------------------------------------------DIFFZ
      IF(NSP.EQ.1 .OR. NODIFF.EQ.1) GO TO 500
C---------------------------------------------------------------DIFFZ
C     USE WT ARRAY FOR TEMPERATURE AT THE TOP CELL FACE.
C     USE DISS ARRAY FOR PRESSURE AT THE TOP CELL FACE.
C---------------------------------------------------------------DIFFZ
      IF(IELC.NE.0) THEN
        DO 100 IC=IC1W,IC2
        WT(IC)=(DZ(IC)*TEMP(IC+NXYT)+DZ(IC+NXYT)*TEMP(IC))*HRDZCT(IC)
        DISS(IC)=(DZ(IC)*PN(IC+NXYT)+DZ(IC+NXYT)*PN(IC))*HRDZCT(IC)
        IF(NLTE.EQ.0) THEN
          RECC=SPDN(IC,IELC)*RGAS*TEMP(IC)/PN(IC)
          RECT=SPDN(IC+NXYT,IELC)*RGAS*TEMP(IC+NXYT)/PN(IC+NXYT)
        ELSE
          RECC=SPDN(IC,IELC)*RGAS*TE(IC)/PN(IC)
          RECT=SPDN(IC+NXYT,IELC)*RGAS*TE(IC+NXYT)/PN(IC+NXYT)
        END IF
        DRECDZ=TWO*(RECT-RECC)*HRDZCT(IC)
        SPDET=(DZ(IC)*SPDN(IC+NXYT,IELC)+DZ(IC+NXYT)*SPDN(IC,IELC))
     1        *HRDZCT(IC)
        CRQMET=DISS(IC)/(RGAS*WT(IC)*QOM(IELC)*MW(IELC)*SPDET+SMALL)
C----- USE FLSQ TO STORE SOME CONSTANTS AT TOP CELL FACE
        FLSQ(IC)=CRQMET*DRECDZ*RT(IC)
        FLAS(IC)=ZERO
        FLAQ(IC)=ZERO
  100   CONTINUE
      ELSE
        DO 110 IC=IC1W,IC2
        WT(IC)=(DZ(IC)*TEMP(IC+NXYT)+DZ(IC+NXYT)*TEMP(IC))*HRDZCT(IC)
        DISS(IC)=(DZ(IC)*PN(IC+NXYT)+DZ(IC+NXYT)*PN(IC))*HRDZCT(IC)
        FLSQ(IC)=ZERO
        FLAS(IC)=ZERO
        FLAQ(IC)=ZERO
  110   CONTINUE
      END IF
C---------------------------------------------------------------DIFFZ
      IF(IELC.EQ.0 .OR. IELC.EQ.1) GO TO 240
      DO 210 ISP=1,IELC-1
        DO 205 IC=IC1W,IC2
        SPDNT=(DZ(IC)*SPDN(IC+NXYT,ISP)+DZ(IC+NXYT)*SPDN(IC,ISP))
     1        *HRDZCT(IC)
        DCOEFT=(DZ(IC)*DCOEF(IC+NXYT,ISP)+DZ(IC+NXYT)*DCOEF(IC,ISP))
     1         *HRDZCT(IC)
        TB=HUNDTH*WT(IC)
        IT=INT(TB)
        IT=MIN(IT,NIT1)
        FR=TB-DBLE(IT)
C---------------------------------------------------------------DIFFZ
        ROCC=SPDN(IC,ISP)*RGAS*TEMP(IC)/PN(IC)
        ROCT=SPDN(IC+NXYT,ISP)*RGAS*TEMP(IC+NXYT)/PN(IC+NXYT)
        DROCDZ=TWO*(ROCT-ROCC)*HRDZCT(IC)
C---------------------------------------------------------------DIFFZ
        FLA(IC)=MW(ISP)*SPDNT*QOM(ISP)*DCOEFT*FLSQ(IC)
        FLAS(IC)=FLAS(IC)-FLA(IC)
        FLAQ(IC)=FLAQ(IC)+FLA(IC)*QOM(ISP)
        FL(IC)=-(DISS(IC)/(RGAS*WT(IC)))*DCOEFT*DROCDZ*RT(IC)
        FLS(IC)=FLS(IC)-FL(IC)
        FLQ(IC)=FLQ(IC)-FL(IC)*QOM(ISP)
C-------------------------------------- ENTHALPY DIFFUSION -----DIFFZ
        ENTHC=(ROEN(IC)+PN(IC))/RON(IC)
        ENTHT=(ROEN(IC+NXYT)+PN(IC+NXYT))/RON(IC+NXYT)
        ENTH=(DZ(IC)*ENTHT+DZ(IC+NXYT)*ENTHC)*HRDZCT(IC)
        HSUBK=(ONE-FR)*EK(IT+1,ISP)+FR*EK(IT+2,ISP)
     1        +RGAS*WT(IC)*RMW(ISP)
        FLH(IC)=(HSUBK-ENTH)*(FL(IC)+FLA(IC))
  205   CONTINUE
        DO 206 IC=IC1,IC2
        SPDR(IC,ISP)=SPDR(IC,ISP)-DT(IC)*
     1              (FL(IC)+FLA(IC)-FL(IC-NXYT)-FLA(IC-NXYT))*RDZ(IC)
        ROER(IC)=ROER(IC)-DT(IC)*(FLH(IC)-FLH(IC-NXYT))*RDZ(IC)
  206   CONTINUE
  210 CONTINUE
C---------------------------------------------------------------DIFFZ
  240 IF(IELC.EQ.NSP) GO TO 270
      DO 260 ISP=IELC+1,NSP
        DO 255 IC=IC1W,IC2
        SPDNT=(DZ(IC)*SPDN(IC+NXYT,ISP)+DZ(IC+NXYT)*SPDN(IC,ISP))
     1        *HRDZCT(IC)
        DCOEFT=(DZ(IC)*DCOEF(IC+NXYT,ISP)+DZ(IC+NXYT)*DCOEF(IC,ISP))
     1         *HRDZCT(IC)
        TB=HUNDTH*WT(IC)
        IT=INT(TB)
        IT=MIN(IT,NIT1)
        FR=TB-DBLE(IT)
C---------------------------------------------------------------DIFFZ
        ROCC=SPDN(IC,ISP)*RGAS*TEMP(IC)/PN(IC)
        ROCT=SPDN(IC+NXYT,ISP)*RGAS*TEMP(IC+NXYT)/PN(IC+NXYT)
        DROCDZ=TWO*(ROCT-ROCC)*HRDZCT(IC)
C---------------------------------------------------------------DIFFZ
        FLA(IC)=MW(ISP)*SPDNT*QOM(ISP)*DCOEFT*FLSQ(IC)
        FLAS(IC)=FLAS(IC)-FLA(IC)
        FLAQ(IC)=FLAQ(IC)+FLA(IC)*QOM(ISP)
        FL(IC)=-(DISS(IC)/(RGAS*WT(IC)))*DCOEFT*DROCDZ*RT(IC)
        FLS(IC)=FLS(IC)-FL(IC)
        FLQ(IC)=FLQ(IC)-FL(IC)*QOM(ISP)
C-------------------------------------- ENTHALPY DIFFUSION -----DIFFZ
        ENTHC=(ROEN(IC)+PN(IC))/RON(IC)
        ENTHT=(ROEN(IC+NXYT)+PN(IC+NXYT))/RON(IC+NXYT)
        ENTH=(DZ(IC)*ENTHT+DZ(IC+NXYT)*ENTHC)*HRDZCT(IC)
        HSUBK=(ONE-FR)*EK(IT+1,ISP)+FR*EK(IT+2,ISP)
     1        +RGAS*WT(IC)*RMW(ISP)
        FLH(IC)=(HSUBK-ENTH)*(FL(IC)+FLA(IC))
  255   CONTINUE
        DO 256 IC=IC1,IC2
        SPDR(IC,ISP)=SPDR(IC,ISP)-DT(IC)*
     1              (FL(IC)+FLA(IC)-FL(IC-NXYT)-FLA(IC-NXYT))*RDZ(IC)
        ROER(IC)=ROER(IC)-DT(IC)*(FLH(IC)-FLH(IC-NXYT))*RDZ(IC)
  256   CONTINUE
  260 CONTINUE
C---------------------------------------------------------------DIFFZ
  270 IF(IELC.EQ.0 .OR. IELC.EQ.1) GO TO 280
      DO 279 ISP=1,IELC-1
        DO 273 IC=IC1W,IC2
        FRMC=SPDN(IC,ISP)/RON(IC)
        FRMT=SPDN(IC+NXYT,ISP)/RON(IC+NXYT)
        FRM=(DZ(IC)*FRMT+DZ(IC+NXYT)*FRMC)*HRDZCT(IC)
        FL(IC)=FRM*(FLS(IC)+FLAS(IC))
  273   CONTINUE
        DO 276 IC=IC1,IC2
        SPDR(IC,ISP)=SPDR(IC,ISP)-DT(IC)*(FL(IC)-FL(IC-NXYT))*RDZ(IC)
  276 CONTINUE
  279 CONTINUE
C---------------------------------------------------------------DIFFZ
  280 IF(IELC.EQ.NSP) GO TO 300
      DO 290 ISP=IELC+1,NSP
        DO 283 IC=IC1W,IC2
        FRMC=SPDN(IC,ISP)/RON(IC)
        FRMT=SPDN(IC+NXYT,ISP)/RON(IC+NXYT)
        FRM=(DZ(IC)*FRMT+DZ(IC+NXYT)*FRMC)*HRDZCT(IC)
        FL(IC)=FRM*(FLS(IC)+FLAS(IC))
  283   CONTINUE
        DO 286 IC=IC1,IC2
        SPDR(IC,ISP)=SPDR(IC,ISP)-DT(IC)*(FL(IC)-FL(IC-NXYT))*RDZ(IC)
  286 CONTINUE
  290 CONTINUE
C-------------------------------------- ELECTRON DIFFUSION -----DIFFZ
  300 IF(IELC.EQ.0) GO TO 500
      DO 410 IC=IC1W,IC2
      FRMC=SPDN(IC,IELC)/RON(IC)
      FRMT=SPDN(IC+NXYT,IELC)/RON(IC+NXYT)
      FRM=(DZ(IC)*FRMT+DZ(IC+NXYT)*FRMC)*HRDZCT(IC)
      FL(IC)=FRM*(FLS(IC)+FLAS(IC))+(FLQ(IC)-FLAQ(IC))/QOM(IELC)
C------------ ENTHALPY DIFFUSION
      WT(IC)=(DZ(IC)*TE(IC+NXYT)+DZ(IC+NXYT)*TE(IC))*HRDZCT(IC)
      TB=HUNDTH*WT(IC)
      IT=INT(TB)
      IT=MIN(IT,NIT1)
      FR=TB-DBLE(IT)
      HSUBK=(ONE-FR)*EK(IT+1,IELC)+FR*EK(IT+2,IELC)
     1      +RGAS*WT(IC)*RMW(IELC)
      FLH(IC)=HSUBK*(FLQ(IC)-FLAQ(IC))/QOM(IELC)
  410 CONTINUE
      DO 420 IC=IC1,IC2
      SPDR(IC,IELC)=SPDR(IC,IELC)-DT(IC)*(FL(IC)-FL(IC-NXYT))*RDZ(IC)
      ROER(IC)=ROER(IC)-DT(IC)*(FLH(IC)-FLH(IC-NXYT))*RDZ(IC)
  420 CONTINUE
      IF(NLTE.NE.0) THEN
        DO 430 IC=IC1,IC2
        ROEER(IC)=ROEER(IC)-DT(IC)*(FLH(IC)-FLH(IC-NXYT))*RDZ(IC)
  430   CONTINUE
      END IF
C============================ FOR TURBULENT KINETIC ENERGY =====DIFFZ
  500 IF(ITURB.NE.0) THEN
        DO 600 IC=IC1W,IC2
        DKDZ=(TKEN(IC+NXYT)-TKEN(IC))*TWO*HRDZCT(IC)
        AMU1=(DZ(IC+NXYT)+DZ(IC))*VISC(IC)*VISC(IC+NXYT)
        AMU2=DZ(IC+NXYT)*VISC(IC)+DZ(IC)*VISC(IC+NXYT)
        AMU=AMU1/AMU2
        FL(IC)=-AMU*DKDZ*RT(IC)
  600   CONTINUE
        DO 610 IC=IC1,IC2
        ROTKER(IC)=ROTKER(IC)-DT(IC)*(FL(IC)-FL(IC-NXYT))*RDZ(IC)
  610   CONTINUE
      END IF
C------------- FOR DISSIPATION OF TURBULENT KINETIC ENERGY -----DIFFZ
      IF(ITURB.EQ.2) THEN
        DO 700 IC=IC1W,IC2
        DEDZ=(EPSN(IC+NXYT)-EPSN(IC))*TWO*HRDZCT(IC)
        AMU1=(DZ(IC+NXYT)+DZ(IC))*(VISC(IC)+(RPRE-ONE)*VIST(IC))
     1       *(VISC(IC+NXYT)+(RPRE-ONE)*VIST(IC+NXYT))
        AMU2=DZ(IC+NXYT)*(VISC(IC)+(RPRE-ONE)*VIST(IC))
     1       +DZ(IC)*(VISC(IC+NXYT)+(RPRE-ONE)*VIST(IC+NXYT))
        AMU=AMU1/AMU2
        FL(IC)=-AMU*DEDZ*RT(IC)
  700   CONTINUE
        DO 710 IC=IC1,IC2
        ROEPSR(IC)=ROEPSR(IC)-DT(IC)*(FL(IC)-FL(IC-NXYT))*RDZ(IC)
  710   CONTINUE
      END IF
C===============================================================DIFFZ
      RETURN
      END
