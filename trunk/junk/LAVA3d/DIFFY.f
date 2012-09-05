*DECK DIFFY
      SUBROUTINE DIFFY
C===============================================================DIFFY
C     THIS ROUTINE CALCULATES DIFFUSION OF HEAT, CHEMICAL SPECIES,
C     AND SWIRL VELOCITY IN THE Y-DIRECTION
C----------
C     CALLED BY LAVA
C===============================================================DIFFY
      INCLUDE 'COML.h'
C----------
      DOUBLE PRECISION DTDY,DTEDY,FLEE,DWDY,
     1     TB,FR,AKAPP1,AKAPP2,AKAPPA,AMU,AMU1,AMU2,
     2     DCOEFF,ROCC,ROCF,DROCDY,
     3     ENTHC,ENTHF,ENTH,HSUBK,FRMC,FRMF,FRM,
     4     DKDY,DEDY,
     5     RECC,RECF,DRECDY,CRQMEF,SPDEF,SPDNF
C
C===============================================================DIFFY
      IF(INCOMP.EQ.1) GO TO 500
C========================================= HEAT CONDUCTION =====DIFFY
      DO 10 IC=IC1V,IC2
      AKAPP1=(DY(IC+NXT)+DY(IC))*COND(IC)*COND(IC+NXT)
      AKAPP2=DY(IC+NXT)*COND(IC)+DY(IC)*COND(IC+NXT)
      AKAPPA=AKAPP1/AKAPP2
      DTDY=(TEMP(IC+NXT)-TEMP(IC))*TWO*HRDYCF(IC)
      FL(IC)=-AKAPPA*DTDY*RF(IC)
      FLS(IC)=ZERO
      FLQ(IC)=ZERO
   10 CONTINUE
      DO 20 IC=IC1,IC2
      ROER(IC)=ROER(IC)-DT(IC)*(FL(IC)-FL(IC-NXT))*RDY(IC)
   20 CONTINUE
C--------------------------------------------- FOR NON-LTE -----DIFFY
      IF(NLTE.NE.0) THEN
        DO 50 IC=IC1V,IC2
        AKAPP1=(DY(IC+NXT)+DY(IC))*CONDE(IC)*CONDE(IC+NXT)
        AKAPP2=DY(IC+NXT)*CONDE(IC)+DY(IC)*CONDE(IC+NXT)
        AKAPPA=AKAPP1/AKAPP2
        DTEDY=(TE(IC+NXT)-TE(IC))*TWO*HRDYCF(IC)
        FL(IC)=-AKAPPA*DTEDY*RF(IC)
   50   CONTINUE
        DO 60 IC=IC1,IC2
        FLEE=-DT(IC)*(FL(IC)-FL(IC-NXT))*RDY(IC)
        ROER(IC)=ROER(IC)+FLEE
        ROEER(IC)=ROEER(IC)+FLEE
   60   CONTINUE
      END IF
C======================================= SPECIES DIFFUSION =====DIFFY
C     NOTE THAT ISP DOES NOT GO OVER IELC
C     SO WE DON'T HAVE TO USE DCOEF OF ELECTRON.
C---------------------------------------------------------------DIFFY
      IF(NSP.EQ.1 .OR. NODIFF.EQ.1) GO TO 500
C---------------------------------------------------------------DIFFY
C     USE WT ARRAY FOR TEMPERATURE AT THE FRONT CELL FACE.
C     USE DISS ARRAY FOR PRESSURE AT THE FRONT CELL FACE.
C---------------------------------------------------------------DIFFY
      IF(IELC.NE.0) THEN
        DO 100 IC=IC1V,IC2
        WT(IC)=(DY(IC)*TEMP(IC+NXT)+DY(IC+NXT)*TEMP(IC))*HRDYCF(IC)
        DISS(IC)=(DY(IC)*PN(IC+NXT)+DY(IC+NXT)*PN(IC))*HRDYCF(IC)
        IF(NLTE.EQ.0) THEN
          RECC=SPDN(IC,IELC)*RGAS*TEMP(IC)/PN(IC)
          RECF=SPDN(IC+NXT,IELC)*RGAS*TEMP(IC+NXT)/PN(IC+NXT)
        ELSE
          RECC=SPDN(IC,IELC)*RGAS*TE(IC)/PN(IC)
          RECF=SPDN(IC+NXT,IELC)*RGAS*TE(IC+NXT)/PN(IC+NXT)
        END IF
        DRECDY=TWO*(RECF-RECC)*HRDYCF(IC)
        SPDEF=(DY(IC)*SPDN(IC+NXT,IELC)+DY(IC+NXT)*SPDN(IC,IELC))
     1        *HRDYCF(IC)
        CRQMEF=DISS(IC)/(RGAS*WT(IC)*QOM(IELC)*MW(IELC)*SPDEF+SMALL)
C----- USE FLSQ TO STORE SOME CONSTANTS AT FRONT CELL FACE
        FLSQ(IC)=CRQMEF*DRECDY*RF(IC)
        FLAS(IC)=ZERO
        FLAQ(IC)=ZERO
  100   CONTINUE
      ELSE
        DO 110 IC=IC1V,IC2
        WT(IC)=(DY(IC)*TEMP(IC+NXT)+DY(IC+NXT)*TEMP(IC))*HRDYCF(IC)
        DISS(IC)=(DY(IC)*PN(IC+NXT)+DY(IC+NXT)*PN(IC))*HRDYCF(IC)
        FLSQ(IC)=ZERO
        FLAS(IC)=ZERO
        FLAQ(IC)=ZERO
  110   CONTINUE
      END IF
C---------------------------------------------------------------DIFFY
      IF(IELC.EQ.0 .OR. IELC.EQ.1) GO TO 240
      DO 210 ISP=1,IELC-1
        DO 205 IC=IC1V,IC2
      SPDNF=(DY(IC)*SPDN(IC+NXT,ISP)+DY(IC+NXT)*SPDN(IC,ISP))*HRDYCF(IC)
        DCOEFF=(DY(IC)*DCOEF(IC+NXT,ISP)+DY(IC+NXT)*DCOEF(IC,ISP))
     1         *HRDYCF(IC)
        TB=HUNDTH*WT(IC)
        IT=INT(TB)
        IT=MIN(IT,NIT1)
        FR=TB-DBLE(IT)
C---------------------------------------------------------------DIFFY
        ROCC=SPDN(IC,ISP)*RGAS*TEMP(IC)/PN(IC)
        ROCF=SPDN(IC+NXT,ISP)*RGAS*TEMP(IC+NXT)/PN(IC+NXT)
        DROCDY=TWO*(ROCF-ROCC)*HRDYCF(IC)
C---------------------------------------------------------------DIFFY
        FLA(IC)=MW(ISP)*SPDNF*QOM(ISP)*DCOEFF*FLSQ(IC)
        FLAS(IC)=FLAS(IC)-FLA(IC)
        FLAQ(IC)=FLAQ(IC)+FLA(IC)*QOM(ISP)
        FL(IC)=-(DISS(IC)/(RGAS*WT(IC)))*DCOEFF*DROCDY*RF(IC)
        FLS(IC)=FLS(IC)-FL(IC)
        FLQ(IC)=FLQ(IC)-FL(IC)*QOM(ISP)
C-------------------------------------- ENTHALPY DIFFUSION -----DIFFY
        ENTHC=(ROEN(IC)+PN(IC))/RON(IC)
        ENTHF=(ROEN(IC+NXT)+PN(IC+NXT))/RON(IC+NXT)
        ENTH=(DY(IC)*ENTHF+DY(IC+NXT)*ENTHC)*HRDYCF(IC)
        HSUBK=(ONE-FR)*EK(IT+1,ISP)+FR*EK(IT+2,ISP)
     1        +RGAS*WT(IC)*RMW(ISP)
        FLH(IC)=(HSUBK-ENTH)*(FL(IC)+FLA(IC))
  205   CONTINUE
        DO 206 IC=IC1,IC2
        SPDR(IC,ISP)=SPDR(IC,ISP)-DT(IC)*
     1               (FL(IC)+FLA(IC)-FL(IC-NXT)-FLA(IC-NXT))*RDY(IC)
        ROER(IC)=ROER(IC)-DT(IC)*(FLH(IC)-FLH(IC-NXT))*RDY(IC)
  206   CONTINUE
  210 CONTINUE
C---------------------------------------------------------------DIFFY
  240 IF(IELC.EQ.NSP) GO TO 270
      DO 260 ISP=IELC+1,NSP
        DO 255 IC=IC1V,IC2
      SPDNF=(DY(IC)*SPDN(IC+NXT,ISP)+DY(IC+NXT)*SPDN(IC,ISP))*HRDYCF(IC)
        DCOEFF=(DY(IC)*DCOEF(IC+NXT,ISP)+DY(IC+NXT)*DCOEF(IC,ISP))
     1         *HRDYCF(IC)
        TB=HUNDTH*WT(IC)
        IT=INT(TB)
        IT=MIN(IT,NIT1)
        FR=TB-DBLE(IT)
C---------------------------------------------------------------DIFFY
        ROCC=SPDN(IC,ISP)*RGAS*TEMP(IC)/PN(IC)
        ROCF=SPDN(IC+NXT,ISP)*RGAS*TEMP(IC+NXT)/PN(IC+NXT)
        DROCDY=TWO*(ROCF-ROCC)*HRDYCF(IC)
C---------------------------------------------------------------DIFFY
        FLA(IC)=MW(ISP)*SPDNF*QOM(ISP)*DCOEFF*FLSQ(IC)
        FLAS(IC)=FLAS(IC)-FLA(IC)
        FLAQ(IC)=FLAQ(IC)+FLA(IC)*QOM(ISP)
        FL(IC)=-(DISS(IC)/(RGAS*WT(IC)))*DCOEFF*DROCDY*RF(IC)
        FLS(IC)=FLS(IC)-FL(IC)
        FLQ(IC)=FLQ(IC)-FL(IC)*QOM(ISP)
C-------------------------------------- ENTHALPY DIFFUSION -----DIFFY
        ENTHC=(ROEN(IC)+PN(IC))/RON(IC)
        ENTHF=(ROEN(IC+NXT)+PN(IC+NXT))/RON(IC+NXT)
        ENTH=(DY(IC)*ENTHF+DY(IC+NXT)*ENTHC)*HRDYCF(IC)
        HSUBK=(ONE-FR)*EK(IT+1,ISP)+FR*EK(IT+2,ISP)
     1        +RGAS*WT(IC)*RMW(ISP)
        FLH(IC)=(HSUBK-ENTH)*(FL(IC)+FLA(IC))
  255   CONTINUE
        DO 256 IC=IC1,IC2
        SPDR(IC,ISP)=SPDR(IC,ISP)-DT(IC)*
     1               (FL(IC)+FLA(IC)-FL(IC-NXT)-FLA(IC-NXT))*RDY(IC)
        ROER(IC)=ROER(IC)-DT(IC)*(FLH(IC)-FLH(IC-NXT))*RDY(IC)
  256   CONTINUE
  260 CONTINUE
C---------------------------------------------------------------DIFFY
  270 IF(IELC.EQ.0 .OR. IELC.EQ.1) GO TO 280
      DO 279 ISP=1,IELC-1
        DO 273 IC=IC1V,IC2
        FRMC=SPDN(IC,ISP)/RON(IC)
        FRMF=SPDN(IC+NXT,ISP)/RON(IC+NXT)
        FRM=(DY(IC)*FRMF+DY(IC+NXT)*FRMC)*HRDYCF(IC)
        FL(IC)=FRM*(FLS(IC)+FLAS(IC))
  273   CONTINUE
        DO 276 IC=IC1,IC2
        SPDR(IC,ISP)=SPDR(IC,ISP)-DT(IC)*(FL(IC)-FL(IC-NXT))*RDY(IC)
  276 CONTINUE
  279 CONTINUE
C---------------------------------------------------------------DIFFY
  280 IF(IELC.EQ.NSP) GO TO 300
      DO 290 ISP=IELC+1,NSP
        DO 283 IC=IC1V,IC2
        FRMC=SPDN(IC,ISP)/RON(IC)
        FRMF=SPDN(IC+NXT,ISP)/RON(IC+NXT)
        FRM=(DY(IC)*FRMF+DY(IC+NXT)*FRMC)*HRDYCF(IC)
        FL(IC)=FRM*(FLS(IC)+FLAS(IC))
  283   CONTINUE
        DO 286 IC=IC1,IC2
        SPDR(IC,ISP)=SPDR(IC,ISP)-DT(IC)*(FL(IC)-FL(IC-NXT))*RDY(IC)
  286 CONTINUE
  290 CONTINUE
C-------------------------------------- ELECTRON DIFFUSION -----DIFFY
  300 IF(IELC.EQ.0) GO TO 500
      DO 410 IC=IC1V,IC2
      FRMC=SPDN(IC,IELC)/RON(IC)
      FRMF=SPDN(IC+NXT,IELC)/RON(IC+NXT)
      FRM=(DY(IC)*FRMF+DY(IC+NXT)*FRMC)*HRDYCF(IC)
      FL(IC)=FRM*(FLS(IC)+FLAS(IC))+(FLQ(IC)-FLAQ(IC))/QOM(IELC)
C------------ ENTHALPY DIFFUSION
      WT(IC)=(DY(IC)*TE(IC+NXT)+DY(IC+NXT)*TE(IC))*HRDYCF(IC)
      TB=HUNDTH*WT(IC)
      IT=INT(TB)
      IT=MIN(IT,NIT1)
      FR=TB-DBLE(IT)
      HSUBK=(ONE-FR)*EK(IT+1,IELC)+FR*EK(IT+2,IELC)
     1      +RGAS*WT(IC)*RMW(IELC)
      FLH(IC)=HSUBK*(FLQ(IC)-FLAQ(IC))/QOM(IELC)
  410 CONTINUE
      DO 420 IC=IC1,IC2
      SPDR(IC,IELC)=SPDR(IC,IELC)-DT(IC)*(FL(IC)-FL(IC-NXT))*RDY(IC)
      ROER(IC)=ROER(IC)-DT(IC)*(FLH(IC)-FLH(IC-NXT))*RDY(IC)
  420 CONTINUE
      IF(NLTE.NE.0) THEN
        DO 430 IC=IC1,IC2
        ROEER(IC)=ROEER(IC)-DT(IC)*(FLH(IC)-FLH(IC-NXT))*RDY(IC)
  430   CONTINUE
      END IF
C============================================ SWIRL STRESS =====DIFFY
  500 IF(ISWIRL.EQ.0) GO TO 590
      DO 510 IC=IC1V,IC2
      DWDY=(WN(IC+NXT)-WN(IC))*TWO*HRDYCF(IC)
      AMU1=(DY(IC+NXT)+DY(IC))*VISC(IC)*VISC(IC+NXT)
      AMU2=DY(IC+NXT)*VISC(IC)+DY(IC)*VISC(IC+NXT)
      AMU=AMU1/AMU2
      FL(IC)=AMU*DWDY*RF(IC)
      DISS(IC)=HALF*DT(IC)*DWDY*FL(IC)
  510 CONTINUE
C----------------- VISCOUS DISSIPATION DUE TO SWIRL STRESS -----DIFFY
      IF(ITURB.EQ.0) THEN
        DO 520 IC=IC1,IC2
        RORW(IC)=RORW(IC)+DT(IC)*(FL(IC)-FL(IC-NXT))*RDY(IC)
        ROER(IC)=ROER(IC)+DISS(IC)+DISS(IC-NXT)
  520   CONTINUE
      ELSEIF(ITURB.EQ.1) THEN
        DO 530 IC=IC1,IC2
        RORW(IC)=RORW(IC)+DT(IC)*(FL(IC)-FL(IC-NXT))*RDY(IC)
        ROTKER(IC)=ROTKER(IC)+DISS(IC)+DISS(IC-NXT)
  530   CONTINUE
      ELSEIF(ITURB.EQ.2) THEN
        DO 540 IC=IC1,IC2
        RORW(IC)=RORW(IC)+DT(IC)*(FL(IC)-FL(IC-NXT))*RDY(IC)
        ROTKER(IC)=ROTKER(IC)+DISS(IC)+DISS(IC-NXT)
        ROEPSR(IC)=ROEPSR(IC)
     1             +(DISS(IC)+DISS(IC-NXT))*CE1*EPSN(IC)/TKEN(IC)
  540   CONTINUE
      END IF
  590 CONTINUE
C============================ FOR TURBULENT KINETIC ENERGY =====DIFFY
      IF(ITURB.NE.0) THEN
        DO 600 IC=IC1V,IC2
        DKDY=(TKEN(IC+NXT)-TKEN(IC))*TWO*HRDYCF(IC)
        AMU1=(DY(IC+NXT)+DY(IC))*VISC(IC)*VISC(IC+NXT)
        AMU2=DY(IC+NXT)*VISC(IC)+DY(IC)*VISC(IC+NXT)
        AMU=AMU1/AMU2
        FL(IC)=-AMU*DKDY*RF(IC)
  600   CONTINUE
        DO 610 IC=IC1,IC2
        ROTKER(IC)=ROTKER(IC)-DT(IC)*(FL(IC)-FL(IC-NXT))*RDY(IC)
  610   CONTINUE
      END IF
C------------- FOR DISSIPATION OF TURBULENT KINETIC ENERGY -----DIFFY
      IF(ITURB.EQ.2) THEN
        DO 700 IC=IC1V,IC2
        DEDY=(EPSN(IC+NXT)-EPSN(IC))*TWO*HRDYCF(IC)
        AMU1=(DY(IC+NXT)+DY(IC))*(VISC(IC)+(RPRE-ONE)*VIST(IC))
     1       *(VISC(IC+NXT)+(RPRE-ONE)*VIST(IC+NXT))
        AMU2=DY(IC+NXT)*(VISC(IC)+(RPRE-ONE)*VIST(IC))
     1       +DY(IC)*(VISC(IC+NXT)+(RPRE-ONE)*VIST(IC+NXT))
        AMU=AMU1/AMU2
        FL(IC)=-AMU*DEDY*RF(IC)
  700   CONTINUE
        DO 710 IC=IC1,IC2
        ROEPSR(IC)=ROEPSR(IC)-DT(IC)*(FL(IC)-FL(IC-NXT))*RDY(IC)
  710   CONTINUE
      END IF
C===============================================================DIFFY
      RETURN
      END
