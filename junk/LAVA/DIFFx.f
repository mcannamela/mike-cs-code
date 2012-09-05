*DECK DIFFX
      SUBROUTINE DIFFX
C===============================================================DIFFX
C     THIS ROUTINE CALCULATES DIFFUSION OF HEAT, CHEMICAL SPECIES,
C     AND SWIRL VELOCITY IN THE X-DIRECTION
C----------
C     CALLED BY LAVA
C===============================================================DIFFX
      INCLUDE 'COML.h'
C----------
      DOUBLE PRECISION DTDX,DTEDX,FLEE,DWOXDX,XR,
     1     TB,FR,AKAPP1,AKAPP2,AKAPPA,AMU,AMU1,AMU2,
     2     DCOEFR,ROCC,ROCR,DROCDX,
     3     ENTHC,ENTHR,ENTH,HSUBK,FRMC,FRMR,FRM,
     4     DKDX,DEDX,RHOR,
     5     RECC,RECR,DRECDX,CRQMER,SPDER,SPDNR
C
C===============================================================DIFFX
      IF(INCOMP.EQ.1) GO TO 500
C========================================= HEAT CONDUCTION =====DIFFX
C     IN CYLINDRICAL COORDINATES, HARMONIC MEANS OF THERMAL
C     CONDUCTIVITIES AS IN SLABS ARE USED (LOG IS NOT USED).
C     IN Y AND Z DIRECTION, THIS APPROACH IS CORRECT.
C---------------------------------------------------------------DIFFX
      DO 10 IC=IC1U,IC2
      AKAPP1=(DX(IC+1)+DX(IC))*COND(IC)*COND(IC+1)
      AKAPP2=DX(IC+1)*COND(IC)+DX(IC)*COND(IC+1)
      AKAPPA=AKAPP1/AKAPP2
      DTDX=(TEMP(IC+1)-TEMP(IC))*TWO*HRDXCR(IC)
      FL(IC)=-AKAPPA*DTDX*RR(IC)
      FLS(IC)=ZERO
      FLQ(IC)=ZERO
   10 CONTINUE
      DO 20 IC=IC1,IC2
      ROER(IC)=ROER(IC)-DT(IC)*(FL(IC)-FL(IC-1))*RDX(IC)
   20 CONTINUE
C--------------------------------------------- FOR NON-LTE -----DIFFX
      IF(NLTE.NE.0) THEN
        DO 50 IC=IC1U,IC2
        AKAPP1=(DX(IC+1)+DX(IC))*CONDE(IC)*CONDE(IC+1)
        AKAPP2=DX(IC+1)*CONDE(IC)+DX(IC)*CONDE(IC+1)
        AKAPPA=AKAPP1/AKAPP2
        DTEDX=(TE(IC+1)-TE(IC))*TWO*HRDXCR(IC)
        FL(IC)=-AKAPPA*DTEDX*RR(IC)
   50   CONTINUE
        DO 60 IC=IC1,IC2
        FLEE=-DT(IC)*(FL(IC)-FL(IC-1))*RDX(IC)
        ROER(IC)=ROER(IC)+FLEE
        ROEER(IC)=ROEER(IC)+FLEE
   60   CONTINUE
      END IF
C======================================= SPECIES DIFFUSION =====DIFFX
C     NOTE THAT ISP DOES NOT GO OVER IELC
C     SO WE DON'T HAVE TO USE DCOEF OF ELECTRON.
C---------------------------------------------------------------DIFFX
      IF(NSP.EQ.1 .OR. NODIFF.EQ.1) GO TO 500
C---------------------------------------------------------------DIFFX
C     USE WT ARRAY FOR TEMPERATURE AT THE RIGHT CELL FACE.
C     USE DISS ARRAY FOR PRESSURE AT THE RIGHT CELL FACE.
C---------------------------------------------------------------DIFFX
      IF(IELC.NE.0) THEN
        DO 100 IC=IC1U,IC2
        WT(IC)=(DX(IC)*TEMP(IC+1)+DX(IC+1)*TEMP(IC))*HRDXCR(IC)
        DISS(IC)=(DX(IC)*PN(IC+1)+DX(IC+1)*PN(IC))*HRDXCR(IC)
        IF(NLTE.EQ.0) THEN
          RECC=SPDN(IC,IELC)*RGAS*TEMP(IC)/PN(IC)
          RECR=SPDN(IC+1,IELC)*RGAS*TEMP(IC+1)/PN(IC+1)
        ELSE
          RECC=SPDN(IC,IELC)*RGAS*TE(IC)/PN(IC)
          RECR=SPDN(IC+1,IELC)*RGAS*TE(IC+1)/PN(IC+1)
        END IF
        DRECDX=TWO*(RECR-RECC)*HRDXCR(IC)
        SPDER=(DX(IC)*SPDN(IC+1,IELC)+DX(IC+1)*SPDN(IC,IELC))*HRDXCR(IC)
        CRQMER=DISS(IC)/(RGAS*WT(IC)*QOM(IELC)*MW(IELC)*SPDER+SMALL)
C----- USE FLSQ TO STORE SOME CONSTANTS AT RIGHT CELL FACE
        FLSQ(IC)=CRQMER*DRECDX*RR(IC)
        FLAS(IC)=ZERO
        FLAQ(IC)=ZERO
  100   CONTINUE
      ELSE
        DO 110 IC=IC1U,IC2
        WT(IC)=(DX(IC)*TEMP(IC+1)+DX(IC+1)*TEMP(IC))*HRDXCR(IC)
        DISS(IC)=(DX(IC)*PN(IC+1)+DX(IC+1)*PN(IC))*HRDXCR(IC)
        FLSQ(IC)=ZERO
        FLAS(IC)=ZERO
        FLAQ(IC)=ZERO
  110   CONTINUE
      END IF
C---------------------------------------------------------------DIFFX
      IF(IELC.EQ.0 .OR. IELC.EQ.1) GO TO 240
      DO 210 ISP=1,IELC-1
        DO 205 IC=IC1U,IC2
        SPDNR=(DX(IC)*SPDN(IC+1,ISP)+DX(IC+1)*SPDN(IC,ISP))*HRDXCR(IC)
        DCOEFR=(DX(IC)*DCOEF(IC+1,ISP)+DX(IC+1)*DCOEF(IC,ISP))
     1         *HRDXCR(IC)
        TB=HUNDTH*WT(IC)
        IT=INT(TB)
        IT=MIN(IT,NIT1)
        FR=TB-DBLE(IT)
C---------------------------------------------------------------DIFFX
        ROCC=SPDN(IC,ISP)*RGAS*TEMP(IC)/PN(IC)
        ROCR=SPDN(IC+1,ISP)*RGAS*TEMP(IC+1)/PN(IC+1)
        DROCDX=TWO*(ROCR-ROCC)*HRDXCR(IC)
C---------------------------------------------------------------DIFFX
        FLA(IC)=MW(ISP)*SPDNR*QOM(ISP)*DCOEFR*FLSQ(IC)
        FLAS(IC)=FLAS(IC)-FLA(IC)
        FLAQ(IC)=FLAQ(IC)+FLA(IC)*QOM(ISP)
        FL(IC)=-(DISS(IC)/(RGAS*WT(IC)))*DCOEFR*DROCDX*RR(IC)
        FLS(IC)=FLS(IC)-FL(IC)
        FLQ(IC)=FLQ(IC)-FL(IC)*QOM(ISP)
C-------------------------------------- ENTHALPY DIFFUSION -----DIFFX
        ENTHC=(ROEN(IC)+PN(IC))/RON(IC)
        ENTHR=(ROEN(IC+1)+PN(IC+1))/RON(IC+1)
        ENTH=(DX(IC)*ENTHR+DX(IC+1)*ENTHC)*HRDXCR(IC)
        HSUBK=(ONE-FR)*EK(IT+1,ISP)+FR*EK(IT+2,ISP)
     1        +RGAS*WT(IC)*RMW(ISP)
        FLH(IC)=(HSUBK-ENTH)*(FL(IC)+FLA(IC))
  205   CONTINUE
        DO 206 IC=IC1,IC2
        SPDR(IC,ISP)=SPDR(IC,ISP)-DT(IC)*
     1               (FL(IC)+FLA(IC)-FL(IC-1)-FLA(IC-1))*RDX(IC)
        ROER(IC)=ROER(IC)-DT(IC)*(FLH(IC)-FLH(IC-1))*RDX(IC)
  206   CONTINUE
  210 CONTINUE
C---------------------------------------------------------------DIFFX
  240 IF(IELC.EQ.NSP) GO TO 270
      DO 260 ISP=IELC+1,NSP
        DO 255 IC=IC1U,IC2
        SPDNR=(DX(IC)*SPDN(IC+1,ISP)+DX(IC+1)*SPDN(IC,ISP))*HRDXCR(IC)
        DCOEFR=(DX(IC)*DCOEF(IC+1,ISP)+DX(IC+1)*DCOEF(IC,ISP))
     1         *HRDXCR(IC)
        TB=HUNDTH*WT(IC)
        IT=INT(TB)
        IT=MIN(IT,NIT1)
        FR=TB-DBLE(IT)
C---------------------------------------------------------------DIFFX
        ROCC=SPDN(IC,ISP)*RGAS*TEMP(IC)/PN(IC)
        ROCR=SPDN(IC+1,ISP)*RGAS*TEMP(IC+1)/PN(IC+1)
        DROCDX=TWO*(ROCR-ROCC)*HRDXCR(IC)
C---------------------------------------------------------------DIFFX
        FLA(IC)=MW(ISP)*SPDNR*QOM(ISP)*DCOEFR*FLSQ(IC)
        FLAS(IC)=FLAS(IC)-FLA(IC)
        FLAQ(IC)=FLAQ(IC)+FLA(IC)*QOM(ISP)
        FL(IC)=-(DISS(IC)/(RGAS*WT(IC)))*DCOEFR*DROCDX*RR(IC)
        FLS(IC)=FLS(IC)-FL(IC)
        FLQ(IC)=FLQ(IC)-FL(IC)*QOM(ISP)
C-------------------------------------- ENTHALPY DIFFUSION -----DIFFX
        ENTHC=(ROEN(IC)+PN(IC))/RON(IC)
        ENTHR=(ROEN(IC+1)+PN(IC+1))/RON(IC+1)
        ENTH=(DX(IC)*ENTHR+DX(IC+1)*ENTHC)*HRDXCR(IC)
        HSUBK=(ONE-FR)*EK(IT+1,ISP)+FR*EK(IT+2,ISP)
     1        +RGAS*WT(IC)*RMW(ISP)
        FLH(IC)=(HSUBK-ENTH)*(FL(IC)+FLA(IC))
  255   CONTINUE
        DO 256 IC=IC1,IC2
        SPDR(IC,ISP)=SPDR(IC,ISP)-DT(IC)*
     1               (FL(IC)+FLA(IC)-FL(IC-1)-FLA(IC-1))*RDX(IC)
        ROER(IC)=ROER(IC)-DT(IC)*(FLH(IC)-FLH(IC-1))*RDX(IC)
  256   CONTINUE
  260 CONTINUE
C---------------------------------------------------------------DIFFX
  270 IF(IELC.EQ.0 .OR. IELC.EQ.1) GO TO 280
      DO 279 ISP=1,IELC-1
        DO 273 IC=IC1U,IC2
        FRMC=SPDN(IC,ISP)/RON(IC)
        FRMR=SPDN(IC+1,ISP)/RON(IC+1)
        FRM=(DX(IC)*FRMR+DX(IC+1)*FRMC)*HRDXCR(IC)
        FL(IC)=FRM*(FLS(IC)+FLAS(IC))
  273   CONTINUE
        DO 276 IC=IC1,IC2
        SPDR(IC,ISP)=SPDR(IC,ISP)-DT(IC)*(FL(IC)-FL(IC-1))*RDX(IC)
  276 CONTINUE
  279 CONTINUE
C---------------------------------------------------------------DIFFX
  280 IF(IELC.EQ.NSP) GO TO 300
      DO 290 ISP=IELC+1,NSP
        DO 283 IC=IC1U,IC2
        FRMC=SPDN(IC,ISP)/RON(IC)
        FRMR=SPDN(IC+1,ISP)/RON(IC+1)
        FRM=(DX(IC)*FRMR+DX(IC+1)*FRMC)*HRDXCR(IC)
        FL(IC)=FRM*(FLS(IC)+FLAS(IC))
  283   CONTINUE
        DO 286 IC=IC1,IC2
        SPDR(IC,ISP)=SPDR(IC,ISP)-DT(IC)*(FL(IC)-FL(IC-1))*RDX(IC)
  286 CONTINUE
  290 CONTINUE
C-------------------------------------- ELECTRON DIFFUSION -----DIFFX
  300 IF(IELC.EQ.0) GO TO 500
      DO 410 IC=IC1U,IC2
      FRMC=SPDN(IC,IELC)/RON(IC)
      FRMR=SPDN(IC+1,IELC)/RON(IC+1)
      FRM=(DX(IC)*FRMR+DX(IC+1)*FRMC)*HRDXCR(IC)
      FL(IC)=FRM*(FLS(IC)+FLAS(IC))+(FLQ(IC)-FLAQ(IC))/QOM(IELC)
C------------ ENTHALPY DIFFUSION
      WT(IC)=(DX(IC)*TE(IC+1)+DX(IC+1)*TE(IC))*HRDXCR(IC)
      TB=HUNDTH*WT(IC)
      IT=INT(TB)
      IT=MIN(IT,NIT1)
      FR=TB-DBLE(IT)
      HSUBK=(ONE-FR)*EK(IT+1,IELC)+FR*EK(IT+2,IELC)
     1      +RGAS*WT(IC)*RMW(IELC)
      FLH(IC)=HSUBK*(FLQ(IC)-FLAQ(IC))/QOM(IELC)
  410 CONTINUE
      DO 420 IC=IC1,IC2
      SPDR(IC,IELC)=SPDR(IC,IELC)-DT(IC)*(FL(IC)-FL(IC-1))*RDX(IC)
      ROER(IC)=ROER(IC)-DT(IC)*(FLH(IC)-FLH(IC-1))*RDX(IC)
  420 CONTINUE
      IF(NLTE.NE.0) THEN
        DO 430 IC=IC1,IC2
        ROEER(IC)=ROEER(IC)-DT(IC)*(FLH(IC)-FLH(IC-1))*RDX(IC)
  430   CONTINUE
      END IF
C============================================ SWIRL STRESS =====DIFFX
  500 IF(ISWIRL.EQ.0) GO TO 590
      DO 510 IC=IC1U,IC2
      DWOXDX=(WN(IC+1)/XC(IC+1)-WN(IC)/XC(IC))*TWO*HRDXCR(IC)
      AMU1=(DX(IC+1)+DX(IC))*VISC(IC)*VISC(IC+1)
      AMU2=DX(IC+1)*VISC(IC)+DX(IC)*VISC(IC+1)
      AMU=AMU1/AMU2
      XR=XC(IC)+HALF*DX(IC)
      FL(IC)=AMU*XR*XR*DWOXDX*RR(IC)
      DISS(IC)=HALF*DT(IC)*DWOXDX*FL(IC)
  510 CONTINUE
C----------------- VISCOUS DISSIPATION DUE TO SWIRL STRESS -----DIFFX
      IF(ITURB.EQ.0) THEN
        DO 520 IC=IC1,IC2
        RORW(IC)=RORW(IC)+DT(IC)*(FL(IC)-FL(IC-1))*RDX(IC)/XC(IC)
        ROER(IC)=ROER(IC)+DISS(IC)+DISS(IC-1)
  520   CONTINUE
      ELSEIF(ITURB.EQ.1) THEN
        DO 530 IC=IC1,IC2
        RORW(IC)=RORW(IC)+DT(IC)*(FL(IC)-FL(IC-1))*RDX(IC)/XC(IC)
        ROTKER(IC)=ROTKER(IC)+DISS(IC)+DISS(IC-1)
  530   CONTINUE
      ELSEIF(ITURB.EQ.2) THEN
        DO 540 IC=IC1,IC2
        RORW(IC)=RORW(IC)+DT(IC)*(FL(IC)-FL(IC-1))*RDX(IC)/XC(IC)
        ROTKER(IC)=ROTKER(IC)+DISS(IC)+DISS(IC-1)
        ROEPSR(IC)=ROEPSR(IC)
     1             +(DISS(IC)+DISS(IC-1))*CE1*EPSN(IC)/TKEN(IC)
  540   CONTINUE
      END IF
  590 CONTINUE
C============================ FOR TURBULENT KINETIC ENERGY =====DIFFX
      IF(ITURB.NE.0) THEN
        DO 600 IC=IC1U,IC2
        DKDX=(TKEN(IC+1)-TKEN(IC))*TWO*HRDXCR(IC)
        AMU1=(DX(IC+1)+DX(IC))*VISC(IC)*VISC(IC+1)
        AMU2=DX(IC+1)*VISC(IC)+DX(IC)*VISC(IC+1)
        AMU=AMU1/AMU2
        FL(IC)=-AMU*DKDX*RR(IC)
  600   CONTINUE
        DO 610 IC=IC1,IC2
        ROTKER(IC)=ROTKER(IC)-DT(IC)*(FL(IC)-FL(IC-1))*RDX(IC)
  610   CONTINUE
      END IF
C------------- FOR DISSIPATION OF TURBULENT KINETIC ENERGY -----DIFFX
      IF(ITURB.EQ.2) THEN
        DO 700 IC=IC1U,IC2
        DEDX=(EPSN(IC+1)-EPSN(IC))*TWO*HRDXCR(IC)
        AMU1=(DX(IC+1)+DX(IC))*(VISC(IC)+(RPRE-ONE)*VIST(IC))
     1       *(VISC(IC+1)+(RPRE-ONE)*VIST(IC+1))
        AMU2=DX(IC+1)*(VISC(IC)+(RPRE-ONE)*VIST(IC))
     1       +DX(IC)*(VISC(IC+1)+(RPRE-ONE)*VIST(IC+1))
        AMU=AMU1/AMU2
        FL(IC)=-AMU*DEDX*RR(IC)
  700   CONTINUE
        DO 710 IC=IC1,IC2
        ROEPSR(IC)=ROEPSR(IC)-DT(IC)*(FL(IC)-FL(IC-1))*RDX(IC)
  710   CONTINUE
      END IF
C===============================================================DIFFX
      RETURN
      END
