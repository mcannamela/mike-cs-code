*DECK XPORT
      SUBROUTINE XPORT
C===============================================================XPORT
C     THIS ROUTINE
C     PROVIDES VISCOSITIY, CONDUCTIVITY, AND DIFFUSION COEFFICIENTS
C     LOCATED AT CELL CENTER.
C----------
C     CALLED BY LAVA
C===============================================================XPORT
      INCLUDE 'COML.h'
C----------
      DOUBLE PRECISION TB,FR,FR1,DCOE1,RMKBAR,PSUBE
C----------
C----- K-EPSILON CORRECTION
      INTEGER ICC,ICS
      DOUBLE PRECISION XXX,XX1,XX2,VCL1,VCL2,DELV1,DELV2,DVCLDY,
     1                 VV,VD1,VD2,XDEL,FCOR,CE2COR
c-----
      double precision DAAi,DAAe,DAiAe,wgha,wghi,wghe,wght,
     1     csuba,csubi,csube,xtotal,xsuba,xsubi,xsube,cbase
c-----
      double precision Airtio,cArb,cHb,cNb,cArp,cHp,cAirp,rctot
c     data Airtio /3.761904762d0/
C=================================== LAMINAR DIFFUSIVITIES =====XPORT
C     PROPERTIES FOR ARGON PLASMA. FOR MIXTURE, WE NEED TO COME UP
C     WITH SOME KIND OF MIXTURE RULE FOR 2-T PLASMAS (KRUGER'S BOOK?)
C     COND IS TOTAL THERMAL CONDUCTIVITY (EITHER TRANSLATIONAL OR
C     REACTIONAL INCLUDED) FOR NLTE=0, WHILE COND IS HEAVY PARTICLE
C     THERMAL CONDUCTIVITY FOR NLTE>0.
C---------------------------------------------------------------XPORT
C     RADIATION LOSS DEFINED HERE SHOULD NOT BE OVERLAPPED WITH THE
C     RADIATION LOSS USED IN CHEMISTRY ROUTINES. I.E., USER SHOULD BE
C     CAREFUL NOT TO COUNT SAME RADIATION LOSS TWICE.
C---------------------------------------------------------------XPORT
      DO 100 IC=1,NXYZT
      TB=HUNDTH*TEMP(IC)
      IT=INT(TB)
      IT=MIN(IT,NIT1)
      FR=TB-DBLE(IT)
      FR1=ONE-FR
C---------------------------------------------------------------XPORT
c---- using argon, hydrogen, and air properties and mole fractions
c---- ignore the differences of N-O ratio than air
      cArb=spd(ic,1)*rmw(1)+spd(ic,2)*rmw(2)
      cHb=spd(ic,3)*rmw(3)+half*(spd(ic,4)*rmw(4)+spd(ic,5)*rmw(5))
      cNb=spd(ic,6)*rmw(6)+spd(ic,7)*rmw(7)
     1    +half*(spd(ic,8)*rmw(8)+spd(ic,9)*rmw(9))
      rctot=one/(cArb+cHb+cNb)
      cArp=cArb*rctot
      cHp=cHb*rctot
      cAirp=cNb*rctot
      VISC(IC)=cArp*(FR1*VIS(IT+1,1)+FR*VIS(IT+2,1))
     1         +cHp*(FR1*VIS(IT+1,2)+FR*VIS(IT+2,2))
     2       +cAirp*(FR1*VIS(IT+1,3)+FR*VIS(IT+2,3))
      COND(IC)=cArp*(FR1*CND(IT+1,1)+FR*CND(IT+2,1))
     1         +cHp*(FR1*CND(IT+1,2)+FR*CND(IT+2,2))
     2       +cAirp*(FR1*CND(IT+1,3)+FR*CND(IT+2,3))
      RADLSS(IC)=cArp*(FR1*SSUBR(IT+1,1)+FR*SSUBR(IT+2,1))
     1           +cHp*(FR1*SSUBR(IT+1,2)+FR*SSUBR(IT+2,2))
     2         +cAirp*(FR1*SSUBR(IT+1,3)+FR*SSUBR(IT+2,3))
C---------------------------------------------------------------XPORT
      RADLSS(IC)=RADLSS(IC)*P(IC)/P1ATM
C---------------------------------------------------------------XPORT
      VIST(IC)=ZERO
  100 CONTINUE
C----------ELECTRON THERMAL CONDUCTIVITY
      IF(NLTE.NE.0) THEN
        DO 105 IC=1,NXYZT
        TB=HUNDTH*TE(IC)
        IT=INT(TB)
        IT=MIN(IT,NIT1)
        FR=TB-DBLE(IT)
        FR1=ONE-FR
C----------
        CONDE(IC)=FR1*CNDE(IT+1,1)+FR*CNDE(IT+2,1)
        RADLSS(IC)=FR1*SSUBR(IT+1,1)+FR*SSUBR(IT+2,1)
        RADLSS(IC)=RADLSS(IC)*P(IC)/P1ATM
        IF(NLTE.NE.1) RADLSS(IC)=ZERO
C---------------------------------------------------------------XPORT
C     FOLLOWING IS NECESSARY TO PREVENT THE DIFFUSION TIME STEP LIMIT
C     FROM BECOMING TOO RESTRICTIVE. WHEN PSUBE IS SUFFICIENTLY
C     SMALL, WE WILL SET CONDE=SMALL.
C     HOW SMALL IS SUFFICIENT? -- TO BE DETERMINED.
C---------------------------------------------------------------XPORT
        PSUBE=SPDN(IC,IELC)*RGAS*RMW(IELC)*TE(IC)
        IF(PSUBE.LT.THOUTH*PAMB) CONDE(IC)=SMALL
  105   CONTINUE
      END IF
C======================================= SPECIES DIFFUSION =====XPORT
      IF(NODIFF.NE.1) THEN
C---------------------------------------------------------------XPORT
C     SET DIFFUSION COEFFICEINT FOR ELECTRON ZERO
C     SET DIFFUSION COEFFICEINT FOR CONDENSED PHASE ZERO (AUTOMATIC
C          APACHE FORMULA DUE TO LARGE MOLECULAR WEIGHT)
C---------------------------------------------------------------XPORT
C     DIFFUSION COEFFICEINT IS CALCULATED USING HSU'S PROGRAM FOR
C     ARGON PLASMA. WHEN USING THIS, APACHE FORMULA PART SHOULD BE
C     DEACTIVATED. SAME LINE LABEL IS USED TO PREVENT DUPLICATE.
C---------------------------------------------------------------XPORT
C       DO 110 IC=1,NXYZT
C       TB=HUNDTH*TEMP(IC)
C       IT=INT(TB)
C       IT=MIN(IT,NIT1)
C       FR=TB-DBLE(IT)
C       FR1=ONE-FR
C       DAAi=(FR1*DCOFF(IT+1,1)+FR*DCOFF(IT+2,1))*P1ATM/P(IC)
C       DCOEF(IC,IELC)=ZERO
c--- for Hoffert-Lien or LTE
C       IF(NLTE.EQ.0) THEN
C         DCOEF(IC,1)=DAAi
C         DCOEF(IC,2)=DAAi
C       ELSEIF(NLTE.EQ.1) THEN
C         DCOEF(IC,1)=DAAi
C         DCOEF(IC,2)=DAAi
c--- for Braun and Kunc
C       ELSE
C         DAAe=0.44d0*DAAi
c         DAiAe=0.59d0*DAAi
c         wgha=spd(ic,1)*sqrt(rmw(1))
c         wghi=spd(ic,2)*sqrt(rmw(2))
c         wghe=spd(ic,4)*sqrt(rmw(4))
c         wght=wgha+wghi+wghe
c         csuba=spd(ic,1)*rmw(1)
c         csubi=spd(ic,2)*rmw(2)
c         csube=spd(ic,4)*rmw(4)
c         xtotal=csuba+csubi+csube
c         xsuba=csuba/xtotal
c         xsubi=csubi/xtotal
c         xsube=csube/xtotal
c         if(xsubi.le.small .and. xsube.le.small) then
c           dcoef(ic,1)=zero
c         else
c           dcoef(ic,1)=(one-wgha/wght)/(xsubi/DAAi+xsube/DAAe)
c         endif
c         dcoef(ic,2)=(one-wghi/wght)/(xsuba/DAAi+xsube/DAiAe)
c         dcoef(ic,4)=(one-wghe/wght)/(xsuba/DAAe+xsubi/DAiAe)
C       END IF
C 110   CONTINUE
C---------------------------------------------------------------XPORT
C     DIFFUSION COEFFICIENTS ARE CALCULATED USING THE FORMULA IN
C     APACHE REPORT P. 109, DCOEF=DCOEF+DTURB
C     ELECTRON DIFFUSION COEFFICIENT IS SET TO BE ZERO HERE,
C     FOR LATER USE IN TIMSTP. AFTER DTURB IS ADDED, IT IS NOT
C     NECESSARY TO SET DCOEF OF IELC TO BE ZERO.
C---------------------------------------------------------------XPORT
C     ABOVE RELATIONS CONTAIN THERMODYNAMIC REALTIONS. NEED TO BE
C     FIXED IN ORDER TO USE FOR 2-T CASES.
C---------------------------------------------------------------XPORT
C     SET DIFFUSION COEFFICIENTS OF CONDENSED SPECIES ZERO
C     EVEN THOUGH GAS DIFFUSION COEFFICIENTS ARE NOT ZERO, BINARY
C     DIFFUSION APPROXIMATION WILL TAKE CARE OF MASS BALANCE
C---------------------------------------------------------------XPORT
      DO 120 ISP=1,NSP
      DO 110 IC=1,NXYZT
C---------------------------------------------------------------XPORT
C     SMALLEST REASONABLE VALUE OF MKBAR IS ONE
C---------------------------------------------------------------XPORT
        RMKBAR=(P(IC)/(RGAS*TEMP(IC))-SPD(IC,ISP)*RMW(ISP))
     1          /MAX((RO(IC)-SPD(IC,ISP)),SMALL)
        RMKBAR=MAX(RMKBAR,ZERO)
        RMKBAR=MIN(RMKBAR,ONE)
        DCOE1=RMW(ISP)+RMKBAR
        DCOEF(IC,ISP)=DBLE(IGAS(ISP))*ETA0*TEMP(IC)
     1                               *SQRT(DCOE1*TEMP(IC))/P(IC)
        DCOEF(IC,IELC)=ZERO
  110 CONTINUE
  120 CONTINUE
C---------------------------------------------------------------XPORT
      ELSE
        DO 160 ISP=1,NSP
        DO 150 IC=1,NXYZT
          DCOEF(IC,ISP)=ZERO
  150   CONTINUE
  160   CONTINUE
      END IF
C=================================== TURBULENT VISCOSITIES =====XPORT
      IF (ITURB.EQ.0) GO TO 400
C----------------------------------------------- SGS MODEL -----XPORT
      IF (ITURB.EQ.1) THEN
        DO 200 IC=1,NXYZT
        VIST(IC)=SGSA*SGSL(IC)*RO(IC)*SQRT(TKE(IC))
  200   CONTINUE
C----------------------------------------- K-EPSILON MODEL -----XPORT
      ELSEIF (ITURB.EQ.2) THEN
      IF(KECOR.EQ.1) THEN
C---------------------------------------------------------------XPORT
C     CORRECTION OF K-EPSILON FOR AXISYMMETRIC JET (LAUNDER ET. AL.)
C---------------------------------------------------------------XPORT
      DO 217 J=J1,J2
      ICC=(J-1)*NXT+2
      XXX=ONE/(XC(ICC+1)**2-XC(ICC)**2)
      XX2=XC(ICC+1)**2*XXX
      XX1=XC(ICC)**2*XXX
C----- CENTERLINE VELOCITIES ARE DETERMINED USING PARABOLIC FIT
      VCL1=V(ICC-NXT)*XX2-V(ICC+1-NXT)*XX1
      VCL2=V(ICC)*XX2-V(ICC+1)*XX1
      IF(VCL1.GE.ONE .AND.VCL2.GE.ONE) THEN
C----- MIXING LENGTH IS SET TO BE 1% OF UCL
      DELV1=AKEML1*(VCL1+VCL2)
      DELV2=AKEML2*(VCL1+VCL2)
      DVCLDY=(VCL2-VCL1)*RDY(ICC)
C----- FIND MIXING LENGTH XDEL
        DO 212 I=I1,I2
        IC=(J-1)*NXT+I
        VV=(V(IC-NXT)+V(IC))*HALF
        ICS=IC
        IF(VV.LE.DELV2) THEN
          GO TO 213
        END IF
  212   CONTINUE
  213 CONTINUE
      VD1=(V(ICS-1)+V(ICS-1-NXT))*HALF
      VD2=(V(ICS)+V(ICS-NXT))*HALF
        IF(VD1.EQ.VD2) THEN
          XDEL=ZERO
        ELSE
          XDEL=XC(ICS-1)+(XC(ICS)-XC(ICS-1))*(VD1-DELV2)/(VD1-VD2)
        END IF
C----- CALCULATE CORRECTION FACTOR
      FCOR=(ABS(XDEL*(DVCLDY-ABS(DVCLDY))/DELV1))**0.2
      ELSE
        FCOR=ZERO
      END IF
C----- CALCULATE TURBULENT VISCOSITY
        DO 215 I=I1,I2
        IC=(J-1)*NXT+I
        VIST(IC)=(CMU-CMUCOR*FCOR)*TKE(IC)*TKE(IC)*RO(IC)/EPS(IC)
  215   CONTINUE
  217 CONTINUE
C----------
      ELSE
        DO 210 IC=1,NXYZT
        VIST(IC)=CMU*TKE(IC)*TKE(IC)*RO(IC)/EPS(IC)
  210   CONTINUE
      END IF
      END IF
C================================= TURBULENT DIFFUSIVITIES =====XPORT
      IF(NLTE.EQ.0) THEN
      
		DO 310 IC=1,NXYZT
        VISC(IC)=VISC(IC)+VIST(IC)
        COND(IC)=COND(IC)+VIST(IC)*RPRT*GAMMA(IC)*P(IC)
     1                    /((GAMMA(IC)-ONE)*RO(IC)*TEMP(IC))
  310   CONTINUE
C--------------------------------------------- FOR NON-LTE -----XPORT
      ELSE
        DO 360 IC=1,NXYZT
        VISC(IC)=VISC(IC)+VIST(IC)
        PSUBE=SPDN(IC,IELC)*RGAS*RMW(IELC)*TE(IC)
        COND(IC)=COND(IC)+VIST(IC)*RPRT*GAMMA(IC)*(P(IC)-PSUBE)
     1                    /((GAMMA(IC)-ONE)*RO(IC)*TEMP(IC))
        CONDE(IC)=CONDE(IC)+VIST(IC)*RPRT*FIVHAF*PSUBE
     1                      /(RO(IC)*TE(IC))
  360   CONTINUE
      END IF
      DO 380 ISP=1,NSP
      DO 370 IC=1,NXYZT
        DCOEF(IC,ISP)=DCOEF(IC,ISP)+RSCT*VIST(IC)/RO(IC)
  370 CONTINUE
  380 CONTINUE
C===============================================================XPORT
  400 CONTINUE
C===============================================================XPORT
C     ALONG THE BOUNDARY. FOR BLOCKAGE, SET PROPERTIES ACCORDINGLY AT
C     THE BLOCKAGE BOUNDARY.
C--------------------------- FOR LEFT AND RIGHT BOUNDARIES -----XPORT
      IF(NX.GT.1) THEN
        DO 410 ICLR=1,NLR
        ICL=1+(ICLR-1)*NXT
        ICR=ICL+NX+1
        VISC(ICL)=VISC(ICL+1)
        VIST(ICL)=VIST(ICL+1)
        COND(ICL)=COND(ICL+1)
        CONDE(ICL)=CONDE(ICL+1)
        VISC(ICR)=VISC(ICR-1)
        VIST(ICR)=VIST(ICR-1)
        COND(ICR)=COND(ICR-1)
        CONDE(ICR)=CONDE(ICR-1)
  410   CONTINUE
        DO 430 ISP=1,NSP
          DO 420 ICLR=1,NLR
          ICL=1+(ICLR-1)*NXT
          ICR=ICL+NX+1
          DCOEF(ICL,ISP)=DCOEF(ICL+1,ISP)
          DCOEF(ICR,ISP)=DCOEF(ICR-1,ISP)
  420     CONTINUE
  430   CONTINUE
      END IF
C----------------------- FOR DERRIERE AND FRONT BOUNDARIES -----XPORT
      IF(NY.GT.1) THEN
        DO 470 K=1,NZT
        ICK=(K-1)*NXYT
          DO 460 I=1,NXT
          ICD=I+ICK
          ICF=ICD+NXYT-NXT
          VISC(ICD)=VISC(ICD+NXT)
          VIST(ICD)=VIST(ICD+NXT)
          COND(ICD)=COND(ICD+NXT)
          CONDE(ICD)=CONDE(ICD+NXT)
          VISC(ICF)=VISC(ICF-NXT)
          VIST(ICF)=VIST(ICF-NXT)
          COND(ICF)=COND(ICF-NXT)
          CONDE(ICF)=CONDE(ICF-NXT)
  460     CONTINUE
  470   CONTINUE
        DO 500 ISP=1,NSP
          DO 490 K=1,NZT
          ICK=(K-1)*NXYT
            DO 480 I=1,NXT
            ICD=I+ICK
            ICF=ICD+NXYT-NXT
            DCOEF(ICD,ISP)=DCOEF(ICD+NXT,ISP)
            DCOEF(ICF,ISP)=DCOEF(ICF-NXT,ISP)
  480       CONTINUE
  490     CONTINUE
  500   CONTINUE
      END IF
C--------------------------- FOR TOP AND BOTTOM BOUNDARIES -----XPORT
      IF(NZ.GT.1) THEN
        DO 510 ICBT=1,NBT
        ICB=ICBT
        ICT=ICBT+NXYZT-NXYT
        VISC(ICB)=VISC(ICB+NXYT)
        VIST(ICB)=VIST(ICB+NXYT)
        COND(ICB)=COND(ICB+NXYT)
        CONDE(ICB)=CONDE(ICB+NXYT)
        VISC(ICT)=VISC(ICT-NXYT)
        VIST(ICT)=VIST(ICT-NXYT)
        COND(ICT)=COND(ICT-NXYT)
        CONDE(ICT)=CONDE(ICT-NXYT)
  510   CONTINUE
        DO 530 ISP=1,NSP
          DO 520 ICBT=1,NBT
          ICB=ICBT
          ICT=ICBT+NXYZT-NXYT
          DCOEF(ICB,ISP)=DCOEF(ICB+NXYT,ISP)
          DCOEF(ICT,ISP)=DCOEF(ICT-NXYT,ISP)
  520     CONTINUE
  530   CONTINUE
      END IF
C===============================================================XPORT
      RETURN
      END
