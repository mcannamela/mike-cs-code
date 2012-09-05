*DECK PDVETC
      SUBROUTINE PDVETC
C==============================================================PDVETC
C     THIS ROUTINE UPDATES ENERGY DUE TO PDV WORK, RADIATION, ETC.
C----------
C     CALLED BY LAVA
C----------
C     INPUT: ROER,ROEN,DT,DIV,RC,OTHER SOURCE/SINK TERMS
C     OUTPUT: ROER
C==============================================================PDVETC
      INCLUDE 'COML.h'
C----------
      DOUBLE PRECISION PEDIVA,PEDIVB,ALFEI,ALFEA,ALFEH,SEA,SS
C----------
C----- K-EPSILON CORRECTION
      INTEGER ICC,ICS
      DOUBLE PRECISION XXX,XX1,XX2,VCL1,VCL2,DELV1,DELV2,DVCLDY,
     1                 VV,VD1,VD2,XDEL,FCOR,CE2COR
C-------------------------------------- FOR NON-LTE MODEL -----PDVETC
      DOUBLE PRECISION COLF,TEIC,TEIC2,DTIC,QEA,QEAE,QEI,SROQOM,ALF
      DOUBLE PRECISION CEA1,CEA2,CEA3,CEA4,CEA5,CEI1,CEI2
      DOUBLE PRECISION DRATI,EQK,RADLSE,ALFOR,DENOM
      DATA CEI1,CEI2 /1.95D-6,1.53D8/
C----------FOR HOFFERT-LIEN MODEL
      DATA CEA1,CEA2 /-3.5D-17,7.75D-21/
      DATA CEA3,CEA4,CEA5 /3.9D-17,5.51D-21,5.95D-25/
C--------------------------------------------------------------PDVETC
C     WSDT IS FROM SUBROUTINE PINT
C     WSDT SHOULD BE NEGATIVE
C         (SEE P. J. O'ROURKE, J. COMP. PHYS. 83, 345 (1989).)
C========================================== FOR SGS MODEL =====PDVETC
      IF(ITURB.EQ.1) THEN
        DO 100 IC=IC1,IC2
        WSDT(IC)=MIN(WSDT(IC),ZERO)
        ROTKER(IC)=(ROTKER(IC)-TWOO3*DIV(IC)*RC(IC)*RON(IC)
     1    *TKEN(IC)*DT(IC))/(ONE+DT(IC)*SGSD*SQRT(TKEN(IC))/SGSL(IC)
     2                       -WSDT(IC)*RRC(IC)/(RON(IC)*TKEN(IC)))
  100   CONTINUE
      END IF
C==================================== FOR K-EPSILON MODEL =====PDVETC
      IF(ITURB.EQ.2) THEN
      IF(KECOR.EQ.1) THEN
C--------------------------------------------------------------PDVETC
C     CORRECTION OF K-EPSILON FOR AXISYMMETRIC JET (LAUNDER ET. AL.)
C--------------------------------------------------------------PDVETC
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
C--------------------------------------------------------------PDVETC
        CE2COR=CE2-CE2CR*FCOR
        DO 215 I=I1,I2
        IC=(J-1)*NXT+I
         WSDT(IC)=MIN(WSDT(IC),ZERO)
         ROTKER(IC)=(ROTKER(IC)-TWOO3*DIV(IC)*RC(IC)*RON(IC)
     1              *TKEN(IC)*DT(IC))*TKEN(IC)
     2           /(TKEN(IC)+DT(IC)*EPSN(IC)-WSDT(IC)*RRC(IC)/RON(IC))
      ROEPSR(IC)=(ROEPSR(IC)-(TWOO3*CE1-CE3)*DIV(IC)*RC(IC)*RON(IC)
     1           *EPSN(IC)*DT(IC))*TKEN(IC)
     2 /(TKEN(IC)+DT(IC)*CE2COR*EPSN(IC)-CES*WSDT(IC)*RRC(IC)/RON(IC))
  215   CONTINUE
  217 CONTINUE
      ELSE
        DO 110 IC=IC1,IC2
         WSDT(IC)=MIN(WSDT(IC),ZERO)
         ROTKER(IC)=(ROTKER(IC)-TWOO3*DIV(IC)*RC(IC)*RON(IC)
     1              *TKEN(IC)*DT(IC))*TKEN(IC)
     2           /(TKEN(IC)+DT(IC)*EPSN(IC)-WSDT(IC)*RRC(IC)/RON(IC))
      ROEPSR(IC)=(ROEPSR(IC)-(TWOO3*CE1-CE3)*DIV(IC)*RC(IC)*RON(IC)
     1           *EPSN(IC)*DT(IC))*TKEN(IC)
     2   /(TKEN(IC)+DT(IC)*CE2*EPSN(IC)-CES*WSDT(IC)*RRC(IC)/RON(IC))
  110   CONTINUE
      END IF
      END IF
C======================== ENERGY SOURCE DUE TO TURBULENCE =====PDVETC
C     ENERGY SOURCE DUE TO TURBULENCE DOES NOT HAVE ANY CONTRIBUTION
C     TO ELECTRON ENERGY EQUATION, SINCE THIS SOURCE IS RELATED TO
C     VISCOUS DISSIPATION.
C--------------------------------------------------------------PDVETC
      IF(ITURB.EQ.1) THEN
        DO 300 IC=IC1,IC2
        ROER(IC)=ROER(IC)+DT(IC)*ROTKER(IC)*SGSD*SQRT(TKEN(IC))/SGSL(IC)
  300   CONTINUE
      END IF
C----------
      IF(ITURB.EQ.2) THEN
        DO 310 IC=IC1,IC2
        ROER(IC)=ROER(IC)+DT(IC)*ROTKER(IC)*EPSN(IC)/TKEN(IC)
  310   CONTINUE
      END IF
C================================================ SOURCES =====PDVETC
      IF(INCOMP.EQ.1) GO TO 90
C==============================================================PDVETC
C     FOR LTE, ENERGY EQUATION WILL BE UPDATED EXPLICITLY
C     RADLSS = RADIATION HEAT LOSS OF FLUID > ZERO
C--------------------------------------------------------------PDVETC
      IF(NLTE.EQ.0) THEN
        DO 400 IC=IC1,IC2
        ROER(IC)=ROER(IC)-DT(IC)*(P(IC)*DIV(IC)+RADLSS(IC))*RC(IC)
  400   CONTINUE
C==============================================================PDVETC
C     FOR NON-LTE, BOTH TOTAL ENERGY EQUATION AND ELECTRON ENERGY
C     EQUATION WILL BE UPDATED LINEARLY IMPLICITLY FOR ALL SOURCE
C     TERMS AT ONCE.  THIS LOGIC WITH CHEMEQ MIGHT NEED FUTURE
C     IMPROVEMENTS.
C     THIS PART BOTH CONTAINS EXPLICIT AND IMPLICIT UPDATE OF THE
C     ENERGY, SHOULD BE PLACED AFTER ALL THE EXPLICIT CHANGE AND
C     BEFORE THE IMPLICIT CHANGE.
C--------------------------------------------------------------PDVETC
      ELSE
C--------------------------------------- FOR TOTAL ENERGY -----PDVETC
        DO 500 IC=IC1,IC2
        SS=SIGN(ONE,DIV(IC))
        ROER(IC)=(ROER(IC)-DT(IC)*HALF*(ONE-SS)*P(IC)*DIV(IC)*RC(IC))
     1/(ONE+DT(IC)*(RADLSS(IC)+HALF*(ONE+SS)*P(IC)*DIV(IC))/ROEN(IC))
  500   CONTINUE
C------------------------------------ FOR ELECTRON ENERGY -----PDVETC
        COLF=THREE*AVOGAD*RGAS*SQRT(EIGHT*RGAS*RMW(IELC)/PI)
        DO 600 IC=IC1,IC2
        SS=SIGN(ONE,DIV(IC))
        PEDIVA=ONEO3*(ONE-SS)*DIV(IC)*RC(IC)*ROEEN(IC)
        PEDIVB=ONEO3*(ONE+SS)*DIV(IC)
C--------------------------------------------------------------PDVETC
C     RADIATION LOSS DUE TO FREE ELECTRONS ARE ASSUMED TO BE ZERO FOR
C     EXCITED-STATE KINETICS.  FOR PLAIN TWO-TEMPERATURE MODEL LIKE
C     HOFFERT-LIEN MODEL, IT'S ASSUMED TO BE 30 % OF THE RADIATION
C     PROPORTIONAL TO THE RATIO OF ELECTRON DENSITY TO THE
C     EQUILIBRIUM ELECTRON DENSITY BASED ON ELECTRON TEMPERATURE.
C--------------------------------------------------------------PDVETC
        IF(NLTE.EQ.1) THEN
          IF(TE(IC).GT.8000.0D0) THEN
            EQK=2.9D16*TE(IC)*SQRT(TE(IC))*EXP(-183100.0D0/TE(IC))
            DRATI=SPDR(IC,IELC)*RMW(IELC)*AVOGAD*RRC(IC)
     1            /(SQRT((EQK+PAMB/BOLTZ/TE(IC))*EQK)-EQK)
            RADLSE=RADLSS(IC)*0.3D0*DRATI/ROEEN(IC)
          ELSE
            RADLSE=ZERO
          END IF
        ELSE
          RADLSE=ZERO
        END IF
C---------------- ELECTRON-HEAVY PARTICLE ENERGY EXCHANGE -----PDVETC
C     RELATIONS FROM
C     M. I. HOFFERT AND H. LIEN, PHYS. FLUIDS 10, P1769 (1967).
C     C. G. BRAUN AND J. A. KUNC, PHYS. FLUIDS 30, P499 (1987).
C     Q ARE IN CGS (CM**2)
C--------------------------------------------------------------PDVETC
        TEIC=TE(IC)
        TEIC2=TE(IC)*TEIC
        DTIC=DT(IC)
C--------------------------------------------------------------PDVETC
C     ELECTRON-ION COULOMB COLLISION (BY LANDAU)
C     THIS PART IS CODED JUST FOR SINGLY IONIZED IONS ONLY.
C--------------------------------------------------------------PDVETC
        IF(SPDN(IC,IELC).LE.SMALL) THEN
          QEI=ZERO
        ELSE
          QEI=(CEI1/TEIC2)
     1     *LOG(ONE+CEI2*TEIC*TEIC2/(AVOGAD*RMW(IELC)*SPDN(IC,IELC)))
        END IF
C----------ELECTRON-ARGON ATOM FROM HOFFERT-LIEN (BY JAFFRIN)
        IF(NLTE.EQ.1) THEN
          IF(TEIC.GT.10000.0D0) THEN
            QEA=CEA1+CEA2*TEIC
          ELSE
            QEA=CEA3+(CEA4+CEA5*TEIC)*TEIC
          END IF
          QEAE=ZERO
C--------------------------------------------------------------PDVETC
C     ELECTRON-ARGON ATOM FROM BRAUN-KUNC
C     (BY FROST-PHELPS, FLETCHER-BURCH)
C--------------------------------------------------------------PDVETC
        ELSEIF(NLTE.EQ.2) THEN
          IF(TEIC.le.2732.0d0) THEN
           QEA=1.229D-13/TEIC
          ELSEIF(TEIC.GT.2732.0D0 .AND. TEIC.LE. 91064.0D0) THEN
           QEA=1.647D-20*TEIC
          ELSEIF(TEIC.GT.91064.0D0) THEN
           QEA=1.366D-10/TEIC
          END IF
C--------------------------------------------------------------PDVETC
C    MORE ELECTRON-NEUTRAL COLLISION CROSS SECTIONS CAN BE CODED HERE
C--------------------------------------------------------------PDVETC
C     ELECTRON-EXCITED STATE ARGON CROSS SECTION CAN BE DETERMINED
C     FROM THE SIZE RATIO OF AR AND AR* BASED ON THE FOLLOWING
C     ASSUMPTIONS AND SHOULD BE CODED HERE.
C     1. AR AND AR* HAVE SIMILAR POTETIAL WITH SIZE DIFFERENCE
C     2. AR* HAS SIZE OF POTASIUM (K)
C         QEAE=QEA*2.25D0
C--------------------------------------------------------------PDVETC
        END IF
C--------------------------------------------------------------PDVETC
C     SROQOM NEEDS TO BE CODED BY HAND
C--------------------------------------------------------------PDVETC
        sroqom=Qea*(spdn(ic,1)*rmw(1)*rmw(1)+spdn(ic,2)*rmw(2)*rmw(2))
     1        +Qei*(spdn(ic,3)*rmw(3)*rmw(3)+spdn(ic,7)*rmw(7)*rmw(7)
     2         +spdn(ic,9)*rmw(9)*rmw(9)+spdn(ic,10)*rmw(10)*rmw(10))
C--------------------------------------------------------------PDVETC
        ALF=COLF*SPDN(IC,IELC)*SQRT(TEIC)*SROQOM
        ALFOR=TWOO3*MW(IELC)*COLF*SQRT(TEIC)*SROQOM/RGAS
        DENOM=ONE+(PEDIVB+RADLSE+ALFOR)*DT(IC)
      ROEER(IC)=(ROEER(IC)+(ALF*TEMP(IC)*RC(IC)-PEDIVA)*DT(IC))/DENOM
  600   CONTINUE
      END IF
C==============================================================PDVETC
   90 CONTINUE
      RETURN
      END
