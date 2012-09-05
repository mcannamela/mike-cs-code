*DECK PINTCD
      SUBROUTINE PINT
C================================================================PINT
C     THIS SUBROUTINE SOLVES PARTICLE MOMENTUM AND ENERGY EQUATIONS.
C     FIRST, RELVEL IS CALCULATED AND PROVIDED TO USERPP.
C     THEN, USERPP CALCULATES CD AND HTCO.
C     ACCUMULATE HEAT SOURCE/SINK DUE TO HEAT TRANSFER TO PARTICLES
C     AND MOMENTUM SOURCE/SINK DUE TO MOMENTUM TRANSFER TO PARTICLES.
C----------
C     CALLED BY LAVA
C     CALLS USERPP FOR CD(IP) AND HTCO(IP), CALLs PTCHEATING
C----------
C     SOME MODIFICATION WERE MADE TO INCLUDE THE NEW PARTICLE HEATING
C     MODEL BY Y.P. WAN 12-10-97
C================================================================PINT
      INCLUDE 'COML.h'
      INCLUDE 'COMMPTC.h'
C---------- 
      DOUBLE PRECISION FRAN
      EXTERNAL FRAN,USERPP
C----------
      DOUBLE PRECISION CDRVDT,UPSAVE,VPSAVE,WPSAVE,UG,VG,WG,
     1     UF,VF,WF,TRMPX,TRMPY,TRMPZ,SINTH,COSTH
      DOUBLE PRECISION EPSAVE,QPDT,ADTOM,HASET3
      DOUBLE PRECISION TSC1,TSC2,ETA,ETA1
      DOUBLE PRECISION TURVEL,DRAGDT,ATD,EXPMTD,EXPMDT,PALP,RDRGDT,
     1     PALPTE,PBET,PGAM,PLAM,TURTE,TURDIS
      DOUBLE PRECISION UPRIMX,VPRIMY,WPRIMZ,UPRIME,VPRIME,WPRIME,
     1     XPRIME,YPRIME,ZPRIME
      DOUBLE PRECISION DTIPSUM,DTIP
C----------
      INTEGER IND
C ---- CONTROL THE TURBULENCE DISPERSION BY Y.P. WAN 4-22-98
      INTEGER IDISPER
      DATA IDISPER/1/
C================================================================PINT
      DO 100 IP=1,NP
      IC=ICP(IP)
C----------------------------------------------------------------PINT
      UF=HALF*(UN(IC)+UN(IC-1))
      UG=UF
      VF=HALF*(VN(IC)+VN(IC-NXT))
      VG=VF
      IF(ISWIRL.EQ.1) THEN
        WF=WN(IC)
      ELSEIF(NZ.GT.1) THEN
        WF=HALF*(WN(IC)+WN(IC-NXYT))
      ELSE
        WF=ZERO
      END IF
      WG=WF
      IF(I3DP2D.EQ.1) THEN
        COSTH=XP(IP)/MAX(SQRT(XP(IP)**2+ZP(IP)**2),SMALL)
        SINTH=ZP(IP)/MAX(SQRT(XP(IP)**2+ZP(IP)**2),SMALL)
        UG=UF*COSTH-WF*SINTH
        WG=UF*SINTH+WF*COSTH
      END IF
      RELVEL(IP)=SQRT((UG+UTRB(IP)-UP(IP))**2
     1           +(VG+VTRB(IP)-VP(IP))**2+(WG+WTRB(IP)-WP(IP))**2)
C==================== TURBULENT DISPERSION WHEN DT < TSCALE =====PINT
C     TURBULENT DISPERSION OF PARTICLES FROM
C          KIVA1 (FOR TIME SCALE FOR SGS MODEL), KIVA2,
C          AND P. J. O'ROURKE, J. COMP. PHYS. 83, 345 (1989).
C----------------------------------------------------------------PINT
      IF(ITURB.NE.0) THEN
        IF(ITURB.EQ.1) THEN
          TSC1=SGSL(IC)/SQRT(TKE(IC))
          TSC2=SGSL(IC)/RELVEL(IP)
        ELSE
          TSC1=TKE(IC)/(EPS(IC)+SMALL)
          TSC2=CPS*SQRT(TKE(IC))*TKE(IC)
     1         /((EPS(IC)+SMALL)*(RELVEL(IP)+SMALL))
        END IF
        TSCALE(IP)=MIN(TSC1,TSC2)
C ---- SET A SWITCH FOR THE TURBULENCE DISPERSION IDISPER=1, TURN ON
        IF(DT(IC1).LE.TSCALE(IP).AND.IDISPER.EQ.1) THEN
          TURBT(IP)=TURBT(IP)+DT(IC1)
          UPRIMX=ZERO
          VPRIMY=ZERO
          WPRIMZ=ZERO
          IF(TURBT(IP).GT.TSCALE(IP)) THEN
C----------------------------------------------------------------PINT
C     CALCULATE NEW TURBULENT VELOCITIES WHICH PARTICLES SEE.
C----------------------------------------------------------------PINT
            TURBT(IP)=ZERO
            TURVEL=SQRT(FOURO3*TKE(IC))
C----- FOR X
            ETA1=ONE-TWO*FRAN(ZERO)
            ETA=ABS(ETA1)
            IND=INT(ONE+TWENTY*ETA)
            UPRIME=RERF(IND)+(ONE+TWENTY*ETA-DBLE(IND))
     1                      *(RERF(IND+1)-RERF(IND))
            UTRB(IP)=UPRIME*SIGN(ONE,ETA1)*TURVEL
C----- FOR Y
            ETA1=ONE-TWO*FRAN(ZERO)
            ETA=ABS(ETA1)
            IND=INT(ONE+TWENTY*ETA)
            VPRIME=RERF(IND)+(ONE+TWENTY*ETA-DBLE(IND))
     1                      *(RERF(IND+1)-RERF(IND))
            VTRB(IP)=VPRIME*SIGN(ONE,ETA1)*TURVEL
C----- FOR Z
C---------- WTRB(IP)=ZERO FOR NZ=1 OR 2-D PARTICLES (SEE INJECT)
            IF(I3DP2D.EQ.1 .OR. NZ.GT.1) THEN
              ETA1=ONE-TWO*FRAN(ZERO)
              ETA=ABS(ETA1)
              IND=INT(ONE+TWENTY*ETA)
              WPRIME=RERF(IND)+(ONE+TWENTY*ETA-DBLE(IND))
     1                        *(RERF(IND+1)-RERF(IND))
              WTRB(IP)=WPRIME*SIGN(ONE,ETA1)*TURVEL
            END IF
          END IF
        ELSE
          TURBT(IP)=ZERO
          UTRB(IP)=ZERO
          VTRB(IP)=ZERO
          WTRB(IP)=ZERO
        END IF
      END IF
  100 CONTINUE
C================================================================PINT
C     OBTAIN DRAG AND HEAT TRANSFER COEFFICIENTS FROM USER SPECIFIED
C     SUBROUTINE USERPP
C----------------------------------------------------------------PINT
      CALL USERPP
C================================================================PINT
      DO 200 IP=1,NP
      IC=ICP(IP)
      UF=HALF*(UN(IC)+UN(IC-1))
      UG=UF
      VF=HALF*(VN(IC)+VN(IC-NXT))
      VG=VF
      IF(ISWIRL.EQ.1) THEN
        WF=WN(IC)
      ELSEIF(NZ.GT.1) THEN
        WF=HALF*(WN(IC)+WN(IC-NXYT))
      ELSE
        WF=ZERO
      END IF
      WG=WF
      IF(I3DP2D.EQ.1) THEN
        COSTH=XP(IP)/MAX(SQRT(XP(IP)**2+ZP(IP)**2),SMALL)
        SINTH=ZP(IP)/MAX(SQRT(XP(IP)**2+ZP(IP)**2),SMALL)
        UG=UF*COSTH-WF*SINTH
        WG=UF*SINTH+WF*COSTH
      END IF
      EPSAVE=EP(IP)
      UPSAVE=UP(IP)
      VPSAVE=VP(IP)
      WPSAVE=WP(IP)
C================================================================PINT
C ----- FOR THE SOLUTION OF PARTICLE TEMP. TWO DIFFERENT MODELS ARE
C       AVAILABLE. SELECT ICONDP=0 FOR THE SIMPLE UNIFORM TEMP
C                                1 FOR THE NEW HEATING MODEL
C                              BY Y.P. WAN, 12-10-97
      IF(ICONDP.EQ.1) THEN
        CALL PTCHEATING(HTCO(IP),DMPTC(IP),TEMP(IC),DT(IC1),
     1       TML(IP),HLTNT(MTYPE(IP)),STEBOL,EMSSP(IP),TEMAMB,
     1       TP(IP),MTYPE(IP),IBOIL(IP),IP,TIME)
      ELSE
c  --- note here the implicit solution of the energy equation
        ADTOM=AREAP(IP)*DT(IC1)/PMASS(IP)
        HASET3=HTCO(IP)+STEBOL*EMSSP(IP)*TP(IP)**3
C        EP(IP)=(EPSAVE+ADTOM*(HTCO(IP)*TEMP(IC)-STEBOL*EMSSP(IP)
C       1                      *TEMAMB**3)
C       1        +ADTOM*HASET3*(DTPDEP(IP)*EPSAVE-TP(IP)))
C       2       /(ONE+ADTOM*HASET3*DTPDEP(IP))
C ---- THE OLD EXPRESSION WAS NOT CORRECT AT THE TERM TEMAMB, IT SHOULD
C      BE:     BY Y.P. WAN, 12-10-97
        IF(IBOIL(IP).EQ.0) THEN
          EP(IP)=(EPSAVE+DMPTC(IP)*DT(IC1)*HLTNT(MTYPE(IP))/PMASS(IP)+
     1          ADTOM*(HTCO(IP)*TEMP(IC)+STEBOL*EMSSP(IP)*
     1                       TEMAMB**4)
     1        +ADTOM*HASET3*(DTPDEP(IP)*EPSAVE-TP(IP)))/
     2       (ONE+ADTOM*HASET3*DTPDEP(IP))
        END IF
      END IF
C================================================================PINT
      CDRVDT=CD(IP)*PIO2*RADP(IP)*RADP(IP)*RON(IC)*RELVEL(IP)*DT(IC1)
C=================== TURBULENT DISPERSION, WHEN DT > TSCALE =====PINT
C     HERE WE ADD CHANGE OF POSITION TO THE PARTICLE POSITION
C          EXPLICITLY.
C     WE ALSO CALCULATE THE PARTICLE VELOCITY CHANGE HERE. THOSE WILL
C          BE USED IN UPDATING PARTICLE VELOCITY IMPLICITLY BELOW.
C----------------------------------------------------------------PINT
      IF(ITURB.NE.0) THEN
       IF(DT(IC1).GT.TSCALE(IP)) THEN
        DRAGDT=CDRVDT/PMASS(IP)
        ATD=DRAGDT/DT(IC1)*TSCALE(IP)
        EXPMTD=EXP(-ATD)
        EXPMDT=EXP(-DRAGDT)
        PALP=(ONE-EXPMTD)/(ONE+EXPMTD)*(ONE-EXPMDT**2)
        RDRGDT=ONE/DRAGDT
        PALPTE=PALP*DT(IC1)*RDRGDT
        PBET=DT(IC1)*(TSCALE(IP)
     1             -TWO*TSCALE(IP)*(ONE-EXPMDT)*RDRGDT+PALPTE*RDRGDT)
        PGAM=TSCALE(IP)*(ONE-EXPMDT)-PALPTE
        PLAM=PGAM/PALP
C------------------------------------ FOR PARTICLE VELOCITY -----PINT
C     VELOCITY IS UPDATED IMPLICITLY BELOW USING UPRIMX, ETC.
C----------------------------------------------------------------PINT
        TURVEL=SQRT(FOURO3*PALP*TKE(IC))
C----- FOR X
        ETA1=ONE-TWO*FRAN(ZERO)
        ETA=ABS(ETA1)
        IND=INT(ONE+TWENTY*ETA)
        UPRIME=RERF(IND)+(ONE+TWENTY*ETA-DBLE(IND))
     1                  *(RERF(IND+1)-RERF(IND))
        UPRIMX=UPRIME*SIGN(ONE,ETA1)*TURVEL
C----- FOR Y
        ETA1=ONE-TWO*FRAN(ZERO)
        ETA=ABS(ETA1)
        IND=INT(ONE+TWENTY*ETA)
        VPRIME=RERF(IND)+(ONE+TWENTY*ETA-DBLE(IND))
     1                  *(RERF(IND+1)-RERF(IND))
        VPRIMY=VPRIME*SIGN(ONE,ETA1)*TURVEL
C----- FOR Z
        IF(I3DP2D.EQ.1 .OR. NZ.GT.1) THEN
        ETA1=ONE-TWO*FRAN(ZERO)
        ETA=ABS(ETA1)
        IND=INT(ONE+TWENTY*ETA)
        WPRIME=RERF(IND)+(ONE+TWENTY*ETA-DBLE(IND))
     1                  *(RERF(IND+1)-RERF(IND))
        WPRIMZ=WPRIME*SIGN(ONE,ETA1)*TURVEL
        END IF
C------------------------------------ FOR PARTICLE POSITION -----PINT
C     POSITION IS UPDATED EXPLICITLY
C----------------------------------------------------------------PINT
        TURTE=MAX((PBET-PLAM**2*PALP)*FOURO3*TKE(IC),ZERO)
        TURDIS=SQRT(TURTE)
C----- FOR X
        ETA1=ONE-TWO*FRAN(ZERO)
        ETA=ABS(ETA1)
        IND=INT(ONE+TWENTY*ETA)
        XPRIME=RERF(IND)+(ONE+20.*ETA-DBLE(IND))*(RERF(IND+1)-RERF(IND))
        XP(IP)=XP(IP)+XPRIME*SIGN(ONE,ETA1)*TURDIS+PLAM*UPRIME
C----- FOR Y
        ETA1=ONE-TWO*FRAN(ZERO)
        ETA=ABS(ETA1)
        IND=INT(ONE+TWENTY*ETA)
        YPRIME=RERF(IND)+(ONE+20.*ETA-DBLE(IND))*(RERF(IND+1)-RERF(IND))
        YP(IP)=YP(IP)+YPRIME*SIGN(ONE,ETA1)*TURDIS+PLAM*VPRIME
C----- FOR Z
        IF(I3DP2D.EQ.1 .OR. NZ.GT.1) THEN
        ETA1=ONE-TWO*FRAN(ZERO)
        ETA=ABS(ETA1)
        IND=INT(ONE+TWENTY*ETA)
        ZPRIME=RERF(IND)+(ONE+20.*ETA-DBLE(IND))*(RERF(IND+1)-RERF(IND))
        ZP(IP)=ZP(IP)+ZPRIME*SIGN(ONE,ETA1)*TURDIS+PLAM*ZPRIME
        END IF
       END IF
      ELSE
       UPRIMX=ZERO
       VPRIMY=ZERO
       WPRIMZ=ZERO
      END IF
C======================= UPDATE PARTICLE VELOCITY IMPLICITLY ====PINT
      UP(IP)=(PMASS(IP)*(UPSAVE+UPRIMX)+CDRVDT*(UG+UTRB(IP)))/
     1       (PMASS(IP)+CDRVDT)
      VP(IP)=(PMASS(IP)*(VPSAVE+VPRIMY)+CDRVDT*(VG+VTRB(IP)))/
     1       (PMASS(IP)+CDRVDT)
      WP(IP)=(PMASS(IP)*(WPSAVE+WPRIMZ)+CDRVDT*(WG+WTRB(IP)))/
     1       (PMASS(IP)+CDRVDT)
C================= MOMENTUM AND ENERGY SOURCE DUE TO PARCEL =====PINT
C     NOTE THAT QP IS NEGATIVE WHEN PARTICLES ARE HEATED.
C      AND THAT TRMP IS NEGATIVE WHEN PARTICLES ARE ACCELERATED.
C     QPDT = HEAT GAIN PER UNIT VOLUME DURING DT = HEATFLUX*DT/VOL
C            (ACTUAL ENERGY TRANSFER FROM FLUID TO PARTICLE ONLY.)
C     TRMP = MOMENTUM GAIN DURING DT
C     ABOVE RELATION IS USED TO AVOID ARITHMATIC OPERATIONS.
C----------------------------------------------------------------PINT
C     MOMENTUM SOURCE IS DISTRIBUTED AMONG ASSOCIATED MOMENTUM CELLS
C     EXCEPT FOR SWIRL FLOW.
C----------------------------------------------------------------PINT
C     SOURCE AND SINK DUE TO PARTICLE DRAG IN INTERNAL ENERGY,
C          TURBULENT KINETIC ENERGY, AND ITS DISSIPATION.
C     SEE P. J. O'ROURKE, KIVA II REPORT, LA-11560-MS (1989).
C     TURBULENT KINETIC ENERGY AND ITS DISSIPATION ARE AFFECTED
C          ONLY FOR THE CASE OF DT < TSCALE.
C----------------------------------------------------------------PINT
C     WSDT IS SUMMATION OF EFFECTS DUE TO INDIVIDUAL PARTICLES
C         WSDT IS SET ZERO IN OLDTIM AND USED IN PDVETC
C--------- SUMMATION OF HEAT AND MOMENTUM LOSS TO PARTICLES -----PINT
C ----  FOR NEW HEATING MODEL, QPDT IS CALC USING TEMP BY Y.P. WAN
      IF(ICONDP.EQ.1) THEN
        QPDT=PARTN(IP)*DT(IC1)*AREAP(IP)*HTCO(IP)
     1     *(TP(IP)-TEMP(IC))/VOL(IC)
      ELSE
        QPDT=PARTN(IP)*DT(IC1)*AREAP(IP)*HTCO(IP)
     1     *(TP(IP)+DTPDEP(IP)*(EP(IP)-EPSAVE)-TEMP(IC))/VOL(IC)
      END IF
      TRMPX=PARTN(IP)*CDRVDT*(UP(IP)-UG-UTRB(IP))
      TRMPY=PARTN(IP)*CDRVDT*(VP(IP)-VG-VTRB(IP))
      TRMPZ=PARTN(IP)*CDRVDT*(WP(IP)-WG-WTRB(IP))
C-------------------------------------- FOR CYLINDRICAL 2-D -----PINT
      IF(ICYL.EQ.1) THEN
C----------------------------------------------------------------PINT
C     AT SYMMETRY AXIS, WE DO NOT NEED SPECIAL FIX, SINCE MOMENTUM
C     CELL FOR RORU OCCUPIES PART OF THE CELL NEXT TO THE AXIS.
C----------------------------------------------------------------PINT
C----------FOR PLANE CYLIDRICAL CASE
       IF(I3DP2D.EQ.0) THEN
        RORU(IC)=RORU(IC)
     1          +TRMPX*(XC(IC)+HALF*DX(IC))/(VOL(IC)+VOL(IC+1))
        RORU(IC-1)=RORU(IC-1)
     1            +TRMPX*(XC(IC-1)+HALF*DX(IC-1))/(VOL(IC)+VOL(IC-1))
C----------FOR PSEUDO 3-D CASE
       ELSE
        RORU(IC)=RORU(IC)+(TRMPX*COSTH+TRMPZ*SINTH)
     1                    *(XC(IC)+HALF*DX(IC))/(VOL(IC)+VOL(IC+1))
        RORU(IC-1)=RORU(IC-1)+(TRMPX*COSTH+TRMPZ*SINTH)
     1                  *(XC(IC-1)+HALF*DX(IC-1))/(VOL(IC)+VOL(IC-1))
        RORW(IC)=RORW(IC)+(TRMPZ*COSTH-TRMPX*SINTH)*XC(IC)/VOL(IC)
       END IF
C----------
       RORV(IC)=RORV(IC)+TRMPY*XC(IC)/(VOL(IC)+VOL(IC+NXT))
       RORV(IC-NXT)=RORV(IC-NXT)+TRMPY*XC(IC)/(VOL(IC)+VOL(IC-NXT))
C----------
       ROER(IC)=ROER(IC)+QPDT*XC(IC)+XC(IC)/VOL(IC)
     1                               *(TRMPX*(UP(IP)-UG-UTRB(IP))
     2                                +TRMPY*(VP(IP)-VG-VTRB(IP))
     3                                +TRMPZ*(WP(IP)-WG-WTRB(IP)))
       IF(ITURB.NE.0) THEN
        WSDT(IC)=WSDT(IC)
     1 +(TRMPX*UTRB(IP)+TRMPY*VTRB(IP)+TRMPZ*WTRB(IP))*XC(IC)/VOL(IC)
       END IF
C-------------------------------- FOR CARTESIAN 2-D AND 3-D -----PINT
      ELSE
       RORU(IC)=RORU(IC)+TRMPX/(VOL(IC)+VOL(IC+1))
       RORU(IC-1)=RORU(IC-1)+TRMPX/(VOL(IC)+VOL(IC-1))
       RORV(IC)=RORV(IC)+TRMPY/(VOL(IC)+VOL(IC+NXT))
       RORV(IC-NXT)=RORV(IC-NXT)+TRMPY/(VOL(IC)+VOL(IC-NXT))
       IF(NZ.GT.1) THEN
        RORW(IC)=RORW(IC)+TRMPZ/(VOL(IC)+VOL(IC+NXYT))
        RORW(IC-NXYT)=RORW(IC-NXYT)+TRMPZ/(VOL(IC)+VOL(IC-NXYT))
       END IF
C----------
       ROER(IC)=ROER(IC)+QPDT+(TRMPX*(UP(IP)-UG-UTRB(IP))
     1                        +TRMPY*(VP(IP)-VG-VTRB(IP))
     2                        +TRMPZ*(WP(IP)-WG-WTRB(IP)))/VOL(IC)
C----------
       WSDT(IC)=WSDT(IC)
     1       +(TRMPX*UTRB(IP)+TRMPY*VTRB(IP)+TRMPZ*WTRB(IP))/VOL(IC)
C----------
      END IF
  200 CONTINUE
C================================================================PINT
      RETURN
      END
