*DECK INJECT
      SUBROUTINE INJECT
C==============================================================INJECT
C     THIS ROUTINE GENERATES PARTICLES
C----------
C     CALLED BY LAVA
C     CALLS FRAN
C==============================================================INJECT
      INCLUDE 'COML.h'
      INCLUDE 'COMMPTC.h'
C--------------------------------------------------------------INJECT
      DOUBLE PRECISION FRAN 
      EXTERNAL FRAN
C--------------------------------------------------------------INJECT
      INTEGER N,NPPTS,NOZ,IRAND
      DOUBLE PRECISION FRANS1,FRANS2,FRAND,RANBAK,FRANV
      DOUBLE PRECISION SPRANG,PHI,VAXI,VNORM,VOTHR,VINJ,PSZDEV
      DOUBLE PRECISION PARVOL,SUMS,SUML
C==============================================================INJECT
      DO 40 NOZ=1,NUMNOZ
      TM2INJ(NOZ)=TSPMAS(NOZ)*(TIME-T1INJ(NOZ))
      NPPTS=INT((TM2INJ(NOZ)-TM1INJ(NOZ))/PMINJ(NOZ))
      IF(NPPTS.LT.1) GO TO 40
      IF(NP+NPPTS.GT.NPAR) GO TO 50
      DO 30 N=1,NPPTS
      NP=NP+1
c ----  FRAN is a Random number generator
      FRANS1=FRAN(ZERO)
      FRANS2=FRAN(ZERO)
      FRAND =FRAN(ZERO)
      FRANV =FRAN(ZERO)
      RANBAK=FRAN(ZERO)*DT(IC1)
c ----  instead of the random of the injection parameter, fix its value
c      by Y.P.WAN
c      frans1=0.5
c      frans2=0.5
c      frand=0.5
c      franv=0.5
c      ranbak=0.5*dt(ic1)
C------------------------------------ SOLID CONE SPRAY -----INJECT
C     TILTXZ IS THE MEAN SPRAY ANGLE IN X-Z SPACE, TILTXY (EMBODIED
C     IN THE UNIT VECTORS) IS THE MEAN SPRAY ANGLE IN X-Y SPACE,
C     AND THE ANGLE PHI ENSURES A CIRCULAR CROSS-SECTION THROUGH
C     THE SPRAY JET, WHOSE THICKNESS IS GIVEN BY THE ANGLE CONE:
C--------------------------------------------------------------INJECT
      VINJ=VELINJ(NOZ)*(ONE+VDISP(NOZ)*(FRANV-HALF))*PAMB/P(ICNOZ(NOZ))
C---- FOR PLANAR 2-D AND PLANE AXISYMMETRIC 2-D W/O SWIRL -----INJECT
      IF(NZ.EQ.1 .AND. ISWIRL.EQ.0 .AND. I3DP2D.EQ.0) THEN
c  --    The SPRANG shoud be [-0.5CONE, 0.5CONE], by Y.P. WAN 10/29/97
c        SPRANG=TILTX(NOZ)+CONE(NOZ)*ABS(FRANS1-HALF)
        SPRANG=TILTX(NOZ)+CONE(NOZ)*(FRANS1-HALF)
        UP(NP)=VINJ*COS(SPRANG)
        VP(NP)=VINJ*SIN(SPRANG)
        WP(NP)=ZERO
C---------- FOR PSEUDO 3-D (W AND W/O SWIRL) AND REAL 3-D -----INJECT
      ELSE
       SPRANG=CONE(NOZ)*ABS(FRANS1-HALF)
       PHI=PI2*FRANS2
       VAXI =VINJ*COS(SPRANG)
       VNORM=VINJ*SIN(SPRANG)*COS(PHI)
       VOTHR=VINJ*SIN(SPRANG)*SIN(PHI)
       UP(NP)=VAXI*EAVEC(NOZ,1)+VNORM*ENVEC(NOZ,1)+VOTHR*EOVEC(NOZ,1)
       VP(NP)=VAXI*EAVEC(NOZ,2)+VNORM*ENVEC(NOZ,2)+VOTHR*EOVEC(NOZ,2)
       WP(NP)=VAXI*EAVEC(NOZ,3)+VNORM*ENVEC(NOZ,3)+VOTHR*EOVEC(NOZ,3)
      END IF
C--------------------------------------------------------------INJECT
C     TURBULENT FLUCTUATION IS SET TO BE ZERO.
C     IF TURBULENT, THESE QUANTITIES ARE CALCULATED AT SUB. PINT.
C--------------------------------------------------------------INJECT
      UTRB(NP)=ZERO
      VTRB(NP)=ZERO
      WTRB(NP)=ZERO
      TURBT(NP)=ZERO
C--------------------------------------------------------------INJECT
C     INITIALIZE XP,YP,ZP AT A RANDOM POINT ALONG THE PARCEL
C     TRAJECTORY --> BEHIND <-- THE INJECTOR, AS PARCEL WILL
C     BE IMMEDIATELY MOVED IN SUBR. PMOVE, WHERE REAL LOCATION
C     WILL BE FOUND.
C--------------------------------------------------------------INJECT
C     XP, YP, AND ZP ARE IN RECTANGULAR COORDINATES IN ALL CASES
C--------------------------------------------------------------INJECT
      XP(NP)=XINJ(NOZ)-RANBAK*UP(NP)
      YP(NP)=YINJ(NOZ)-RANBAK*VP(NP)
      ZP(NP)=ZINJ(NOZ)-RANBAK*WP(NP)
      TP(NP)=TPINJ(NOZ)
C--------------------------------------------------------------INJECT
C     PARTICLE SIZE IS DETERMINED FROM GIVEN DISTRIBUTION
C--------------------------------------------------------------INJECT
      IRAND=2
      DO 100 I=2,NPSDIS(NOZ)+1
      IF(FRAND.GT.PSZINT(I-1,NOZ).AND.FRAND.LE.PSZINT(I,NOZ)) IRAND=I
  100 CONTINUE
      PSZDEV=(FRAND-PSZINT(IRAND-1,NOZ))/
     1       (PSZINT(IRAND,NOZ)-PSZINT(IRAND-1,NOZ))
      RADP(NP)=RPMIN(NOZ)+(DBLE(IRAND-2)+PSZDEV)*PSZDEL(NOZ)
C--------------------------------------------------------------INJECT
C     PARTICLES ARE ASSUMED TO BE SPHERICAL SHAPE.
C     THIS PART SHOULD BE MODIFIED IF PARTICLES ARE NOT SPHERICAL.
C--------------------------------------------------------------INJECT
      AREAP(NP)=PI4*RADP(NP)*RADP(NP)
      PARVOL=PI4O3*RADP(NP)**3
C--------------------------------------------------------------INJECT
C     THERMODYNAMIC VARIABLES OF PARTICLES
C     THIS LOGIC IS NOT FULLY GENERALIZED FOR MULTICOMPONENT
C     PARTICLES YET, EVEN THOUGH IT APPEARS TO HAVE THAT CAPABILITY
C     GENERALIZED PSTATE NEEDS TO BE DEVELOPED
C--------------------------------------------------------------INJECT
C     PARTICLE THERMOPHYSICAL PROPERTIES FOR INJECTED PARTICLES:
C     RSPHS, RSPHL, TMM, TML, EPM, AND EPL TO BE USED AT SUBROUTINE
C     PSTATE
C     EP IS SET TO BE ZERO AT 0 K.
C     EPM = EP AT TMM (ERG/GRAM)
C     EPL = EP AT TML (ERG/GRAM)
C--------------------------------------------------------------INJECT
      PMASS(NP)=ZERO
      EMSSP(NP)=ZERO
      TMM(NP)=LARGE
      TML(NP)=SMALL
      SUMS=ZERO
      SUML=ZERO
      DO 200 IPSP=1,NPSP
        PMSP(NP,IPSP)=PARVOL*ROPI(IPSP)*PXINJ(NOZ,IPSP)
        PMASS(NP)=PMASS(NP)+PMSP(NP,IPSP)
        EMSSP(NP)=EMSSPI(IPSP)*PXINJ(NOZ,IPSP)
C  ****** Warning: error in the lava TMM(np) will be the minum of
c         all species.  found by Y.P. WAN 12-30-97
c        TMM(NP)=MIN(TMM(NP),TMMI(IPSP))
c        TML(NP)=MAX(TML(NP),TMLI(IPSP))
        SUMS=SUMS+PSPHS(IPSP)*PYINJ(NOZ,IPSP)
        SUML=SUML+PSPHL(IPSP)*PYINJ(NOZ,IPSP)
  200 CONTINUE
C   INTRODUCING PARTICLE TYPE FOR MORE MATERIALS. HERE JUST FOR TWO MATERIALS
C   AND THE NEW DETERMINATION OF TMM AND TML OF PARTICLE NP
C                                     BY Y.P. WAN  12-30-97
      MTYPE(NP)=1
      IF(PMSP(NP,1)/PMASS(NP).LE.HALF) MTYPE(NP)=2
      TMM(NP)=TMMI(MTYPE(NP))
      TML(NP)=TMLI(MTYPE(NP))
      RPSPHS(NP)=ONE/SUMS
      RPSPHL(NP)=ONE/SUML
      EPM(NP)=SUMS*TMM(NP)
      EPL(NP)=ZERO
      DO 300 IPSP=1,NPSP
        EPL(NP)=EPL(NP)+(PSPHS(IPSP)*TMMI(IPSP)+HLTNT(IPSP)
     1           +PSPHL(IPSP)*(TML(NP)-TMLI(IPSP)))*PYINJ(NOZ,IPSP)
  300 CONTINUE
C--------------------------------------------------------------INJECT
      IF(TP(NP).LT.TMM(NP)) THEN
        DTPDEP(NP)=RPSPHS(NP)
        EP(NP)=SUMS*TP(NP)
      ELSEIF(TP(NP).GE.TMM(NP) .AND. TP(NP).LE.TML(NP)) THEN
        DTPDEP(NP)=(TML(NP)-TMM(NP))/(EPL(NP)-EPM(NP))
        EP(NP)=EPM(NP)+(TP(NP)-TMM(NP))/DTPDEP(NP)
      ELSE
        DTPDEP(NP)=RPSPHL(NP)
        EP(NP)=EPL(NP)+SUML*(TP(NP)-TML(NP))
      END IF
C--------------------------------------------------------------INJECT
      PARTN(NP)=PMINJ(NOZ)/PMASS(NP)
      ICP(NP)=ICNOZ(NOZ)
C  -----  INITIALIZING THE PARTICLE TEMP FIELD IF NEW HEATING MODEL USED
C                             BY Y.P. WAN, 12-10-97
      IF(ICONDP.EQ.1) THEN
        SOLVS(NP)=1
        SOLVL(NP)=0
        SOLVRS(NP)=0
        RD(NP)=RADP(NP)*1.D-2
        DRD(NP)=0.D0
        RM(NP)=RD(NP)
        RRS(NP)=RD(NP)
        VI(NP)=0.D0
        VIRS(NP)=0.D0
        DO 310 I=1,MGD
        TPCD(NP,I)=TP(NP)
 310    CONTINUE
        DO 320 I=2,NGD
        TPCD(NP,I+MGD-1)=TP(NP)
 320    CONTINUE
        DO 330 I=2,LGD
        TPCD(NP,I+MGD+NGD-2)=TP(NP)
 330    CONTINUE
      END IF
C 
   30 CONTINUE
      TM1INJ(NOZ)=TM1INJ(NOZ)+DBLE(NPPTS)*PMINJ(NOZ)
   40 CONTINUE
C==============================================================INJECT
      RETURN
   50 WRITE(NFOUT,900) TIME,NCYC
      WRITE(*,900)     TIME,NCYC
      CALL EXITA(11)
C
  900 FORMAT('PARTICLE STORAGE EXCEEDED AT TIME=',1PE12.5,' CYCLE',I7)
C==============================================================INJECT
      END
