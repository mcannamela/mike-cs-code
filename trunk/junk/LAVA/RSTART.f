*DECK RSTART
      SUBROUTINE RSTART
C==============================================================RSTART
C
C
C
C----------
C     CALLED BY LAVA
C==============================================================RSTART
      INCLUDE 'COML.h'
      INCLUDE 'COMMPTC.h'
C==============================================================RSTART
      OPEN (UNIT=NFRST,FILE='lavadump',STATUS='OLD',
     1      ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND NFRST
C------------------------------------- CONTROL PARAMETERS -----RSTART
      READ(NFRST) ADC,BDC,TIME,TIMOVI,DTOLD
      READ(NFRST) INCOMP,ISWIRL,ICYL,NODIFF,NLTE,NCYC,NMJ
C-------------------------------------- MESH INFORMATIONS -----RSTART
      READ(NFRST) XL,YL,ZL,DX,DY,DZ,RDX,RDY,RDZ,XC,X,Y,Z,AREA,VOL,
     1     HRDXCR,HRDYCF,HRDZCT
      READ(NFRST) I1,I2,J1,J2,K1,K2,IC1U,IC1V,IC1W,IC1,IC2,
     1     NX,NY,NZ,NXYT,NXYZT,IGRIDX,IGRIDY,IGRIDZ,NIT1,
     2     NXS,NXL,NYS,NYL,NZS,NZL
C---------------------------------------------- R FACTORS -----RSTART
      READ(NFRST) RC,RR,RF,RT,RRC
C----------------------------------------- MAIN VARIABLES -----RSTART
      READ(NFRST) U,V,W,UN,VN,WN,RORU,RORV,RORW,
     1     ROE,ROEN,ROER,TEMP,ROEE,ROEEN,ROEER,TE,
     2     RO,RON,SPD,SPDR,SPDN
C------------------------------------ BOUNDARY CONDITIONS -----RSTART
      READ(NFRST) PBCL,PBCR,PBCD,PBCF,PBCB,PBCT,
     1     SDNL,SDNR,SDND,SDNF,SDNB,SDNT,
     2     EDENL,EDENR,EDEND,EDENF,EDENB,EDENT,
     3     TEMPL,TEMPR,TEMPD,TEMPF,TEMPB,TEMPT,
     4     EDEEL,EDEER,EDEED,EDEEF,EDEEB,EDEET,
     5     TEEL,TEER,TEED,TEEF,TEEB,TEET,
     6     PRESL,PRESR,PRESD,PRESF,PRESB,PREST,
     7     UVELL,UVELR,UVELD,UVELF,UVELB,UVELT,
     8     VVELL,VVELR,VVELD,VVELF,VVELB,VVELT,
     9     WVELL,WVELR,WVELD,WVELF,WVELB,WVELT,
     9     TKEL,TKER,TKED,TKEF,TKEB,TKET,
     1     EPSL,EPSR,EPSD,EPSF,EPSB,EPST
      READ(NFRST) IPBC
C----------------- PRESSURE, ACCELERATION, AND DIVERGENCE -----RSTART
      READ(NFRST) P,PN,Q,QN,DIV,DIVOLD,DIVE,UDELR
C---------------------------------------- STATE VARIABLES -----RSTART
      READ(NFRST) MW,RMW,GAMMA,EK,HK,QOM
      READ(NFRST) IELC
C------------------------ EQUILIBRIUM CHEMISTRY VARIABLES -----RSTART
      READ(NFRST) OMGDOT
C------------------------- VARIABLES IN TURBULENCE MODELS -----RSTART
      READ(NFRST) ROTKER,TKE,TKEN,ROEPSR,EPS,EPSN,SGSL,WSDT
C------------------------------------------------ PARTICLES ---RSTART
C THIS PART HAS BEEN CHANGED BY Y.P. WAN 1-7-98 IN ORDER TO 
C REDUCE THE SIZE OF DUMP FILE
      READ(NFRST) NP,ICNOZ
      READ(NFRST) RANB,RANS,TM1INJ,TM2INJ
      READ(NFRST) (ICP(I),I=1,NP)
      READ(NFRST) (RADP(I),I=1,NP)
      READ(NFRST) (PARTN(I),I=1,NP)
      READ(NFRST) (PMASS(I),I=1,NP)
      READ(NFRST) (AREAP(I),I=1,NP)
      READ(NFRST) (XP(I),I=1,NP)
      READ(NFRST) (YP(I),I=1,NP)
      READ(NFRST) (ZP(I),I=1,NP)
      READ(NFRST) (UP(I),I=1,NP)
      READ(NFRST) (VP(I),I=1,NP)
      READ(NFRST) (WP(I),I=1,NP)
      READ(NFRST) (UTRB(I),I=1,NP)
      READ(NFRST) (VTRB(I),I=1,NP)
      READ(NFRST) (WTRB(I),I=1,NP)
      READ(NFRST) (TURBT(I),I=1,NP)
      READ(NFRST) (TP(I),I=1,NP)
      READ(NFRST) (EP(I),I=1,NP)
      READ(NFRST) (DTPDEP(I),I=1,NP)
      READ(NFRST) (TMM(I),I=1,NP)
      READ(NFRST) (TML(I),I=1,NP)
      READ(NFRST) (RPSPHS(I),I=1,NP)
      READ(NFRST) (RPSPHL(I),I=1,NP)
      READ(NFRST) (EPM(I),I=1,NP)
      READ(NFRST) (EPL(I),I=1,NP)
      READ(NFRST) (EMSSP(I),I=1,NP)
      READ(NFRST) (MTYPE(I),I=1,NP)
      DO 10 J=1,NPSP
      READ(NFRST) (PMSP(I,J),I=1,NP)
 10   CONTINUE
C------------------------------------- TRANSPORT PROPERTIES ---RSTART
      READ(NFRST) VISC,COND,CONDE,RADLSS,VIST,DCOEF
C-------------------------------- VARIABLES FOR PARTICLE HEATING
C----------------BY Y.P. WAN 1-6-98  --------------------------RSTART
      READ(NFRST) NMGD,MGD,NGD,LGD,NCS
      READ(NFRST) (Y1C(I),I=1,NHALF)
      READ(NFRST) (Y2C(I),I=1,NHALF)
      READ(NFRST) (Y3C(I),I=1,NHALF)
      READ(NFRST) (Y1B(I),I=1,NHALF)
      READ(NFRST) (Y2B(I),I=1,NHALF)
      READ(NFRST) (Y3B(I),I=1,NHALF)
      READ(NFRST) (DLY1(I),I=1,NHALF)
      READ(NFRST) (DLY2(I),I=1,NHALF)
      READ(NFRST) (DLY3(I),I=1,NHALF)
      READ(NFRST) (RD(I),I=1,NP)
      READ(NFRST) (RM(I),I=1,NP)
      READ(NFRST) (RRS(I),I=1,NP)
      READ(NFRST) (DRD(I),I=1,NP)
      READ(NFRST) (VI(I),I=1,NP)
      READ(NFRST) (VIRS(I),I=1,NP)
      READ(NFRST) (SOLVS(I),I=1,NP)
      READ(NFRST) (SOLVL(I),I=1,NP)
      READ(NFRST) (SOLVRS(I),I=1,NP)
      DO 20 J=1,NMGD
      READ(NFRST) (TPCD(I,J),I=1,NP)
 20   CONTINUE
C--------------------------------------------------------------RSTART
      CLOSE(UNIT=NFRST)
      WRITE(NFOUT,100) NCYC,TIME
  100 FORMAT(5X,'-----------------------------------------------',/,
     1       5X,'RESTART FROM DUMP AT NCYC = ',I8,', TIME =',E12.5)
C==============================================================RSTART
      RETURN
      END
