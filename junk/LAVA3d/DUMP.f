*DECK DUMP 
      SUBROUTINE DUMP 
C================================================================DUMP
C
C
C
C----------
C     CALLED BY LAVA
C================================================================DUMP
      INCLUDE 'COML.h'
      INCLUDE 'COMMPTC.h'
c      INTEGER I,J
C================================================================DUMP
      OPEN (UNIT=NFDMP,FILE='lavadump',STATUS='UNKNOWN',
     1      ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND NFDMP
C--------------------------------------- CONTROL PARAMETERS -----DUMP
      WRITE(NFDMP) ADC,BDC,TIME,TIMOVI,DTOLD
      WRITE(NFDMP) INCOMP,ISWIRL,ICYL,NODIFF,NLTE,NCYC,NMJ
C---------------------------------------- MESH INFORMATIONS -----DUMP
      WRITE(NFDMP) XL,YL,ZL,DX,DY,DZ,RDX,RDY,RDZ,XC,X,Y,Z,AREA,VOL,
     1     HRDXCR,HRDYCF,HRDZCT
      WRITE(NFDMP) I1,I2,J1,J2,K1,K2,IC1U,IC1V,IC1W,IC1,IC2,
     1     NX,NY,NZ,NXYT,NXYZT,IGRIDX,IGRIDY,IGRIDZ,NIT1,
     2     NXS,NXL,NYS,NYL,NZS,NZL
C------------------------------------------------ R FACTORS -----DUMP
      WRITE(NFDMP) RC,RR,RF,RT,RRC
C------------------------------------------- MAIN VARIABLES -----DUMP
      WRITE(NFDMP) U,V,W,UN,VN,WN,RORU,RORV,RORW,
     1     ROE,ROEN,ROER,TEMP,ROEE,ROEEN,ROEER,TE,
     2     RO,RON,SPD,SPDR,SPDN
C-------------------------------------- BOUNDARY CONDITIONS -----DUMP
      WRITE(NFDMP) PBCL,PBCR,PBCD,PBCF,PBCB,PBCT,
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
      WRITE(NFDMP) IPBC
C------------------- PRESSURE, ACCELERATION, AND DIVERGENCE -----DUMP
      WRITE(NFDMP) P,PN,Q,QN,DIV,DIVOLD,DIVE,UDELR
C------------------------------------------ STATE VARIABLES -----DUMP
      WRITE(NFDMP) MW,RMW,GAMMA,EK,HK,QOM
      WRITE(NFDMP) IELC
C-------------------------- EQUILIBRIUM CHEMISTRY VARIABLES -----DUMP
      WRITE(NFDMP) OMGDOT
C--------------------------- VARIABLES IN TURBULENCE MODELS -----DUMP
      WRITE(NFDMP) ROTKER,TKE,TKEN,ROEPSR,EPS,EPSN,SGSL,WSDT
C------------------------------------------------ PARTICLES -----DUMP
C THIS PART HAS BEEN CHANGED BY Y.P. WAN 1-7-98 IN ORDER TO 
C REDUCE THE SIZE OF DUMP FILE
      WRITE(NFDMP) NP,ICNOZ
      WRITE(NFDMP) RANB,RANS,TM1INJ,TM2INJ
      WRITE(NFDMP) (ICP(I),I=1,NP)
      WRITE(NFDMP) (RADP(I),I=1,NP)
      WRITE(NFDMP) (PARTN(I),I=1,NP)
      WRITE(NFDMP) (PMASS(I),I=1,NP)
      WRITE(NFDMP) (AREAP(I),I=1,NP)
      WRITE(NFDMP) (XP(I),I=1,NP)
      WRITE(NFDMP) (YP(I),I=1,NP)
      WRITE(NFDMP) (ZP(I),I=1,NP)
      WRITE(NFDMP) (UP(I),I=1,NP)
      WRITE(NFDMP) (VP(I),I=1,NP)
      WRITE(NFDMP) (WP(I),I=1,NP)
      WRITE(NFDMP) (UTRB(I),I=1,NP)
      WRITE(NFDMP) (VTRB(I),I=1,NP)
      WRITE(NFDMP) (WTRB(I),I=1,NP)
      WRITE(NFDMP) (TURBT(I),I=1,NP)
      WRITE(NFDMP) (TP(I),I=1,NP)
      WRITE(NFDMP) (EP(I),I=1,NP)
      WRITE(NFDMP) (DTPDEP(I),I=1,NP)
      WRITE(NFDMP) (TMM(I),I=1,NP)
      WRITE(NFDMP) (TML(I),I=1,NP)
      WRITE(NFDMP) (RPSPHS(I),I=1,NP)
      WRITE(NFDMP) (RPSPHL(I),I=1,NP)
      WRITE(NFDMP) (EPM(I),I=1,NP)
      WRITE(NFDMP) (EPL(I),I=1,NP)
      WRITE(NFDMP) (EMSSP(I),I=1,NP)
      WRITE(NFDMP) (MTYPE(I),I=1,NP)
      DO 10 J=1,NPSP
      WRITE(NFDMP) (PMSP(I,J),I=1,NP)
 10   CONTINUE
C------------------------------------- TRANSPORT PROPERTIES -----DUMP
      WRITE(NFDMP) VISC,COND,CONDE,RADLSS,VIST,DCOEF
C-------------------------------- VARIABLES FOR PARTICLE HEATING
C----------------BY Y.P. WAN 1-6-98  ----------------------------DUMP
      WRITE(NFDMP) NMGD,MGD,NGD,LGD,NCS 
      WRITE(NFDMP) (Y1C(I),I=1,NHALF)
      WRITE(NFDMP) (Y2C(I),I=1,NHALF)
      WRITE(NFDMP) (Y3C(I),I=1,NHALF)
      WRITE(NFDMP) (Y1B(I),I=1,NHALF)
      WRITE(NFDMP) (Y2B(I),I=1,NHALF)
      WRITE(NFDMP) (Y3B(I),I=1,NHALF)
      WRITE(NFDMP) (DLY1(I),I=1,NHALF)
      WRITE(NFDMP) (DLY2(I),I=1,NHALF)
      WRITE(NFDMP) (DLY3(I),I=1,NHALF)
      WRITE(NFDMP) (RD(I),I=1,NP)
      WRITE(NFDMP) (RM(I),I=1,NP)
      WRITE(NFDMP) (RRS(I),I=1,NP)
      WRITE(NFDMP) (DRD(I),I=1,NP)
      WRITE(NFDMP) (VI(I),I=1,NP)
      WRITE(NFDMP) (VIRS(I),I=1,NP)
      WRITE(NFDMP) (SOLVS(I),I=1,NP)
      WRITE(NFDMP) (SOLVL(I),I=1,NP)
      WRITE(NFDMP) (SOLVRS(I),I=1,NP)
      DO 20 J=1,NMGD
      WRITE(NFDMP) (TPCD(I,J),I=1,NP)
 20   CONTINUE
C----------------------------------------------------------------DUMP
      CLOSE(UNIT=NFDMP)
      WRITE(NFOUT,100) NCYC,TIME
  100 FORMAT(5X,'-----------------------------------------------',/,
     1       5X,'TAPE DUMPED AT NCYC = ',I8,', TIME =',E12.5)
C================================================================DUMP
      RETURN
      END
