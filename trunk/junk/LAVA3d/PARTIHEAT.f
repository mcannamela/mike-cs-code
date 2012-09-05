      SUBROUTINE PTCHEATING(HTCO1,DMPTC1,TG,DT,TMELT,LATM1,STBOL1,
     1           EMISS,TAMB,TPS,MTYPE,IBOIL,IP,TIME)
C =======================================================================
C    THIS SUBROUTINE IS USED TO CALCULATE THE TEMPERATURE DISTRIBUTION OF
C    A SINGLE PARTICLE. HEAT CONDUCTION AND MELTING INTERFACE ARE
C    CONSIDERED. NOTE, ALL THE UNIT USED HERE ARE SI UNIT. THEREFORE ALL
C    THE UNITS INPUTED FROM LAVA SHOULD BE TRANSFERRED TO SI
C    12-10-97   BY Y.P. WAN
C ======
C    Since the melting interface is mostly controlled by the heat transfer
C    and the moving velocity is not so significantly dependent on the
C    melting kinetics (which considers the overheating), therefore, the
C    temperature at the interface can be simply assumed to be melting temp.
C    The iteration based on the melting kinetics can then be spared.
C    for this purpose, we introduce a index which controls how to calculate
C    the melting interface moving. IMKINE=0, NO KINETIC, =1, KINETIC
C =======================================================================
C
C ----- INPUT PARAMETERS INCLUDING:
C    HTCO:   HEAT TRANSFER COEFFICIENT
C    TG  :   TEMPERATURE OF GASPHASE AT THE LOCAL CELL
C    DT  :   TIME STEP USED IN CFD CODE
C    IP  :   THE NUMBER OF THE CALCULATED PARTICLE
C
      IMPLICIT NONE
      INCLUDE "COMMPTC.h"
      INTEGER IP,I,J,ITVI,MTYPE,IBOIL,IMKINE
      DOUBLE PRECISION DT,TPS,TG,HTCO,QCON,QRAD,DMPTC,QEVAP
      DOUBLE PRECISION HTCO1,LATM1,STBOL1,TAMB,DMPTC1
      DOUBLE PRECISION TMIDS(NHALF),TMIDL(NHALF),AP(NHALF),AE(NHALF),
     &                 AW(NHALF),BB(NHALF),TMIDRS(NHALF),
     &                 SP(NHALF),SC(NHALF),TMELT,DRM,RMTMP,DRS,RSTMP,
     &                 YYC,RJ,RJ1
      DOUBLE PRECISION FLUXL,FLUXS,VIMID,EPS,DVI,RELAX,LATM,MUI,TMI,
     &                 TGRS,DTT,DT0,DRMMIN,STBOL,EMISS,PI,TIME
      DOUBLE PRECISION CONDSOL,CONDLIQ,DENSSOL,DENSLIQ,CPLIQ,CPSOL,
     &                 MUINTFACE,LATVAP
      CHARACTER*1 TAB
      EXTERNAL CONDSOL,CONDLIQ,DENSSOL,DENSLIQ,CPLIQ,CPSOL,MUINTFACE,
     &         LATVAP
      DATA EPS,RELAX/1.D-3,0.3D0/
      DATA PI/3.1415926/
C  IMKINE=0, NO KINETIC, =1 WITH KINETIC
      DATA IMKINE/0/
C ----- SETTING UP THE BOUNDARY CONDITIONS, MATERIALS PROPERTIES AND
C       CONTROL PARAMETERS
C     NCS=0 FOR PLANE COORDINATE
C         1 FOR CYLINDRICAL CO.
C         2 FOR SPHERICAL CO.
C     SOLVS(IP): 1 SOLVE THE SOLID DOMAIN, 0 NOT SOLVE
C     SOLVL(IP): 1 SOLVE THE LIQUID DOMAIN, 0 NOT SOLVE
C     HTCO  : HEAT TRANSFER COEFFICIENT, OBTAINED FROM LAVA
C     TMIDS  : INTERMIDIATE TEMP KEEP THE 1-D TEMP OF SOLID
C     TMIDL  : INTERMIDIATE TEMP KEEP THE 1-D TEMP OF LIQUID
c
      MUI=MUINTFACE(MTYPE)
      DRMMIN=DLY1(MGD-1)*RD(IP)
C   ---- TRANSFER THE UNIT TO SI UNIT HTCO (erg/cm^2sK -> J/m^2sK)
C                                    DMPTC (g/s -> kg/s)
C                                     LATM (erg/g -> J/kg)
C                                    STBOL (erg/cm^2K^4 -> J/m^2K^4)
      HTCO=HTCO1*1.D-3
      DMPTC=DMPTC1*1.D-3
      LATM=LATM1*1.D-4
      STBOL=STBOL1*1.D-3
C  --- ITERATION WITHIN DT FOR PARTICLE IP
      DT0=0.D0
      DTT=DT
  3   CONTINUE
C =========== (1,0,0)
      IF(SOLVS(IP).EQ.1.AND.SOLVL(IP).EQ.0.AND.
     &   SOLVRS(IP).EQ.0) THEN
C  --   NO MELTING TAKES PLACE. AND NO VAPORIZATION WILL BE CONSIDERED
        DO 21 I=1,MGD
        TMIDS(I)=TPCD(IP,I)
        SP(I)=0.D0
        SC(I)=0.D0
  21    CONTINUE
        TPS=TPCD(IP,MGD)
        QCON=HTCO*(TG-TPS)
        QRAD=STBOL*EMISS*(TAMB**4-TPS**4)
        QEVAP=DMPTC*LATVAP(TPS,MTYPE)/(4.0D0*PI*RD(IP)*RD(IP))
c           write(*,*)'---- tg,tpcd',tg,TPCD(IP,mgd),
c     &                TPCD(IP,1),qcon,qrad,qevap
        CALL COEFF(AP,AE,AW,BB,MGD,TMIDS,Y1C,Y1B,SP,SC,
     &           DTT,RD(IP),0.D0,DRD(IP),0.D0,NCS,1,MTYPE)
        AP(1)=1.D0
        AE(1)=1.D0
        AW(1)=0.D0
        BB(1)=0.D0
        AP(MGD)=1.D0
        AE(MGD)=0.D0
        AW(MGD)=1.D0
        BB(MGD)=(QCON+QRAD+QEVAP)*DLY1(MGD-1)*RD(IP)
     &          /CONDSOL(TPS,MTYPE)
        IF(IBOIL.EQ.1) BB(MGD)=0.D0
        CALL TDMA(AP,AE,AW,BB,TMIDS,MGD)
C      SET BACK THE TEMPERATURE INTO TP
        DO 22 I=1,MGD
        TPCD(IP,I)=TMIDS(I)
   22   CONTINUE
        DO 23 I=2,NGD
        TPCD(IP,I+MGD-1)=TMIDS(MGD)
   23   CONTINUE
        DO 24 I=2,LGD
        TPCD(IP,I+MGD+NGD-2)=TMIDS(MGD)
   24   CONTINUE
C      CHECK IF THE SURFACE IS TO BE MELTED
        IF(TPCD(IP,MGD).GT.TMELT) THEN
          WRITE(*,*)'NOTE, FROM (1,0,0) TO (1,1,0),TIME,IP=',TIME,IP
          SOLVL(IP)=1
          VI(IP)=MUI*(TMELT-TPCD(IP,MGD))
C      GET THE COMPUTATIONAL DOMAIN FOR LIQUID
C          DRM=VI(IP)*DTT
C      LIMIT THE TIMESTEP SO THAT THE INTERFACE WILL NOT MOVE TOO MUCH
C          IF(ABS(DRM).GT.DRMMIN) THEN
C            DRM=-DRMMIN*0.5D0
C          END IF
C      SET THE INITIAL TEMP OF LIQUID. FIX THE GRADIENT
c         TGRS=(TMIDS(MGD)-TMIDS(MGD-1))/(DLY1(MGD-1)*(RM(IP)+DRM))
c         DO 25 I=2,NGD
c         TPCD(IP,I+MGD-1)=TPCD(IP,I+MGD-2)+
c    &                     TGRS*DLY2(I-1)*(RD(IP)-RM(IP))
c  25     CONTINUE
C      SET THE INITIAL TEMP OF THE RESOLID. SET TO BE THE TEMP AT THE SURFACE
C          DO 26 I=2,LGD
C          TPCD(IP,I+MGD+NGD-2)=TPCD(IP,MGD+NGD-1)
C   26     CONTINUE
        END IF
C
C =========== (1,1,0)
      ELSE IF(SOLVS(IP).EQ.1.AND.SOLVL(IP).EQ.1.AND.
     &        SOLVRS(IP).EQ.0) THEN
C  --   PARTICLE IS IN THE STATUS OF LIQUID AND SOLID COEXISTANCE,
C       NO RESOLIDIFICATION
        DRM=VI(IP)*DTT
C      LIMIT THE TIMESTEP SO THAT THE INTERFACE WILL NOT MOVE TOO MUCH
        IF(ABS(DRM).GT.DRMMIN) THEN
          DRM=-DRMMIN*0.5D0
          DTT=DRM/VI(IP)
        END IF
        IF(DTT+DT0.GT.DT) THEN
          DTT=DT-DT0
          DRM=VI(IP)*DTT
        END IF
        RMTMP=RM(IP)+DRM
        IF(RMTMP.LE.0.D0) THEN
C    SHUT DOWN THE SOLID COMPUTATIONAL DOMAIN
          WRITE(*,*)'NOTE, CHANGE FROM (1,1,0) TO (0,1,0),TIME,IP=',
     &               TIME,IP
          SOLVS(IP)=0
          RM(IP)=0.D0
          GO TO 20
        END IF
C     SOLID
        ITVI=0
  2     CONTINUE
        ITVI=ITVI+1
        TMI=TMELT+VI(IP)/MUI
        IF(IMKINE.EQ.0) TMI=TMELT
        DRM=VI(IP)*DTT
        RMTMP=RM(IP)+DRM
        DO 31 I=1,MGD
        TMIDS(I)=TPCD(IP,I)
        SP(I)=0.D0
        SC(I)=0.D0
  31    CONTINUE
        CALL COEFF(AP,AE,AW,BB,MGD,TMIDS,Y1C,Y1B,SP,SC,
     &           DTT,RMTMP,0.D0,DRM,0.D0,NCS,1,MTYPE)
        AP(1)=1.D0
        AE(1)=1.D0
        AW(1)=0.D0
        BB(1)=0.D0
        AP(MGD)=1.D0
        AE(MGD)=0.D0
        AW(MGD)=0.D0
        BB(MGD)=TMI
        CALL TDMA(AP,AE,AW,BB,TMIDS,MGD)
C     LIQUID
        DO 32 I=1,NGD
        TMIDL(I)=TPCD(IP,I+MGD-1)
        SP(I)=0.D0
        SC(I)=0.D0
  32    CONTINUE
        TPS=TPCD(IP,NGD+MGD-1)
        QCON=HTCO*(TG-TPS)
        QRAD=STBOL*EMISS*(TAMB**4-TPS**4)
        QEVAP=DMPTC*LATVAP(TPS,MTYPE)/(4.0D0*PI*RD(IP)*RD(IP))
        CALL COEFF(AP,AE,AW,BB,NGD,TMIDL,Y2C,Y2B,SP,SC,
     &           DTT,RD(IP),RMTMP,DRD(IP),DRM,NCS,2,MTYPE)
        AP(1)=1.D0
        AE(1)=0.D0
        AW(1)=0.D0
        BB(1)=TMI
        AP(NGD)=1.D0
        AE(NGD)=0.D0
        AW(NGD)=1.D0
        BB(NGD)=(QCON+QRAD+QEVAP)*DLY2(NGD-1)*(RD(IP)-RMTMP)/
     &           CONDLIQ(TPS,MTYPE)
        QCON=QCON+QRAD+QEVAP
        IF(IBOIL.EQ.1) BB(NGD)=0.D0
        CALL TDMA(AP,AE,AW,BB,TMIDL,NGD)
C     CHECK IF THE INTERFACE MOVING VELOCITY CORRECT
        FLUXS=CONDSOL(TMI,MTYPE)*(TMIDS(MGD)-TMIDS(MGD-1))/
     &        (RMTMP*DLY1(MGD-1))
        FLUXL=CONDLIQ(TMI,MTYPE)*(TMIDL(2)-TMIDL(1))/
     &        ((RD(IP)-RMTMP)*DLY2(1))
        VIMID=(FLUXS-FLUXL)/(LATM*DENSSOL(TMI,MTYPE))
        DVI=VIMID-VI(IP)
c         WRITE(*,*)'FLUXS,FLUXL,VIMID,DVI,VI',
c     &            FLUXS,FLUXL,VIMID,DVI,VI(IP)
c          write(*,*)'*** Tg,Tps',tg,tps,tmids(ngd),tmids(ngd-1),tmelt
c          write(*,*)'+++ ITVI,Q',ITVI,qcon,qrad,qevap
c        if(ip.eq.1632) then
c       write(*,*)'drm,rm(ip),rd(ip),tmtmp',drm,rm(ip),rd(ip),rmtmp
c        write(*,*)'drm,rm,rmtmp,dtt,dt0=',
c    &      drm,rm(ip),rmtmp,dtt,dt0
c        write(*,*)'dvi,vi,vimid,fluxs,fluxl,ip,itvi=',
c    &      dvi,vi(ip),vimid,fluxs,fluxl,ip,itvi
c        write(*,*)'tmidl(2),(1),tmids(mgd),(mgd-1)=',
c    &      tmidl(2),tmidl(1),tmids(mgd),tmids(mgd-1),ip,itvi
c        write(*,*)'dvi/vi',dvi/vi(ip)
c        end if
        IF(VI(IP).EQ.0.D0) THEN
          VI(IP)=VI(IP)+RELAX*DVI
          GO TO 2
        ELSE
          IF(ABS(DVI/VI(IP)).GT.EPS) THEN
            VI(IP)=VI(IP)+RELAX*DVI
            IF(ITVI.LT.500) THEN
              GO TO 2
            ELSE
              WRITE(*,*)'*** VI NOT CONVT,1,',ITVI,ip
              write(*,*)'dvi,vi,vimid,fluxs,fluxl=',
     &                  dvi,vi(ip),vimid,fluxs,fluxl
c             STOP
            END IF
          ELSE
            IF(VI(IP).GT.0.D0) THEN
C   vi(ip)>0 indicates resolidify of particle
c             WRITE(*,*)'WARNING VI>0 ', VIMID,VI(IP),DVI,FLUXS,FLUXL
              IF(RM(IP).EQ.RD(IP)) THEN
c               WRITE(*,*)'NOTE, NO VI (1,1,0) TO (1,0,0),TIME=',TIME
c           write(*,*)'Tg,Tps,TPCD',tg,tps,TPCD(IP,mgd),TPCD(IP,mgd-1),
c     &                TPCD(IP,1)
                VI(IP)=0.D0
                SOLVL(IP)=0
                GOTO 3
              END IF
            ELSE
              RM(IP)=RMTMP
              DO 34 I=1,MGD
              TPCD(IP,I)=TMIDS(I)
   34         CONTINUE
              DO 35 I=2,NGD
              TPCD(IP,I+MGD-1)=TMIDL(I)
   35         CONTINUE
              DO 36 I=2,LGD
              TPCD(IP,I+MGD+NGD-2)=TMIDL(NGD)
   36         CONTINUE
C ---- SEE IF THE RESOLIDIFICATION TAKES PLACE
              IF(TMIDL(NGD).LT.TMELT) THEN
                WRITE(*,*)'NOTE, RESOLID (1,1,0) TO (1,1,1),TIME=',
     &                     TIME,IP
                SOLVRS(IP)=1
                VIRS(IP)=MUI*(TMIDL(NGD)-TMELT)
                DO 37 I=2,LGD
                TPCD(IP,I+MGD+NGD-2)=TMIDL(NGD)
   37           CONTINUE
              END IF
            END IF
          END IF
        END IF
C
C =========== (0,1,0)
      ELSE IF(SOLVS(IP).EQ.0.AND.SOLVL(IP).EQ.1.AND.
     &        SOLVRS(IP).EQ.0) THEN
C  --   PARTICLE IS TOTALLY MELTED
        DO 42 I=1,NGD
        TMIDL(I)=TPCD(IP,I+MGD-1)
        SP(I)=0.D0
        SC(I)=0.D0
  42    CONTINUE
        TPS=TPCD(IP,NMGD-2)
        QCON=HTCO*(TG-TPS)
        QRAD=STBOL*EMISS*(TAMB**4-TPS**4)
        QEVAP=DMPTC*LATVAP(TPS,MTYPE)/(4.0D0*PI*RD(IP)*RD(IP))
        CALL COEFF(AP,AE,AW,BB,NGD,TMIDL,Y2C,Y2B,SP,SC,
     &           DTT,RD(IP),0.D0,DRD(IP),0.D0,NCS,2,MTYPE)
        AP(1)=1.D0
        AE(1)=1.D0
        AW(1)=0.D0
        BB(1)=0.D0
        AP(NGD)=1.D0
        AE(NGD)=0.D0
        AW(NGD)=1.D0
        BB(NGD)=(QCON+QRAD+QEVAP)*DLY2(NGD-1)*(RD(IP)-RM(IP))/
     &          CONDLIQ(TPS,MTYPE)
        IF(IBOIL.EQ.1) BB(NGD)=0.D0
        CALL TDMA(AP,AE,AW,BB,TMIDL,NGD)
C   PUT BACK THE TEMPERATURE
        DO 44 I=1,MGD
        TPCD(IP,I)=TMIDL(1)
   44   CONTINUE
        DO 45 I=2,NGD
        TPCD(IP,I+MGD-1)=TMIDL(I)
   45   CONTINUE
        DO 46 I=2,LGD
        TPCD(IP,I+MGD+NGD-2)=TMIDL(NGD)
   46   CONTINUE
C  ----  SEE IF RESOLIDIFICATION TAKES PLACE
        IF(TMIDL(NGD).LT.TMELT) THEN
          WRITE(*,*)'NOTE, RESOLID(0,1,0) TO (0,1,1),TIME=',TIME,IP
          SOLVRS(IP)=1
          VIRS(IP)=MUI*(TMIDL(NGD)-TMELT)
        END IF
C
C =========== (0,1,1)
      ELSE IF(SOLVS(IP).EQ.0.AND.SOLVL(IP).EQ.1.AND.
     &        SOLVRS(IP).EQ.1) THEN
C  --   PARTICLE (0,1,1) IS PARTIALLY RESOLIDIFIED
        DRS=VIRS(IP)*DTT
C      LIMIT THE TIMESTEP SO THAT THE INTERFACE WILL NOT MOVE TOO MUCH
        IF(ABS(DRS).GT.DRMMIN) THEN
          DRS=-DRMMIN*0.5D0
          DTT=DRS/VIRS(IP)
        END IF
        IF(DTT+DT0.GT.DT) THEN
          DTT=DT-DT0
          DRS=VIRS(IP)*DTT
        END IF
        RSTMP=RRS(IP)+DRS
        IF(RSTMP.LE.0.D0) THEN
C    SHUT DOWN THE LIQUID COMPUTATIONAL DOMAIN
          WRITE(*,*)'NOTE, ALL RESOLID (0,1,1) TO (0,0,1),TIME=',TIME,IP
          SOLVL(IP)=0
          RRS(IP)=0.D0
          GO TO 20
        END IF
C     LIQUID
        ITVI=0
  5     CONTINUE
        ITVI=ITVI+1
        TMI=TMELT+VIRS(IP)/MUI
        IF(IMKINE.EQ.0) TMI=TMELT
        DRS=VIRS(IP)*DTT
        RSTMP=RRS(IP)+DRS
        DO 51 I=1,NGD
        TMIDL(I)=TPCD(IP,I+MGD-1)
        SP(I)=0.D0
        SC(I)=0.D0
  51    CONTINUE
        CALL COEFF(AP,AE,AW,BB,NGD,TMIDL,Y2C,Y2B,SP,SC,
     &           DTT,RSTMP,0.D0,DRS,0.D0,NCS,2,MTYPE)
        AP(1)=1.D0
        AE(1)=1.D0
        AW(1)=0.D0
        BB(1)=0.D0
        AP(NGD)=1.D0
        AE(NGD)=0.D0
        AW(NGD)=0.D0
        BB(NGD)=TMI
        CALL TDMA(AP,AE,AW,BB,TMIDL,NGD)
C     RESOLID
        DO 52 I=1,LGD
        TMIDRS(I)=TPCD(IP,I+MGD+NGD-2)
        SP(I)=0.D0
        SC(I)=0.D0
  52    CONTINUE
        TPS=TPCD(IP,NMGD-2)
        QCON=HTCO*(TG-TPS)
        QRAD=STBOL*EMISS*(TAMB**4-TPS**4)
        QEVAP=0.D0
        CALL COEFF(AP,AE,AW,BB,LGD,TMIDRS,Y3C,Y3B,SP,SC,
     &           DTT,RD(IP),RSTMP,DRD(IP),DRS,NCS,1,MTYPE)
        AP(1)=1.D0
        AE(1)=0.D0
        AW(1)=0.D0
        BB(1)=TMI
        AP(LGD)=1.D0
        AE(LGD)=0.D0
        AW(LGD)=1.D0
        BB(LGD)=(QCON+QRAD+QEVAP)*DLY3(LGD-1)*(RD(IP)-RSTMP)/
     &           CONDSOL(TPS,MTYPE)
        CALL TDMA(AP,AE,AW,BB,TMIDRS,LGD)
C     CHECK IF THE INTERFACE MOVING VELOCITY CORRECT
        FLUXS=CONDLIQ(TMI,MTYPE)*(TMIDL(NGD)-TMIDL(NGD-1))/
     &        (RSTMP*DLY2(NGD-1))
        FLUXL=CONDSOL(TMI,MTYPE)*(TMIDRS(2)-TMIDRS(1))/
     &        ((RD(IP)-RSTMP)*DLY3(1))
        VIMID=(FLUXS-FLUXL)/(LATM*DENSLIQ(TMI,MTYPE))
        DVI=VIMID-VIRS(IP)
        IF(VIRS(IP).EQ.0.D0) THEN
          VIRS(IP)=VIRS(IP)+RELAX*DVI
          GO TO 5
        ELSE
          IF(ABS(DVI/VIRS(IP)).GT.EPS) THEN
            VIRS(IP)=VIRS(IP)+RELAX*DVI
            IF(VIRS(IP).GT.0.D0) THEN
c              WRITE(*,*)'ERROR IN THE INTERFACE VELOCITY=',
c     &                 VIMID,VI(IP),DVI,FLUXS,FLUXL
c              STOP
            END IF
            IF(ITVI.LT.200) THEN
              GO TO 5
            ELSE
              WRITE(*,*)'*** ITERATION ON VI NOT CONVE,2,',
     &                   IP,ITVI
              write(*,*)'dvi,vi,vimid,fluxs,fluxl=',
     &                  dvi,virs(ip),vimid,fluxs,fluxl
c             STOP
            END IF
          ELSE
            RRS(IP)=RSTMP
            DO 54 I=1,NGD
            TPCD(IP,I+MGD-1)=TMIDL(I)
   54       CONTINUE
            DO 55 I=2,LGD
            TPCD(IP,I+MGD+NGD-2)=TMIDRS(I)
   55       CONTINUE
          END IF
        END IF
C
C =========== (0,0,1)
      ELSE IF(SOLVS(IP).EQ.0.AND.SOLVL(IP).EQ.0.AND.
     &        SOLVRS(IP).EQ.1) THEN
C  --   PARTICLE IS TOTALLY RESOLIDIFIED
        DO 62 I=1,LGD
        TMIDRS(I)=TPCD(IP,I+MGD+NGD-2)
        SP(I)=0.D0
        SC(I)=0.D0
  62    CONTINUE
        TPS=TPCD(IP,NMGD-2)
        QCON=HTCO*(TG-TPS)
        QRAD=STBOL*EMISS*(TAMB**4-TPS**4)
        QEVAP=0.D0
        CALL COEFF(AP,AE,AW,BB,LGD,TMIDRS,Y3C,Y3B,SP,SC,
     &           DTT,RD(IP),0.D0,DRD(IP),0.D0,NCS,1,MTYPE)
        AP(1)=1.D0
        AE(1)=1.D0
        AW(1)=0.D0
        BB(1)=0.D0
        AP(LGD)=1.D0
        AE(LGD)=0.D0
        AW(LGD)=1.D0
        BB(LGD)=(QCON+QRAD+QEVAP)*DLY3(LGD-1)*(RD(IP)-RRS(IP))/
     &          CONDSOL(TPS,MTYPE)
        CALL TDMA(AP,AE,AW,BB,TMIDRS,LGD)
C   PUT BACK THE TEMPERATURE
        DO 64 I=1,MGD
        TPCD(IP,I)=TMIDRS(1)
   64   CONTINUE
        DO 65 I=2,NGD
        TPCD(IP,I+MGD-1)=TMIDRS(1)
   65   CONTINUE
        DO 66 I=2,LGD
        TPCD(IP,I+MGD+NGD-2)=TMIDRS(I)
   66   CONTINUE
C
C =========== (1,1,1)
      ELSE IF(SOLVS(IP).EQ.1.AND.SOLVL(IP).EQ.1.AND.
     &        SOLVRS(IP).EQ.1) THEN
C  --   PARTICLE (0,1,1) IS PARTIALLY RESOLIDIFIED AND THERE ARE STILL
C       SOLID AND LIQUID REMAIN THERE, THREE DOMAINS AT THE SAME TIME
        DRM=VI(IP)*DTT
        DRS=VIRS(IP)*DTT
C      LIMIT THE TIMESTEP SO THAT THE INTERFACE WILL NOT MOVE TOO MUCH
        IF(ABS(DRM).GT.DRMMIN) THEN
          DRM=-DRMMIN*0.5D0
          DTT=DRM/VI(IP)
          DRS=VIRS(IP)*DTT
        END IF
        IF(ABS(DRS).GT.DRMMIN) THEN
          DRS=-DRMMIN*0.5D0
          DTT=DRS/VIRS(IP)
          DRM=VI(IP)*DTT
        END IF
        IF(DTT+DT0.GT.DT) THEN
          DTT=DT-DT0
          DRM=VI(IP)*DTT
          DRS=VIRS(IP)*DTT
        END IF
        RMTMP=RM(IP)+DRM
        RSTMP=RRS(IP)+DRS
        IF(RMTMP.LE.0.D0) THEN
C    SHUT DOWN THE SOLID COMPUTATIONAL DOMAIN
          WRITE(*,*)'NOTE, NO SOLID CORE (1,1,1) TO (0,1,1),TIME=',
     &               TIME,IP
          SOLVS(IP)=0
          RM(IP)=0.D0
          GO TO 20
        END IF
        IF(RSTMP.LE.RMTMP) THEN
C    SHUT DOWN THE LIQUID COMPUTATIONAL DOMAIN
          WRITE(*,*)'NOTE, NO MELT (1,1,1) TO (1,0,1),TIME=',TIME,IP
          SOLVL(IP)=0
          RRS(IP)=RM(IP)
          GO TO 20
        END IF
C     SOLID
        ITVI=0
  7     CONTINUE
        ITVI=ITVI+1
c        TMI=TMELT+VI(IP)/MUI
        TMI=TMELT-VI(IP)/MUI
        IF(IMKINE.EQ.0) TMI=TMELT
        DRM=VI(IP)*DTT
        RMTMP=RM(IP)+DRM
        DO 71 I=1,MGD
        TMIDS(I)=TPCD(IP,I)
        SP(I)=0.D0
        SC(I)=0.D0
  71    CONTINUE
        CALL COEFF(AP,AE,AW,BB,MGD,TMIDS,Y1C,Y1B,SP,SC,
     &           DTT,RMTMP,0.D0,DRM,0.D0,NCS,1,MTYPE)
        AP(1)=1.D0
        AE(1)=1.D0
        AW(1)=0.D0
        BB(1)=0.D0
        AP(MGD)=1.D0
        AE(MGD)=0.D0
        AW(MGD)=0.D0
        BB(MGD)=TMI
        CALL TDMA(AP,AE,AW,BB,TMIDS,MGD)
C     LIQUID
        DO 72 I=1,NGD
        TMIDL(I)=TPCD(IP,I+MGD-1)
        SP(I)=0.D0
        SC(I)=0.D0
  72    CONTINUE
        CALL COEFF(AP,AE,AW,BB,NGD,TMIDL,Y2C,Y2B,SP,SC,
     &           DTT,RSTMP,RMTMP,DRS,DRM,NCS,2,MTYPE)
        AP(1)=1.D0
        AE(1)=0.D0
        AW(1)=0.D0
        BB(1)=TMI
        AP(NGD)=1.D0
        AE(NGD)=0.D0
        AW(NGD)=0.D0
        BB(NGD)=TMI
        CALL TDMA(AP,AE,AW,BB,TMIDL,NGD)
C     CHECK IF THE INTERFACE MOVING VELOCITY CORRECT
        FLUXS=CONDSOL(TMI,MTYPE)*(TMIDS(MGD)-TMIDS(MGD-1))/
     &        (RMTMP*DLY1(MGD-1))
        FLUXL=CONDLIQ(TMI,MTYPE)*(TMIDL(2)-TMIDL(1))/
     &        ((RRS(IP)-RMTMP)*DLY2(1))
        VIMID=(FLUXS-FLUXL)/(LATM*DENSSOL(TMI,MTYPE))
        DVI=VIMID-VI(IP)
        IF(VI(IP).EQ.0.D0) THEN
          VI(IP)=VI(IP)+RELAX*DVI
          GO TO 7
        ELSE
          IF(ABS(DVI/VI(IP)).GT.EPS) THEN
            VI(IP)=VI(IP)+RELAX*DVI
            IF(VI(IP).GT.0.D0) THEN
c              WRITE(*,*)'ERROR IN THE INTERFACE VELOCITY=',
c     &                 VIMID,VI(IP),DVI,FLUXS,FLUXL
c              STOP
            END IF
            IF(ITVI.LT.500) THEN
              GO TO 7
            ELSE
              WRITE(*,*)'*** ITERATION ON VI NOT CONV,3,',
     &                   ITVI,IP
              write(*,*)'dvi,vi,vimid,fluxs,fluxl=',
     &                  dvi,virs(ip),vimid,fluxs,fluxl
c             STOP
            END IF
          ELSE
            RM(IP)=RMTMP
            DO 74 I=1,MGD
            TPCD(IP,I)=TMIDS(I)
   74       CONTINUE
            DO 75 I=2,NGD
            TPCD(IP,I+MGD-1)=TMIDL(I)
   75       CONTINUE
          END IF
        END IF
C  RESOLID + LIQUID
C     LIQUID
        ITVI=0
  8     CONTINUE
        ITVI=ITVI+1
        TMI=TMELT+VIRS(IP)/MUI
        IF(IMKINE.EQ.0) TMI=TMELT
        DRS=VIRS(IP)*DTT
        RSTMP=RRS(IP)+DRS
        DO 81 I=1,NGD
        TMIDL(I)=TPCD(IP,I+MGD-1)
        SP(I)=0.D0
        SC(I)=0.D0
  81    CONTINUE
        CALL COEFF(AP,AE,AW,BB,NGD,TMIDL,Y2C,Y2B,SP,SC,
     &           DTT,RSTMP,RM(IP),DRS,DRM,NCS,2,MTYPE)
        AP(1)=1.D0
        AE(1)=0.D0
        AW(1)=0.D0
        BB(1)=TMI
        AP(NGD)=1.D0
        AE(NGD)=0.D0
        AW(NGD)=0.D0
        BB(NGD)=TMI
        CALL TDMA(AP,AE,AW,BB,TMIDL,NGD)
C     RESOLID
        DO 82 I=1,LGD
        TMIDRS(I)=TPCD(IP,I+MGD+NGD-2)
        SP(I)=0.D0
        SC(I)=0.D0
  82    CONTINUE
        TPS=TPCD(IP,NMGD-2)
        QCON=HTCO*(TG-TPS)
        QRAD=STBOL*EMISS*(TAMB**4-TPS**4)
        QEVAP=0.D0
        CALL COEFF(AP,AE,AW,BB,LGD,TMIDRS,Y3C,Y3B,SP,SC,
     &           DTT,RD(IP),RSTMP,0.D0,DRS,NCS,1,MTYPE)
        AP(1)=1.D0
        AE(1)=0.D0
        AW(1)=0.D0
        BB(1)=TMI
        AP(LGD)=1.D0
        AE(LGD)=0.D0
        AW(LGD)=1.D0
        BB(LGD)=(QCON+QRAD+QEVAP)*DLY3(LGD-1)*(RD(IP)-RSTMP)/
     &           CONDSOL(TPS,MTYPE)
        CALL TDMA(AP,AE,AW,BB,TMIDRS,LGD)
C     CHECK IF THE INTERFACE MOVING VELOCITY CORRECT
        FLUXS=CONDLIQ(TMI,MTYPE)*(TMIDL(NGD)-TMIDL(NGD-1))/
     &        ((RSTMP-RM(IP))*DLY2(NGD-1))
        FLUXL=CONDSOL(TMI,MTYPE)*(TMIDRS(2)-TMIDRS(1))/
     &        ((RD(IP)-RSTMP)*DLY3(1))
        VIMID=-(FLUXS-FLUXL)/(LATM*DENSLIQ(TMI,MTYPE))
        DVI=VIMID-VIRS(IP)
        IF(VIRS(IP).EQ.0.D0) THEN
          VIRS(IP)=VIRS(IP)+RELAX*DVI
          GO TO 8
        ELSE
          IF(ABS(DVI/VIRS(IP)).GT.EPS) THEN
            VIRS(IP)=VIRS(IP)+RELAX*DVI
            IF(VIRS(IP).GT.0.D0) THEN
c              WRITE(*,*)'ERROR IN THE INTERFACE VELOCITY=',
c     &                 VIMID,VI(IP),DVI,FLUXS,FLUXL
c              STOP
            END IF
            IF(ITVI.LT.200) THEN
              GO TO 8
            ELSE
              WRITE(*,*)'*** ITERATION ON VI NOT CONT,4,',ITVI
              STOP
            END IF
          ELSE
            RRS(IP)=RSTMP
            DO 84 I=1,NGD
            TPCD(IP,I+MGD-1)=TMIDL(I)
   84       CONTINUE
            DO 85 I=2,LGD
            TPCD(IP,I+MGD+NGD-2)=TMIDRS(I)
   85       CONTINUE
          END IF
        END IF
C
C =========== (1,0,1)
      ELSE IF(SOLVS(IP).EQ.1.AND.SOLVL(IP).EQ.0.AND.
     &        SOLVRS(IP).EQ.1) THEN
C  --   PARTICLE IS SOLID PLUS RESOLID, NO LIQUID INSIDE
C  --   NO MELTING TAKES PLACE. AND NO VAPORIZATION WILL BE CONSIDERED
C        TMIDS(1)=TPCD(IP,1)
C        TMIDS(MGD)=TPCD(IP,MGD+NGD+LGD-2)
        DO 91 I=1,MGD
        YYC=Y1C(I)*RD(IP)
        IF(YYC.LE.RM(IP)) THEN
C  IN THE ORIGINAL SOLID PART (CAN BE REWRITE FOR SHORTER CPU TIME)
          DO 911 J=1,MGD-1
          RJ=Y1C(J)*RM(IP)
          RJ1=Y1C(J+1)*RM(IP)
          IF(YYC.GE.RJ.AND.YYC.LE.RJ1) THEN
            TMIDS(I)=TPCD(IP,J)+(TPCD(IP,J+1)-TPCD(IP,J))*
     &               (YYC-RJ)/(RJ1-RJ)
          END IF
 911      CONTINUE
        ELSE
          DO 912 J=1,LGD-1
          RJ=Y3C(J)*(RD(IP)-RM(IP))+RM(IP)
          RJ1=Y3C(J+1)*(RD(IP)-RM(IP))+RM(IP)
          IF(YYC.GE.RJ.AND.YYC.LE.RJ1) THEN
            TMIDS(I)=TPCD(IP,J+MGD+NGD-2)+
     &               (TPCD(IP,J+MGD+NGD-1)-TPCD(IP,J+MGD+NGD-2))*
     &               (YYC-RJ)/(RJ1-RJ)
          END IF
 912      CONTINUE
        END IF
        SP(I)=0.D0
        SC(I)=0.D0
  91    CONTINUE
        TPS=TPCD(IP,MGD+NGD+LGD-2)
        QCON=HTCO*(TG-TPS)
        QRAD=STBOL*EMISS*(TAMB**4-TPS**4)
        QEVAP=0.D0
        CALL COEFF(AP,AE,AW,BB,MGD,TMIDS,Y1C,Y1B,SP,SC,
     &           DTT,RD(IP),0.D0,0.D0,0.D0,NCS,1,MTYPE)
        AP(1)=1.D0
        AE(1)=1.D0
        AW(1)=0.D0
        BB(1)=0.D0
        AP(MGD)=1.D0
        AE(MGD)=0.D0
        AW(MGD)=1.D0
        BB(MGD)=(QCON+QRAD+QEVAP)*DLY1(MGD-1)*RD(IP)
     &          /CONDSOL(TPS,MTYPE)
        CALL TDMA(AP,AE,AW,BB,TMIDS,MGD)
C      SET BACK THE TEMPERATURE INTO TP
        DO 93 J=1,MGD
        YYC=Y1C(J)*RM(IP)
        DO 931 I=1,MGD-1
        RJ=Y1C(I)*RD(IP)
        RJ1=Y1C(I+1)*RD(IP)
        IF(YYC.GE.RJ.AND.YYC.LE.RJ1) THEN
        TPCD(IP,J)=TMIDS(I)+(TMIDS(I+1)-TMIDS(I))*
     &             (YYC-RJ)/(RJ1-RJ)
        END IF
 931    CONTINUE
   93   CONTINUE
        DO 94 J=MGD+1,MGD+NGD-1
        TPCD(IP,J)=TPCD(IP,MGD)
   94   CONTINUE
        DO 95 J=MGD+NGD,MGD+NGD+LGD-2
        YYC=Y3C(J)*(RD(IP)-RM(IP))+RM(IP)
        DO 951 I=1,MGD-1
        RJ=Y1C(I)*RD(IP)
        RJ1=Y1C(I+1)*RD(IP)
        IF(YYC.GE.RJ.AND.YYC.LE.RJ1) THEN
        TPCD(IP,J)=TMIDS(I)+(TMIDS(I+1)-TMIDS(I))*
     &             (YYC-RJ)/(RJ1-RJ)
        END IF
 951    CONTINUE
   95   CONTINUE
      END IF
C   CHECK IF THE ITERATION WITHIN PARTICLE IP IS FINISHED
      DT0=DT0+DTT
      IF(DT0.LT.DT) GO TO 3
  20  CONTINUE
C
C -----
      RETURN
      END
      SUBROUTINE SETUPPTC
C =======================================================================
C    THE INITIAL SETUP OF THE PARTICLE HEATING MODELING BY Y.P. WAN
C    12-10-97
C =======================================================================
      IMPLICIT NONE
      INCLUDE "COMMPTC.h"
C
C     NCS=0 FOR PLANE COORDINATE
C         1 FOR CYLINDRICAL CO.
C         2 FOR SPHERICAL CO.
C     PT1, PT2,PT3:   POWER INDEX USED FOR GRID FORMATION, CF. G.X. WANG 1997
C     NGD,MGD,LGD:    NUMBER OF GRID USED FOR LIQUID, SOLID AND RESOLID
C        REGION, RES.
      NCS=2
      PT1=1.0D0
      PT2=1.0D0
      PT3=1.0D0
      NGD=10
      MGD=10
      LGD=10
      NMGD=NGD+MGD+LGD
C  ----- SETTING UP THE GRID MESH
      CALL GRIDPTC
C
      RETURN
      END
      SUBROUTINE GRIDPTC
C =======================================================================
C    SETTING UP THE GRID SYSTEM FOR SOLID AND LIQUID COMPUTATIONAL DOMAIN
C =======================================================================
      IMPLICIT NONE
      INCLUDE "COMMPTC.h"
      INTEGER I
C   DIVIDE THE SOLID DOMAIN Y1=0-1 INTO MGD-2 ELEMENTS, MGD-1 FACES, M NODE
      DO 10 I=1,MGD-1
      Y1B(I)=(FLOAT(I-1)/FLOAT(MGD-2))**PT1
 10   CONTINUE
      Y1C(1)=0.D0
      Y1C(MGD)=Y1B(MGD-1)
C   -- CORRECT THE GEOMETRY FOR NODE 2 AND MGD-1
      DO 12 I=2,MGD-1
      Y1C(I)=Y1B(I-1)+0.5D0*(Y1B(I)-Y1B(I-1))
 12   CONTINUE
      DO 14 I=1,MGD-1
      DLY1(I)=Y1C(I+1)-Y1C(I)
 14   CONTINUE
C
C   DIVIDE THE LIQUID DOMAIN Y2=0-1 INTO NGD-2 ELEMENTS, NGD-1 FACES, N NODES
      DO 20 I=1,NGD-1
      Y2B(I)=(FLOAT(I-1)/FLOAT(NGD-2))**PT2
 20   CONTINUE
      Y2C(1)=0.D0
      Y2C(NGD)=Y2B(NGD-1)
C   -- CORRECT THE GEOMETRY FOR NODE 2 AND NGD-1
      DO 22 I=2,NGD-1
      Y2C(I)=Y2B(I-1)+0.5D0*(Y2B(I)-Y2B(I-1))
 22   CONTINUE
      DO 24 I=1,NGD-1
      DLY2(I)=Y2C(I+1)-Y2C(I)
 24   CONTINUE
C
C   DIVIDE THE RESOLIDIFICATION DOMAIN Y3=0-1 INTO LGD-2 ELEMENTS, LGD-1
C     FACES, LGD NODES
      DO 30 I=1,LGD-1
      Y3B(I)=(FLOAT(I-1)/FLOAT(LGD-2))**PT3
 30   CONTINUE
      Y3C(1)=0.D0
      Y3C(LGD)=Y3B(LGD-1)
C   -- CORRECT THE GEOMETRY FOR NODE 2 AND NGD-1
      DO 32 I=2,LGD-1
      Y3C(I)=Y3B(I-1)+0.5D0*(Y3B(I)-Y3B(I-1))
 32   CONTINUE
      DO 34 I=1,LGD-1
      DLY3(I)=Y3C(I+1)-Y3C(I)
 34   CONTINUE
C
C
      RETURN
      END
      SUBROUTINE COEFF(AP,AE,AW,BB,N,TOLD,YC,YB,SP,SC,
     &           DT,RPLUS,RMINUS,DRPLU,DRMIN,NCS,SOLID,MTYPE)
C =======================================================================
C    CALCULATE THE COEFFICIENCES FOR THE FINITE DIFFERENTIAL EQS
C    AP,AE,AW ARE THE COEF. OF Ti, Ti+1, and Ti-1, RESPECTIVELY
C    BB IS THE SOURE TERM
C    N IS THE TOTAL NUMBER OF GRID NODES
C    TOLD(N):  TEMP AT LAST TIME STEP
C    NCS    :  =0,1,2, FOR PLANE, CYLINDRICAL AND SPHERICAL COORDINATE
C    F(N)   :  FACTOR FOR COORDINATE SYSTEM
C    DV(N)  :  VOLUME OF THE CONTROL VOLUME ELEMENT
C    YC(N)  :  COORDINATE OF NODE
C    YB(N)  :  COORDINATE OF CONTROL VOLUME FACE
C    SP,SC  :  SOURCE TERM
C    DT     :  TIMESTEP
C    RPLUS  :  RIGHTMOST COORDINATE OF THE PHYSICAL DOMAIN, e.g., for
C              solid part, RPLUS=rm, RMINUS=0.0, for liquid part, RPLUS=rd,
C              RMINUS=rm.
C    RMINUS :  LEFTMOST COORDINATE OF THE PHYSICAL DOMAIN
C    DRPLU,DRMIN : THE CHANGE OF THE CORRESPONDING POSITION
C    SOLID  :  =1 FOR SOLID, =2 FOR LIQUID, MAINLY FOR THE DIFFERENT PROPERTIES
C =======================================================================
      IMPLICIT NONE
      INTEGER N,I,SOLID,NCS,MTYPE
      DOUBLE PRECISION AP(N),AE(N),AW(N),BB(N)
      DOUBLE PRECISION TOLD(N),YC(N),YB(N),SP(N),SC(N),
     &                 DT,DRPLU,DRMIN,RPLUS,RMINUS,RATIO,DR,DTDR2
      DOUBLE PRECISION AP0,KE,KW,RHO,CP,TI,TE,TW,DYE,DYW,
     &                 DPDYE,DPDYW,FE,FW,DVOL,FACE,FACW,DELTAE,DELTAW
      DOUBLE PRECISION CONDSOL,CONDLIQ,DENSSOL,DENSLIQ,CPSOL,CPLIQ
      EXTERNAL CONDSOL,CONDLIQ,DENSSOL,DENSLIQ,CPSOL,CPLIQ
C
C ----- SET UP GEOMETRY FACTOR
      DR=RPLUS-RMINUS
      RATIO=RMINUS/DR
      DTDR2=DT/(DR*DR)
C ----- THE BOUNDARY NODE I=1 AND N ARE CALCULATED OUTSIDE THIS SUBROUTINE
C       FOR THE SPECIAL NODE, I=2 AND N-1 COEFFICIENTS MUST BE SPECIALLY
C       TREATED
C ----- NODE 1,  ACCORDING ZERO GRADIENT
      AP(1)=1.D0
      AE(1)=1.D0
      AW(1)=0.D0
      BB(1)=0.D0
C ----- NODE N,  ACCORDING ZERO GRADIENT
      AP(N)=1.D0
      AE(N)=0.D0
      AW(N)=1.D0
      BB(N)=0.D0
C ----- NODE 2 TO N-1,  USING THE WHOLE CONTROL VOLUME AND LINEAR INTERPOLATION
      DO 10 I=2,N-1
      DYE=YC(I+1)-YC(I)
      DYW=YC(I)-YC(I-1)
      DPDYE=(YC(I+1)-YB(I))/DYE
      DPDYW=(YC(I)-YB(I-1))/DYW
      FACE=YB(I)+RATIO
      FACW=YB(I-1)+RATIO
      DELTAE=YB(I)*DRPLU+(1.0D0-YB(I))*DRMIN
      DELTAW=YB(I-1)*DRPLU+(1.0D0-YB(I-1))*DRMIN
      IF(NCS.EQ.0) THEN
        FE=1.D0
        FW=1.D0
        DVOL=YB(I)-YB(I-1)
      ELSE IF(NCS.EQ.1) THEN
        FE=FACE
        FW=FACW
        DVOL=(YB(I)-YB(I-1))*0.5D0*(FE+FW)
      ELSE IF(NCS.EQ.2) THEN
        FE=FACE*FACE
        FW=FACW*FACW
        DVOL=(YB(I)-YB(I-1))*1.D0/3.D0*
     &      (FE+FW+FACE*FACW)
      ELSE
        WRITE(*,*) 'ERROR WITH NCS, NO COORDINATE SYSTEM IDENTIFIED'
        STOP
      END IF
      TI=TOLD(I)
      TE=TOLD(I+1)*(1.D0-DPDYE)+TI*DPDYE
      TW=TI*(1.D0-DPDYW)+TOLD(I-1)*DPDYW
      IF(SOLID.EQ.1) THEN
        KE=CONDSOL(TE,MTYPE)
        KW=CONDSOL(TW,MTYPE)
        RHO=DENSSOL(TI,MTYPE)
        CP=CPSOL(TI,MTYPE)
      ELSE IF(SOLID.EQ.2) THEN
        KE=CONDLIQ(TE,MTYPE)
        KW=CONDLIQ(TW,MTYPE)
        RHO=DENSLIQ(TI,MTYPE)
        CP=CPLIQ(TI,MTYPE)
      END IF
      AP0=RHO*CP*DVOL
      AE(I)=DTDR2*KE*FE/DYE+
     &      RHO*CP/DR*DELTAE*FE*(1.D0-DPDYE)
      AW(I)=DTDR2*KW*FW/DYW-
     &      RHO*CP/DR*DELTAW*FW*DPDYW
      AP(I)=AE(I)+AW(I)+AP0-SP(I)*DVOL*DT
      BB(I)=SC(I)*DVOL*DT+AP0*TI
 10   CONTINUE
      RETURN
      END
      Subroutine tdma(a,b,c,d,T,n)
      Implicit none
      integer i,j,n
      double precision a,b,c,d,t,p,q
      dimension a(800),b(800),c(800),d(800),T(800)
      dimension p(800),q(800)
c------------------------------------------------------------------
c      From GuoXiang Wang, deleting Ti which is not used
c      Calculating the coefficient of recurrent relation Pj and Qj
c------------------------------------------------------------------
      p(1) = b(1)/a(1)
      q(1) = d(1)/a(1)
c
      do 111 i=2, n
      p(i) = b(i)/(a(i) -c(i)*p(i-1))
      q(i) = (d(i)+c(i)*q(i-1))/(a(i)-c(i)*p(i-1))
111      continue
c
c---------------------------------------------------
c      Successive back-substitution
c---------------------------------------------------
      T(n) = q(n)
      do 211 j=2,n
      i = n + 1 - j
      T(i) = p(i)*T(i+1) + q(i)
211      continue
c-----------------------------------------------------
c
      return
      end
      SUBROUTINE PRINTOUT(IP,TIME)
C =======================================================================
C    SETTING UP THE GRID SYSTEM FOR SOLID AND LIQUID COMPUTATIONAL DOMAIN
C =======================================================================
      IMPLICIT NONE
      INCLUDE "COMMPTC.h"
      INTEGER I,IP,IFILE,TIME
      DOUBLE PRECISION XCCD(NGRID)
      character*1 tab
      CHARACTER*6 CHATIME
      tab=char(9)
C   DIVIDE THE SOLID DOMAIN Y1=0-1 INTO MGD-2 ELEMENTS, MGD-1 FACES, M NODE
      IFILE=IP+10
      CALL CHARREAL(CHATIME,TIME)
      OPEN(IFILE,FILE='PT'//CHATIME//'.dout')
      DO 10 I=1,MGD
      XCCD(I)=RM(IP)*Y1C(I)
 10   CONTINUE
      DO 20 I=2,NGD
      XCCD(I+MGD-1)=((RRS(IP)-RM(IP))*Y2C(I)+RM(IP))
 20   CONTINUE
      DO 30 I=2,LGD
      XCCD(I+MGD+NGD-2)=((RD(IP)-RRS(IP))*Y3C(I)+RRS(IP))
 30   CONTINUE
      WRITE(IFILE,100)TAB
      DO 40 I=1,NMGD-2
      WRITE(IFILE,101)XCCD(I)*1.d6,TAB,TPCD(IP,I)
 40   CONTINUE
 100  FORMAT('*'/'r(mum)',a1,'Temp(K)')
 101  FORMAT(e13.5,1(a1,e13.5))
      RETURN
      END
      SUBROUTINE CHARREAL(CHATIME,REALTIME)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   convert a real value to a chracter string        C
C   which corresponds to this real value             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
      INTEGER MTIME,MTIME1,N1,I,K
      CHARACTER*6 CHATIME
      DOUBLE PRECISION REALTIME
C
      CHATIME(3:3)='_'
      MTIME=REALTIME*10000000
      MTIME1=MTIME/10
      N1=MTIME-MTIME1*10
      IF(N1.GE.5) MTIME1=MTIME1+1
      MTIME=MTIME1
      IF(MTIME.LT.1) THEN
        WRITE(6,*)'TIME IS TOO SMALL=',REALTIME
        STOP
      END IF
C
      DO 10 I=1,5
      K=7-I
      IF(I.GE.4) K=K-1
      MTIME1=MTIME/10
      N1=MTIME-MTIME1*10
      MTIME=MTIME1
      CHATIME(K:K)=CHAR(N1+48)
10    CONTINUE
      RETURN
      END
      SUBROUTINE PRINTPTC(IP)
C
C BY. Y.P.WAN 12-18-97
C ======================================================================
C
C      PRINTPTC IS USED TO WRITE OUT A DATA FILE FOR KALAIDERGRAPHE
C               THIS DATA FILE CONTAINS THE PARTICLE BEHAVIOR AT
C               DIFFERENT TIME
C      PRINTPTC CALLS THE FOLLOWING SUBROUTINES AND FUNCTIONS: NONE
C
C ======================================================================
C
C
      INCLUDE 'COML.h'
      INCLUDE 'COMMPTC.h'
      CHARACTER*1 TAB
      INTEGER KPRINT,KPRIPDF,KPFREQ,NPFREQ
      SAVE KPRINT,KPRIPDF,KPFREQ
      DATA KPRINT,KPRIPDF /0,0/
      DATA NPFREQ/200/
C
C <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
      TAB=CHAR(9)
      IF(KPRINT.EQ.0) THEN
        OPEN(07,FILE='PTCINFO.dout',STATUS='UNKNOWN',
     &                         FORM='FORMATTED')
        KPRINT=1
        KPFREQ=NPFREQ
        WRITE(07,100)TAB,TAB,TAB,TAB,TAB,TAB,TAB,TAB,TAB,TAB,TAB
      ELSE
c        OPEN(07,FILE='PTCINFO.dout',ACCESS='APPEND',
c     &                         FORM='FORMATTED')
      END IF
      IF(KPFREQ.EQ.NPFREQ) THEN
C
C ----- OUTPUT OF PTCINFO.dout
c        WRITE(07,102)TIME,TAB,XP(IP),TAB,YP(IP),TAB,ZP(IP),TAB,
c     &     UP(IP),TAB,VP(IP),TAB,WP(IP),TAB,TP(IP),TAB,TPCD(IP,1),TAB,
c     &     RM(IP)*1.D6,TAB,RADP(IP)*1.D4
        WRITE(07,102)TIME,TAB,XP(IP),TAB,YP(IP),TAB,ZP(IP),TAB,
     &   UP(IP),TAB,VP(IP),TAB,DMPTC(IP),TAB,TP(IP),TAB,TPCD(IP,1),TAB,
     &   RM(IP)*1.D6,TAB,RADP(IP)*1.D4,TAB,RRS(IP)*1.D6
        KPFREQ=0
      END IF
      KPFREQ=KPFREQ+1
100   FORMAT('*'/'Time(s)',A1,'XP(cm)',A1,'YP(cm)',A1,'ZP(cm)',
     &  A1,'UP(cm/s)',A1,'VP(cm/s)',A1,'dmvap(g/s)',
     &  A1,'TPS(K)',A1,'TPCNT(K)',A1,'Rm(mum)',A1,'Rd(mum)',
     &  A1,'RRs(mum)')
102   FORMAT(1PE10.3,11(A1,1PE10.3))
c      CLOSE (07)
      RETURN
C ================  END OF PRINTPTC
      END
      SUBROUTINE PRINTAMB(TIME,TG,TP,REYNUM,ANUSSL,SH,PRANUM,SCD,
     &     ROFILM,VISFLM,CPFILM,CONFLM,RHODIFF,HTCO,CD,
     &     CORRNC,RELVEL,QCON,QRAD,FACB)
C
C BY. Y.P.WAN 12-18-97
C ======================================================================
C
C      PRINTAMB IS USED TO WRITE OUT A DATA FILE FOR KALAIDERGRAPHE
C               THIS DATA FILE CONTAINS INFORMATION ABOUT THE
C               PLASMA PROPERTIES SURROUNDING PARTICLE IP
C      PRINTAMB CALLS THE FOLLOWING SUBROUTINES AND FUNCTIONS: NONE
C
C ======================================================================
C
C
      IMPLICIT NONE
      CHARACTER*1 TAB
      DOUBLE PRECISION TG,TP,REYNUM,ANUSSL,SH,PRANUM,SCD,
     &     ROFILM,VISFLM,CPFILM,CONFLM,RHODIFF,HTCO,CD,
     &     CORRNC,RELVEL,QCON,QRAD,TIME,FACB
      INTEGER KPRINT,KPRIPDF,KPFREQ,NPFREQ,IC
      SAVE KPRINT,KPRIPDF,KPFREQ
      DATA KPRINT,KPRIPDF /0,0/
      DATA NPFREQ/200/
C
C <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
      TAB=CHAR(9)
      IF(KPRINT.EQ.0) THEN
        OPEN(08,FILE='AMBINFO.dout',STATUS='UNKNOWN',
     &                         FORM='FORMATTED')
        KPRINT=1
        KPFREQ=NPFREQ
        WRITE(08,100)TAB,TAB,TAB,TAB,TAB,TAB,TAB,TAB,TAB,
     &               TAB,TAB,TAB,TAB,TAB,TAB,TAB,TAB,TAB,TAB
      ELSE
      END IF
      IF(KPFREQ.EQ.NPFREQ) THEN
C
C ----- OUTPUT OF AMBINFO.dout
        WRITE(08,102)TIME,TAB,TG,TAB,TP,TAB,REYNUM,TAB,ANUSSL,TAB,
     &     SH,TAB,PRANUM,TAB,SCD,TAB,ROFILM,TAB,VISFLM,TAB,CPFILM,TAB,
     &     CONFLM,TAB,RHODIFF,TAB,HTCO,TAB,CD,TAB,CORRNC,TAB,RELVEL,TAB,
     &     QCON,TAB,QRAD,TAB,FACB
        KPFREQ=0
      END IF
      KPFREQ=KPFREQ+1
100   FORMAT('*'/'Time(s)',A1,'Tg(K)',A1,'Tp(K)',A1,'ReNum',
     &  A1,'NuNum',A1,'ShNum',A1,'PrNum',
     &  A1,'ScNum',A1,'Rhof(g/cm3)',A1,'Muf()',A1,'Cp()',A1,'Cond()',
     &  A1,'RhoD()',A1,'Hcov()',A1,'Cd',A1,'CKnud',A1,'Vrel(cm/s)',
     &  A1,'QCON()',A1,'QRAD()',A1,'FACB')
102   FORMAT(1PE10.3,19(A1,1PE10.3))
c      CLOSE (07)
      RETURN
C ================  END OF PRINTPTC
      END
      FUNCTION CONDSOL(T,N)
C***********************************************************************
C
C------  THERMAL CONDUCTIVITY OF SOLID MO AS FUNCTION OF TEMPERATURE.
C         W/M.K
C        N=1 : NiCrAlY
C          2 : ZrO2
C  More detailed see data for ZrO2 from G.V. Samsonov, The Oxide handbook,
C      (Translation from Russian), Plenum, New York, 1973
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION CONDSOL,T
C
      IF(N.EQ.1) THEN
c  --- for Moly
c        CONDSOL=84.0D0
c  --- for Nickel
        CONDSOL=74.0D0
      ELSE IF(N.EQ.2) THEN
c  --- for ZrO2
        CONDSOL=2.0D0
      END IF
      RETURN
      END
      FUNCTION CONDLIQ(T,N)
C***********************************************************************
C
C------  THERMAL CONDUCTIVITY OF LIQUID MO AS FUNCTION OF TEMPERATURE.
C         W/M.K
C        N=1 : NiCrAlY
C          2 : ZrO2
C  More detailed see data for ZrO2 from G.V. Samsonov, The Oxide handbook,
C      (Translation from Russian), Plenum, New York, 1973
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION CONDLIQ,T
C
      IF(N.EQ.1) THEN
c  --- for Moly
c        CONDLIQ=46.0D0
c  --- for Nickel
        CONDLIQ=43.0D0
      ELSE IF(N.EQ.2) THEN
c  --- for ZrO2
        CONDLIQ=3.0D0
      END IF
      RETURN
      END
      FUNCTION DENSSOL(T,N)
C***********************************************************************
C
C------  DENSITY OF SOLID [KG/M3].
C        N=1 : NiCrAlY
C          2 : ZrO2
C
C from Samsonov, 1973, ZrO2 cubic system rho=6.27D3
C                           monclinic system rho=5.56D3
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION DENSSOL,T
C
      IF(N.EQ.1) THEN
c  --- for NiCrAlY
C        DENSSOL=8.11D3
c  --- for Nickel
        DENSSOL=8.90D3
      ELSE IF(N.EQ.2) THEN
c  --- for ZrO2
        DENSSOL=5.89D3
      END IF
      RETURN
      END
      FUNCTION DENSLIQ(T,N)
C***********************************************************************
C
C------  DENSITY OF LIQUID [KG/M3].
C        N=1 : NiCrAlY
C          2 : ZrO2
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION DENSLIQ,T
C
      IF(N.EQ.1) THEN
c  --- for NiCrAlY
C        DENSLIQ=8.11D3
c  --- for Nickel
        DENSLIQ=7.90D3
      ELSE IF(N.EQ.2) THEN
c  --- for ZrO2
        DENSLIQ=5.89D3
      END IF
      RETURN
      END
      FUNCTION CPSOL(T,N)
C***********************************************************************
C
C------  CP OF SOLID MO AS FUNCTION OF TEMPERATURE.
C         J/kg.K
C        N=1 : NiCrAlY
C          2 : ZrO2
C
C from Samsonov, Cp=a+bT-cT^-2 [J/kgmole.K]
C        IF(T.LT.1000.) THEN
C          CPSOL=6.95D2
C        ELSE IF(T.GE.1000.AND.T.LT.1420.) THEN
C          CPSOL=(57797.7+16678.8D-3*T)/107.2
C        ELSE IF(T.GE.1420.AND.T.LT.2500.) THEN
C          CPSOL=7.34D2
C        ELSE IF(T.GE.2500.AND.T.LT.2963.) THEN
C          CPSOL=6.95D2
C        ELSE IF(T.GE.2963.AND.T.LT.6000.) THEN
C          CPSOL=9.38D2
C        ELSE
C          WRITE(*,*)'*** TEMP OF ZrO2 OVER BOILING POINT',T
C        END IF
c
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION CPSOL,T
C
      IF(N.EQ.1) THEN
c  --- for Moly
C        CPSOL=4.48D2
c  --- for Nickel
        CPSOL=5.95D2
      ELSE IF(N.EQ.2) THEN
c  --- for ZrO2
        CPSOL=5.79D2
      END IF
      RETURN
      END
      FUNCTION CPLIQ(T,N)
C***********************************************************************
C
C------  CP OF LIQUID MO AS FUNCTION OF TEMPERATURE.
C         J/kg.K
C        N=1 : NiCrAlY
C          2 : ZrO2
C  More detailed see data for ZrO2 from G.V. Samsonov, The Oxide handbook,
C      (Translation from Russian), Plenum, New York, 1973
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION CPLIQ,T
C
      IF(N.EQ.1) THEN
c  --- for Moly
C        CPLIQ=4.48D2
c  --- for Nickel
        CPLIQ=6.20D2
      ELSE IF(N.EQ.2) THEN
c  --- for ZrO2
        CPLIQ=7.13D2
      END IF
      RETURN
      END
      FUNCTION LATVAP(T,N)
C***********************************************************************
C
C------  LATENT HEAT OF VAPORIZATION (AVERAGED VALUE NOT DEPENDENT ON T HERE).
C         J/kg
C        N=1 : NiCrAlY
C          2 : ZrO2
C  More detailed see data for ZrO2 from G.V. Samsonov, The Oxide handbook,
C      (Translation from Russian), Plenum, New York, 1973
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION LATVAP,T
C
      IF(N.EQ.1) THEN
c  --- for Moly
C        LATVAP=4.48D2
c  --- for Nickel
        LATVAP=7.33D6
      ELSE IF(N.EQ.2) THEN
c  --- for ZrO2
        LATVAP=6.0D6
      END IF
      RETURN
      END
      FUNCTION MUINTFACE(N)
C***********************************************************************
C
C------  CALCULATION OF THE LINEAR KINETICS COEFFICIENT OF THE MELTING
C        INTERFACE M/S.K
C        N=1 : NiCrAlY
C          2 : ZrO2
C
C    WHERE CAN ONE GET SUCH DATA
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION MUINTFACE
C
      IF(N.EQ.1) THEN
c  --- for Moly
C        MUINTFACE=0.23
c  --- for Nickel
        MUINTFACE=0.85
      ELSE IF(N.EQ.2) THEN
c  --- for ZrO2
        MUINTFACE=1.d-2
      END IF
      RETURN
      END
      FUNCTION PVAPOR(T,NSP)
C***********************************************************************
C
C------  CALCULATION OF VAPOR PRESSURE AT THE DROP SURFACE, T IS DROP
C        SURFACE TEMPERATURE.   [N/M^2]
C        WE ARE LACKING THE VAPOR DATA FOR ALLOY, using Mo's data
C        NSP=1 : NiCrAlY
C            2 : ZrO2
C       Data for ZrO2 from G.V. Samsonov, The Oxide handbook,
C      (Translation from Russian), Plenum, New York, 1973
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER NI1,NI2,NLO,NHI,N,ILOW,NSP
      PARAMETER (NI1=24,NI2=25)
      DOUBLE PRECISION TE1,PV1,TE2,PV2,T,TC,PVAPOR,DDC,DRAT
      DIMENSION TE1(NI1),PV1(NI1),TE2(NI2),PV2(NI2)
c  --- for Moly
C      DATA TE1/298.15,500.,1000.,1500.,1944.,
C     &        2000.,2063.,2199.,2200.,2354.,
C     &        2400.,2533.,2600.,2744.,2800.,
C     &        2890.,3000.,3200.,3325.,3400.,
C     &        3600.,3728.,3800.,4000.,4241.,
C     &        4500.,4912.,5000./
C      DATA PV1/5.D-108,1.7D-61,1.3D-27,7.2D-16,1.D-10,
C     &        3.1D-10,1.D-9,1.D-8,1.02D-8,1.D-7,
C     &        1.87D-7,1.D-6,2.16D-6,1.D-5,1.73D-5,
C     &        4.00D-5,9.88D-5,4.34D-4,1.D-3,1.6D-3,
C     &        5.09D-3,1.D-2,1.43D-2,3.65D-2,1.D-1,
C     &        0.263,1.0,1.30/
C
c  --- for Nickel
      DATA TE1/298.15,500.,1000.,1200.,1400.,
     &        1500.,1600.,1700.,1726.,1800.,
     &        1900.,2000.,2100.,2200.,2300.,
     &        2400.,2500.,2600.,2700.,2800.,
     &        2900.,3000.,3187.,3200./
      DATA PV1/3.9D-68,9.8D-38,2.1D-15,1.1D-11,4.96D-9,
     &        5.03D-8,4.12D-7,2.63D-6,4.08D-6,1.28D-5,
     &        5.20D-5,1.82D-4,5.60D-4,1.55D-3,3.92D-3,
     &        9.11D-3,1.97D-2,4.00D-2,7.70D-2,0.141,
     &        0.246,0.414,1.00,1.06/
C
c  --- for ZrO2
      DATA TE2/20.,527.,1027.,1527.,2027.,
     &        2227.,2427.,2527.,2677.,2827.,
     &        3027.,3227.,3427.,3627.,3827.,
     &        4027.,4127.,4227.,4327.,4427.,
     &        4527.,4627.,4827.,5027.,5227./
      DATA PV2/4.78D-118,2.02D-34,5.1D-16,5.65D-8,1.74D-3,
     &        3.35D-2,4.12D-1,1.26D0,5.82D0,19.4D0,
     &        8.01D1,2.77D2,8.26D2,2.17D3,5.14D3,
     &        1.11D4,1.58D4,2.22D4,3.06D4,4.15D4,
     &        5.55D4,7.31D4,1.22D5,1.96D5,3.00D5/
C
      TC=T
      IF(NSP.EQ.1) THEN
C  FOR MATERIALS 1, NICRALY -- USING DATA OF MOLY
        IF(TC.LT.TE1(1)) THEN
          WRITE(6,*)'*** PVAPOR, TEMP OF PARTICLE TOO LOW=',T
          STOP
        END IF
        IF(TC.GT.TE1(NI1)) THEN
          WRITE(6,*)'*** PVAPOR, TEMP. OF PARTICLE TOO HIGH',T
c          STOP
          tc=TE1(NI1)
        END IF
c ---  with bisect interval method
        nlo=1
        nhi=NI1
3718    n=(nlo+nhi)/2
        if(TC-TE1(n)) 3720,3750,3740
c  -- TC is .lt. TE1(n)
3720    nhi=n
        if(TC-TE1(n-1)) 3718,3730,3730
c  -- TC is in (n,n-1)
3730    ilow=n-1
        goto 3760
c  --  TC is .gt. TE1(n)
3740    nlo=n
        if(TC-TE1(n+1)) 3750,3750,3718
c  --  TC is in (n,n+1)
3750    ilow=n
3760    continue
        ddc=TE1(ilow+1)-TE1(ilow)
        if(ddc.lt.1.d-40) then
          drat=0.d0
        else
          drat=(TC-TE1(ilow))/ddc
        end if
        PVAPOR=PV1(ilow)+drat*(PV1(ilow+1)-PV1(ilow))
c ---- convert atm to pasica
        PVAPOR=PVAPOR/0.9869D-5
      END IF
C
      IF(NSP.EQ.2) THEN
C  FOR MATERIALS 2, ZrO2 --
        TC=TC-273.0
        IF(TC.LT.TE2(1)) THEN
          WRITE(6,*)'*** PVAPOR, TEMP OF PARTICLE TOO LOW=',T
          STOP
        END IF
        IF(TC.GT.TE2(NI2)) THEN
          WRITE(6,*)'*** PVAPOR, TEMP. OF PARTICLE TOO HIGH',T
c         STOP
          tc=TE2(NI2)
        END IF
c ---  with bisect interval method
        nlo=1
        nhi=NI2
4718    n=(nlo+nhi)/2
        if(TC-TE2(n)) 4720,4750,4740
c  -- TC is .lt. TE2(n)
4720    nhi=n
        if(TC-TE2(n-1)) 4718,4730,4730
c  -- TC is in (n,n-1)
4730    ilow=n-1
        goto 4760
c  --  TC is .gt. TE2(n)
4740    nlo=n
        if(TC-TE2(n+1)) 4750,4750,4718
c  --  TC is in (n,n+1)
4750    ilow=n
4760    continue
        ddc=TE2(ilow+1)-TE2(ilow)
        if(ddc.lt.1.d-40) then
          drat=0.d0
        else
          drat=(TC-TE2(ilow))/ddc
        end if
        PVAPOR=PV2(ilow)+drat*(PV2(ilow+1)-PV2(ilow))
      END IF
      RETURN
      END
      FUNCTION WPOWDER(N)
C***********************************************************************
C
C------  MOLECULAR WEIGHT OF POWDER MAERIALS
C        [G/MOL]
C        NSP=1 : NiCrAlY
C            2 : ZrO2
C
C***********************************************************************
C
      IMPLICIT NONE
      DOUBLE PRECISION WPOWDER
      INTEGER N
      IF(N.EQ.1) THEN
c  --- for NiCrAlY
C        WPOWDER=226.6D0
c  --- for NICKEL
        WPOWDER=58.69D0
      ELSE IF(N.EQ.2) THEN
c  --- for ZrO2
        WPOWDER=107.2D0
      END IF
      RETURN
      END
      FUNCTION ROPTC(N)
C***********************************************************************
C
C------  DENSITY  OF POWDER MAERIALS
C        [KG/M^3]
C        NSP=1 : NiCrAlY
C            2 : ZrO2
C
C***********************************************************************
C
      IMPLICIT NONE
      DOUBLE PRECISION ROPTC
      INTEGER N
C
      IF(N.EQ.1) THEN
c  --- for NiCrAlY
C        ROPTC=8.11D3
c  --- for NICKEL
        ROPTC=7.90D3
      ELSE IF(N.EQ.2) THEN
c  --- for ZrO2
        ROPTC=5.89D3
      END IF
      RETURN
      END

