*DECK USERPP
      SUBROUTINE USERPP
C==============================================================USERPP
C     THIS SUBROUTINE PROVIDES PARTICLE DRAG AND HEAT TRANSFER
C     COEFFICIENTS. REQUIRES RELVEL FROM SOUBROUTINE PINT
C----------
C     CALLED BY PINT 
C----------
C     INPUT: RELVEL
C     OUTPUT: CD,HTCO
C==============================================================USERPP
      INCLUDE 'COML.h'
      INCLUDE 'COMMPTC.h'
C----------
      DOUBLE PRECISION TB,FR,FR1,RCTOT,AMWA,RROT
      DOUBLE PRECISION ROFILM,TFILM,CPFILM,VISFLM,CONFLM
      DOUBLE PRECISION GAMSUR,VISURF,CPSURF,ROSURF,VISFRE,CPFREE
      DOUBLE PRECISION REYNUM,PRANUM,ANUSSL,VISKRA,CDF
      DOUBLE PRECISION CORRNC,CORRF,CDC1,CDC2,CDC3
C----------
      INTEGER N
C----------
      double precision Airtic,Airtiy,cArb,cHb,cNb,cArp,cHp,cAirp,
     1     yArb,yHb,yNb,yArp,yHp,yAirp
C==============================================================USERPP
C     PLASMA PROPERTIES FOR PARTICLE-PLASMA INTERACTIONS.
C     LAMINAR VISCOSITIES AND TOTAL THERMAL CONDUCTIVITIES INCLUDING
C     REACTIONAL PART.
C--------------------------------------------------------------USERPP
      DOUBLE PRECISION VISCPP(21,NCMP),CONDPP(21,NCMP),CSUBPP(21,NCMP)
      INTEGER IKIVA
      DOUBLE PRECISION SCD,RHODIFF,SH,BNUMB,LN1PB,PVAP,
     &       PVAPOR,WPOWDER,WMETAL,YMS,YMIC,RADPTMP,HE,HM,HMIX,
     &       QCON,QRAD,DMHEAT,DMDTTMP,PGDPV,ROPTC,FACA,FACB,LATVAP
      EXTERNAL PVAPOR,WPOWDER,ROPTC,LATVAP
C--------------------------------------------------------------USERPP
C     ARGON-HELIUM (100:47) PLASMA TRANSPORT PROPERTY DATA,
C     OBTAINED FROM THE UNIVERSITY OF MINNESOTA
C     INTERVALS ARE T=1000(N-1). UNITS ARE IN CGS.
C--------------------------------------------------------------USERPP
c     DATA (VISCPP(N,1),N=1,21)
c    1  /5.00D-04, 1.59D-03, 2.47D-03, 3.13D-03, 3.70D-03,
c    2   4.22D-03, 4.70D-03, 5.15D-03, 5.57D-03, 5.96D-03,
c    3   6.24D-03, 6.28D-03, 5.97D-03, 5.44D-03, 4.97D-03,
c    4   4.69D-03, 4.58D-03, 4.58D-03, 4.59D-03, 4.55D-03, 4.37D-03/
C---------- INCLUDING REACTIONAL PART
c     DATA (CONDPP(N,1),N=1,21)
c    1  /4.00D+04, 6.90D+04, 9.94D+04, 1.22D+05, 1.42D+05,
c    2   1.60D+05, 1.76D+05, 1.91D+05, 2.10D+05, 2.36D+05,
c    3   2.77D+05, 3.35D+05, 4.12D+05, 5.03D+05, 5.83D+05,
c    4   6.19D+05, 6.17D+05, 6.09D+05, 6.15D+05, 6.35D+05, 6.65D+05/
C---------- INCLUDING REACTION ENERGY
c     DATA (CSUBPP(N,1),N=1,21)
c    1  /7.30D+06, 7.30D+06, 7.30D+06, 7.30D+06, 7.30D+06,
c    2   7.31D+06, 7.33D+06, 7.56D+06, 8.55D+06, 1.16D+07,
c    3   1.85D+07, 3.17D+07, 5.25D+07, 7.76D+07, 9.43D+07,
c    4   8.85D+07, 6.56D+07, 4.48D+07, 3.57D+07, 3.89D+07, 5.37D+07/
C==============================================================USERPP
C     PROPERTIES OF ARGON PLASMA
C--------------------------------------------------------------USERPP
      DATA (VISCPP(N,1),N=1,21)
     1  /6.91D-05, 5.35D-04, 8.07D-04, 1.13D-03, 1.46D-03,
     2   1.82D-03, 2.06D-03, 2.29D-03, 2.51D-03, 2.69D-03,
     3   2.80D-03, 2.71D-03, 2.31D-03, 1.69D-03, 1.09D-03,
     4   6.66D-04, 4.07D-04, 2.97D-04, 2.51D-04, 2.36D-04, 2.33D-04/
C---------- INCLUDING REACTIONAL PART
      DATA (CONDPP(N,1),N=1,21)
     1  /5.21D+02, 4.27D+03, 8.92D+03, 9.38D+03, 1.18D+04,
     2   1.43D+04, 1.85D+04, 1.97D+04, 2.58D+04, 3.83D+04,
     3   5.88D+04, 9.13D+04, 1.48D+05, 1.72D+05, 2.08D+05,
     4   2.20D+05, 2.12D+05, 2.06D+05, 2.11D+05, 2.25D+05, 2.44D+05/
C---------- INCLUDING REACTION ENERGY
      DATA (CSUBPP(N,1),N=1,21)
     1  /5.20D+06, 5.20D+06, 5.20D+06, 5.20D+06, 5.20D+06,
     2   5.20D+06, 5.21D+06, 5.42D+06, 6.33D+06, 9.02D+06,
     3   1.55D+07, 2.73D+07, 4.69D+07, 7.23D+07, 9.39D+07,
     4   9.31D+07, 7.07D+07, 4.63D+07, 3.11D+07, 2.55D+07, 2.78D+07/
C==============================================================USERPP
C     PROPERTIES OF HYDROGEN PLASMA
C--------------------------------------------------------------USERPP
      DATA (VISCPP(N,2),N=1,21)
     1  /6.80D-05, 1.94D-04, 3.21D-04, 4.50D-04, 4.54D-04,
     2   5.20D-04, 6.25D-04, 7.33D-04, 8.40D-04, 9.37D-04,
     3   1.01D-03, 1.03D-03, 9.89D-04, 8.78D-04, 7.17D-04,
     4   5.42D-04, 3.85D-04, 2.68D-04, 1.93D-04, 1.44D-04, 1.14D-04/
C---------- INCLUDING REACTIONAL PART
      DATA (CONDPP(N,2),N=1,21)
     1  /1.40D+04, 3.85D+04, 9.11D+04, 1.06D+06, 9.05D+05,
     2   2.34D+05, 2.04D+05, 2.39D+05, 3.04D+05, 4.34D+05,
     3   6.90D+05, 1.15D+06, 1.90D+06, 2.95D+06, 4.12D+06,
     4   4.90D+06, 4.84D+06, 3.96D+06, 2.81D+06, 1.93D+06, 1.34D+06/
C---------- INCLUDING REACTION ENERGY
      DATA (CSUBPP(N,2),N=1,21)
     1  /1.40D+08, 1.50D+08, 5.64D+08, 1.11D+09, 2.13D+08,
     2   2.16D+08, 2.70D+08, 5.92D+08, 1.67D+09, 4.02D+09,
     3   5.62D+09, 3.22D+09, 1.27D+09, 5.58D+08, 4.91D+08,
     4   4.41D+08, 4.24D+08, 4.18D+08, 4.15D+08, 4.14D+08, 4.13D+08/
C==============================================================USERPP
C     PROPERTIES OF NITROGEN PLASMA
C--------------------------------------------------------------USERPP
c     DATA (VISCPP(N,2),N=1,21)
c    1  /1.06D-04, 4.06D-04, 6.78D-04, 9.23D-04, 1.15D-03,
c    2   1.38D-03, 1.65D-03, 1.96D-03, 2.17D-03, 2.32D-03,
c    3   2.40D-03, 2.30D-03, 1.96D-03, 1.47D-03, 9.90D-04,
c    4   6.28D-04, 3.99D-04, 2.69D-04, 2.02D-04, 1.68D-04, 1.53D-04/
C---------- INCLUDING REACTIONAL PART
c     DATA (CONDPP(N,2),N=1,21)
c    1  /1.07D+03, 5.00D+03, 1.20D+04, 2.30D+04, 4.50D+04,
c    2   9.50D+04, 2.66D+05, 4.90D+05, 1.77D+05, 1.29D+05,
c    3   1.18D+05, 1.26D+05, 1.57D+05, 2.10D+05, 2.46D+05,
c    4   2.48D+05, 2.36D+05, 2.18D+05, 2.10D+05, 2.17D+05, 2.34D+05/
C---------- INCLUDING REACTION ENERGY
c     DATA (CSUBPP(N,2),N=1,21)
c    1  /2.83D+07, 2.83D+07, 2.83D+07, 2.83D+07, 2.83D+07,
c    2   2.83D+07, 8.56D+07, 1.76D+08, 1.11D+08, 5.36D+07,
c    3   5.55D+07, 8.03D+07, 1.24D+08, 1.76D+08, 2.19D+08,
c    4   2.24D+08, 1.88D+08, 1.37D+08, 9.63D+07, 7.09D+07, 5.80D+07/
C==============================================================USERPP
C     PROPERTIES OF AIR PLASMA
C--------------------------------------------------------------USERPP
      DATA (VISCPP(N,3),N=1,21)
     1  /1.06D-04, 4.00D-04, 6.60D-04, 8.80D-04, 1.06D-03,
     2   1.30D-03, 1.64D-03, 2.00D-03, 2.30D-03, 2.56D-03,
     3   2.68D-03, 2.64D-03, 2.48D-03, 2.18D-03, 1.76D-03,
     4   1.22D-03, 5.60D-04, 2.69D-04, 2.02D-04, 1.68D-04, 1.53D-04/
C---------- INCLUDING REACTIONAL PART
      DATA (CONDPP(N,3),N=1,21)
     1  /1.00D+03, 3.18D+03, 1.06D+04, 3.77D+04, 4.97D+04,
     2   6.44D+04, 2.27D+05, 4.03D+05, 2.34D+05, 1.47D+05,
     3   1.76D+05, 2.51D+05, 3.61D+05, 4.97D+05, 6.28D+05,
     4   6.91D+05, 6.54D+05, 5.65D+05, 4.86D+05, 4.42D+05, 4.29D+05/
C---------- INCLUDING REACTION ENERGY
      DATA (CSUBPP(N,3),N=1,21)
     1  /1.14D+07, 1.14D+07, 1.46D+07, 2.79D+07, 3.40D+07,
     2   2.76D+07, 7.60D+07, 1.41D+08, 8.44D+07, 4.41D+07,
     3   4.93D+07, 7.38D+07, 1.12D+08, 1.61D+08, 2.04D+08,
     4   2.17D+08, 1.89D+08, 1.43D+08, 1.02D+08, 7.42D+07, 6.03D+07/
C==============================================================USERPP
      DATA CDC1,CDC2,CDC3 /24.0D0,0.4D0,0.45D0/
C--------------------------------------------------------------USERPP
C     CORRNC = CORRECTION FACTOR FOR NON-CONTINUUM EFFECT = 1/(1+AA)
C     AA=((2-A)/A)*(GAMMA/(1+GAMMA))*(4/PRS)*(KNUDSEN NUMBER)
C       =CORRF*((9*GAMMA-5)/(1+GAMMA))*SQRT(TP/MW)*VISURF/(P*RADP)
C     A = THERMAL ACCOMODATION COEFFICIENT = 0.8
C     GAMMA = SPECIFIC HEAT RATIO
C     PRS = PRANDTL NUMBER AT PARTICLE SURFACE = 4*GAMMA/(9*GAMMA-5)
C     CORRF = 0 IN ORDER TO NEGLECT NON-CONTINUUM EFFECT
C--------------------------------------------------------------USERPP
      DATA CORRF /8.571D3/
C ---- INTRODUCING SCHMIDT NUMBER BY Y.P. WAN 12-16-97
C      IEVAP=0, NO EVAP., 1 EVAP.   
C      IKIVA, INDEX FOR THE SELECTION OF VAPORIZATION MODEL
C           =1  USE THE MODEL USED IN KIVA
C            0                        THE PAPER OF WESTHOFF 1992
      DATA SCD/0.7D0/
      DATA IEVAP,IKIVA/1,1/
C--------------------------------------------------------------USERPP
c     data Airtic,Airtiy /3.761904762d0,3.291666667d0/
C==============================================================USERPP
      DO 200 IP=1,NP
      IC=ICP(IP)  
C--------------------------------------------------------------USERPP
c---- using argon, hydrogen, and air properties and mole fractions
c---- ignore the differences of N-O ratio than air
      cArb=spd(ic,1)*rmw(1)+spd(ic,2)*rmw(2)
      cHb=spd(ic,3)*rmw(3)+half*(spd(ic,4)*rmw(4)+spd(ic,5)*rmw(5))
      cNb=spd(ic,6)*rmw(6)+spd(ic,7)*rmw(7)
     1    +half*(spd(ic,8)*rmw(8)+spd(ic,9)*rmw(9))
      rctot=one/(cArb+cHb+cNb)
c----- molar frcations
      cArp=cArb*rctot
      cHp=cHb*rctot
      cAirp=cNb*rctot
c----- mass fractions
      yArb=spd(ic,1)+spd(ic,2)
      yHb=spd(ic,3)+spd(ic,4)+spd(ic,5)
      yNb=spd(ic,6)+spd(ic,7)+spd(ic,8)+spd(ic,9)
      rrot=one/(yArb+yHb+yNb)
      yArp=yArb*rrot
      yHp=yHb*rrot
      yAirp=yNb*rrot
c----- average molecular weight at Tsurf -- may not be correct at Tfilm
      amwa=(MW(1)*cArb+MW(3)*cHb+MW(6)*cNb)*rctot
C========================================= FILM PROPERTIES=====USERPP
C     NOTE THAT TRANSPORT PROPERTIES ARE LAMINAR PROPERTIES.
C        (SEE RINPUT TO SEE IF TRANSLATIONAL PART IS USED)
C     CONFLM = THERMAL CONDUCTIVITY AT FILM TEMPERATURE
C     VISFLM = VISCOSITY
C     CPFILM = SPECIFIC HEAT INCLUDING REACTIONAL PART
C--------------------------------------------------------------USERPP
      TFILM=HALF*(TEMP(IC)+TP(IP))
      TB=THOUTH*TFILM
      IT=INT(TB)
      IT=MIN(IT,19)
      FR=TB-DBLE(IT)
      FR1=ONE-FR
        CONFLM=cArp*(FR1*CONDPP(IT+1,1)+FR*CONDPP(IT+2,1))
     1         +cHp*(FR1*CONDPP(IT+1,2)+FR*CONDPP(IT+2,2))
     2       +cAirp*(FR1*CONDPP(IT+1,3)+FR*CONDPP(IT+2,3))
        VISFLM=cArp*(FR1*VISCPP(IT+1,1)+FR*VISCPP(IT+2,1))
     1         +cHp*(FR1*VISCPP(IT+1,2)+FR*VISCPP(IT+2,2))
     2       +cAirp*(FR1*VISCPP(IT+1,3)+FR*VISCPP(IT+2,3))
C--------------------------------------------------------------USERPP
C     SPECIFIC HEAT IS BASED ON THE DENSITY FRACTIONS.
C--------------------------------------------------------------USERPP
        CPFILM=yArp*(FR1*CSUBPP(IT+1,1)+FR*CSUBPP(IT+2,1))
     1         +yHp*(FR1*CSUBPP(IT+1,2)+FR*CSUBPP(IT+2,2))
     2       +yAirp*(FR1*CSUBPP(IT+1,3)+FR*CSUBPP(IT+2,3))
C--------------------------------------------------------------USERPP
C     ROFILM ARE CALCULATED FROM IDEAL GAS EQUATION OF STATE ASSUMING
C          THAT THERE IS NO DISSOCIATION/IONIZATION.
C--------------------------------------------------------------USERPP
      ROFILM=AMWA*P(IC)/(RGAS*TFILM)
C================================== FREE STREAM PROPERTIES=====USERPP
      TB=THOUTH*TEMP(IC)
      IT=INT(TB)
      IT=MIN(IT,19)
      FR=TB-DBLE(IT)
      FR1=ONE-FR
        VISFRE=cArp*(FR1*VISCPP(IT+1,1)+FR*VISCPP(IT+2,1))
     1         +cHp*(FR1*VISCPP(IT+1,2)+FR*VISCPP(IT+2,2))
     2       +cAirp*(FR1*VISCPP(IT+1,3)+FR*VISCPP(IT+2,3))
        CPFREE=yArp*(FR1*CSUBPP(IT+1,1)+FR*CSUBPP(IT+2,1))
     1         +yHp*(FR1*CSUBPP(IT+1,2)+FR*CSUBPP(IT+2,2))
     2       +yAirp*(FR1*CSUBPP(IT+1,3)+FR*CSUBPP(IT+2,3))
C============================= PARTICLE SURFACE PROPERTIES=====USERPP
      TB=THOUTH*TP(IP)
      IT=INT(TB)
      IT=MIN(IT,19)
      FR=TB-DBLE(IT)
      FR1=ONE-FR
        VISURF=cArp*(FR1*VISCPP(IT+1,1)+FR*VISCPP(IT+2,1))
     1         +cHp*(FR1*VISCPP(IT+1,2)+FR*VISCPP(IT+2,2))
     2       +cAirp*(FR1*VISCPP(IT+1,3)+FR*VISCPP(IT+2,3))
        CPSURF=yArp*(FR1*CSUBPP(IT+1,1)+FR*CSUBPP(IT+2,1))
     1         +yHp*(FR1*CSUBPP(IT+1,2)+FR*CSUBPP(IT+2,2))
     2       +yAirp*(FR1*CSUBPP(IT+1,3)+FR*CSUBPP(IT+2,3))
      ROSURF=AMWA*P(IC)/(RGAS*TP(IP))
      GAMSUR=CPSURF/(CPSURF-RGAS/AMWA)
C==============================================================USERPP
C     CALCULATE REYNOLDS NUMBER AND NUSSELT NUMBER BASED ON
C          PROPERTIES AT THE FILM TEMPERATURE.
C--------------------------------------------------------------USERPP
      REYNUM=ROFILM*RELVEL(IP)*TWO*RADP(IP)/VISFLM
C---------- PRANUM = 0.6*(PRANDTL NUMBER)**(1/3)
      PRANUM=PNTSIX*(VISFLM*CPFILM/CONFLM)**ONEO3
C========================================= KNUDSEN EFFECT =====USERPP
C     CORRECTION FACTOR --- LTE CORRELATION FROM
C          X. CHEN AND E. PFENDER, PLASMA CHEM. PLASMA PROCESS.
C               VOL. 3, 97 (1983); VOL. 3, 351 (1983).
C--------------------------------------------------------------USERPP
C     MEAN FREE PATH LENGTH IS CALCULATED USING HARD SPHERE MODEL AND
C          VISCOSITY AT THE SURFACE.
C--------------------------------------------------------------USERPP
      CORRNC=ONE/(ONE+CORRF*((NINE*GAMSUR-FIVE)/(ONE+GAMSUR))
     1               *SQRT(TP(IP)/AMWA)*VISURF/(P(IC)*RADP(IP)))
C==================== DRAG AND HEAT TRANSFER COEFFICIENTS =====USERPP
C     LTE CORRELATION FROM
C        E. PFENDER AND Y. C. LEE, PLASMA CHEM. PLASMA PROCESS.,
C           5, 211 (1985).
C        Y. C. LEE, Y. P. CHYOU, AND E. PFENDER, IBID, 5, 391 (1985).
C--------------------------------------------------------------USERPP
C     MASS TRANSFER EFFECTS HAS BEEN ADDED BY Y.P. WAN
C     TO THE HTCO
C--------------------------------------- PARTICLE HEATING -----USERPP
      VISKRA=RON(IC)*VISFRE/(ROSURF*VISURF)
      ANUSSL=(TWO+PRANUM*SQRT(REYNUM))*(VISKRA**0.6D0)
     1       *((CPFREE/CPSURF)**0.38D0)
      HTCO(IP)=CORRNC*CONFLM*ANUSSL*HALF/RADP(IP)
C------------------------------ PARTICLE DRAG COEFFICIENT -----USERPP
      CDF=CDC1/REYNUM+SIX/(ONE+SQRT(REYNUM))+CDC2
      CD(IP)=CORRNC**CDC3*CDF/(VISKRA**0.45D0)
C=========== THIS PART IS FOR PARTICLE EVAPORATION BY Y.P. WAN
C                                                     12-16-97
C     NOW JUST FOR ONE MATERIALS
C    MASS DIFFUSIVITY IS BASED ON THE SCHMIDT NUMBER WHICH IS ASSUMED 0.7
      DMPTC(IP)=0.D0
      IBOIL(IP)=0
      IF(IEVAP.EQ.1) THEN
        RHODIFF=VISFLM/SCD
        PVAP=PVAPOR(TP(IP),MTYPE(IP))
        WMETAL=WPOWDER(MTYPE(IP))
        PGDPV=P(IC)*0.1D0/PVAP
        IF(PGDPV.GT.1.D0) THEN
          YMS=WMETAL/(WMETAL+AMWA*(PGDPV-1.D0))
        ELSE
c          WRITE(*,*)'WARNING: PVAPOR GT PGAS',PGDPV
          YMS=1.D0
        END IF
        YMIC=0.D0
C    IF NOT CONSIDER THE KNUDSEN EFFECT
c        SH=2.0D0+0.6D0*SQRT(REYNUM)*SCD**(1.0/3.0)
c    IF CONSIDERING KNUDSEN EFFECT ON MASS TRANSFER 
C       COMPARABLE TO HEAT TRANSFER (CORRECT? NOT CLEAR)
        SH=(2.0D0+0.6D0*SQRT(REYNUM)*SCD**(1.0/3.0))*CORRNC
C
        IF(YMS.LT.1.0D0) THEN
          BNUMB=(YMS-YMIC)/(1.D0-YMS)
          LN1PB=DLOG(1.D0+BNUMB)
        ELSE
        END IF
        IF(IKIVA.EQ.1) THEN
C    USE THE VAPORIZATION MODEL IN KIVA
C    MASS CONCENTRATION OF METAL IS ASSUMED TO BE ZERO
C    NEW RADIUS OF THE PARTICLE IP
          IF(YMS.LT.1.D0) THEN
            DMDTTMP=-PI*2.D0*RADP(IP)*RHODIFF*SH*LN1PB
          ELSE
            DMDTTMP=-1.D30
          END IF
        ELSE IF(IKIVA.EQ.0) THEN
C    USE THE VAPORIZATION MODEL IN WESTHOFF'S PAPER, CONSIDERING THE
C        SO CALLED LANGMIUR EVAPORATION RATE
          HE=SQRT(RGAS*TP(IP)/(PI2*WMETAL))
          HM=RHODIFF*SH/(RADP(IP)*2.D0*ROFILM)
          HMIX=1.D0/(1.D0/HE+1.D0/HM)
          DMDTTMP=-ROFILM*PI4*RADP(IP)*RADP(IP)*HMIX*(YMS-YMIC)
        ELSE
          WRITE(*,*) 'ERROR IN IKIVA SELECTION'
          STOP
        END IF
        RADPTMP=DMDTTMP*DT(IC1)/(PI4O3*ROPTC(MTYPE(IP))*1.D-3)
     &          +RADP(IP)**3
        RADPTMP=RADPTMP**(1./3.)
c ---- NOTE: THE LAST COMMAND WILL INTRODUCE ERROR IF DMDTTMP VERY SMALL
C            WHICH CAUSES RADPTMP.GT.RADP(IP) ALTHOUGH DMDTTMP<0
        IF(DMDTTMP.GT.0.D0) THEN
          WRITE(*,*)'ERROR IN USERPP, DMEVAP>0',DMDTTMP
          STOP
        END IF
        IF(RADPTMP.GT.RADP(IP).AND.DMDTTMP.LE.0.D0) THEN
          RADPTMP=RADP(IP)
        END IF
C 
        DMPTC(IP)=DMDTTMP
      END IF
C==============================================================USERPP
C ---- CHECK IF THE EVAPORATION RATE TOO LARGE OR IF THE HEAT TRANSFER
C      RATE IS THE CONTROLING FACTOR
      IF(IEVAP.EQ.1) THEN
        QCON=HTCO(IP)*(TEMP(IC)-TP(IP))
        QRAD=STEBOL*EMSSP(IP)*(TEMAMB**4-TP(IP)**4)
c        DMHEAT=-(QCON+QRAD)*PI4*RADP(IP)*RADP(IP)/HLTNT(MTYPE(IP))
        DMHEAT=-(QCON+QRAD)*PI4*RADP(IP)*RADP(IP)/
     &         (LATVAP(TP(IP),MTYPE(IP))*1.D4)
C       IF USING CONVENTIONAL HEAT TRANSFER CONTROLLED EVAPORIZATION
c        IF(TP(IP).LE.5100.d0) THEN
c          DMPTC(IP)=0.D0
c          RADPTMP=RADP(IP)
c            write(*,*)'*** time=',time,tp(ip)
c        ELSE
c          DMPTC(IP)=-1.D30
c          RADPTMP=RADP(IP)
c        END IF
C       MIXING THE CONTROLLING FACTORS (DIFFUSION, HEAT TRANSFER)
        IF(DMHEAT.GT.DMPTC(IP)) THEN
          IF(DMHEAT.LT.0.D0) THEN
c           write(*,*)'+++++ heat controled'
C   --- EVAPORATION RATE CONTROLED BY HEAT TRANSFER
            RADPTMP=DMHEAT*DT(IC1)/(PI4O3*ROPTC(MTYPE(IP))*1.D-3)
     &            +RADP(IP)**3
            RADPTMP=RADPTMP**(1./3.)
c           write(*,*)'+++ dmheat,tp=',dmheat,tp(ip),radptmp
            IBOIL(IP)=1
          ELSE
            DMHEAT=0.D0
            RADPTMP=RADP(IP)
          END IF
          DMPTC(IP)=DMHEAT
        END IF
C    SEE IF RADIUS <0, THIS PARTICLE WILL BE REMOVED IN SUBROUTINE REPACK
C    SET SMALL VALUE TO THIS RADIUS TO AVOID DIVISION BY ZERO
        IF(RADPTMP.LE.0.D0) THEN
          RADPTMP=1.D-8
          WRITE(*,*)'*** PARTICLE IP=',IP,' WILL BE 
     & REMOVED AT TIME',TIME
        END IF
C   --- SET THE PARTICLE RADIUS AND THE CHANGE OF THE RADIUS
        DRD(IP)=(RADPTMP-RADP(IP))*1.D-2
        IF(ICONDP.EQ.1) THEN
          RD(IP)=RADPTMP*1.D-2
        END IF
        RADP(IP)=RADPTMP
        PMASS(IP)=PI4O3*ROPTC(MTYPE(IP))*1.D-3*RADP(IP)**3
C ---- IN CASE RM,RRS > RD
        IF(RM(IP).GT.RD(IP)) THEN
c          WRITE(*,*)'*** WARNING, RM>RD',RM(IP),RD(IP)
          RM(IP)=RD(IP)
        END IF
        IF(RRS(IP).GT.RD(IP)) THEN
c          WRITE(*,*)'*** WARNING, RRS>RD',RRS(IP),RD(IP)
          RRS(IP)=RD(IP)
        END IF
      END IF
C ---- NUSSL NUMBER SHOULD BE MODIFIED BY MASS TRANSFER NUMBER 
C                 BY Y.P. WAN  12-16-97 IS THIS CORRECT?
      FACB=1.D0
      IF(IEVAP.EQ.1.AND.DMPTC(IP).LT.0.D0) THEN
        FACA=-DMPTC(IP)*CPFILM/(2.D0*PI*RADP(IP)*CONFLM)
C  --  TO AVOID NUMERICAL ERROR CAUSED BY SMALL FACB
        IF(FACA.GT.1.D-8) THEN
          FACB=FACA/(EXP(FACA)-1.D0)
          HTCO(IP)=HTCO(IP)*FACB
        END IF
      END IF
C -- For single particle, print out particle information
C    not recommanded for many particles, by Y.P. WAN
c     CALL PRINTPTC(IP)
c     CALL PRINTAMB(TIME,TEMP(IC),TP(IP),REYNUM,ANUSSL,SH,PRANUM,SCD,
c    &     ROFILM,VISFLM,CPFILM,CONFLM,RHODIFF,HTCO(IP),CD(IP),
c    &     CORRNC,RELVEL(IP),QCON,QRAD,FACB)
  200 CONTINUE
      RETURN
      END
