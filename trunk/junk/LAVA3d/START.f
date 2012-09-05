*DECK START
      SUBROUTINE START
C===============================================================START
C     THIS ROUTINE PROVIDES BOUNDARY CONDITIONS AND INITIAL VALUES
C     INITIALIZE EPS ONLY IF ITURB=2
C----------
C     CALLED BY LAVA
C===============================================================START
      INCLUDE 'COML.h'
C---------------------------------------------------------------START
C     TKE MUST BE SET TO BE ZERO FOR LAMINAR CASES
C     FOR BOTH INITIAL AND BOUNDARY CONDITIONS
C---------------------------------------------------------------START
C     IN CASE OF NLTE>=1 AND ELECTRON DENSITY IS NOT ZERO, SPECIAL
C     CODING IS REQUIRED TO DO THIS HERE AND BCCC. THIS SUBROUTINE IS
C     NOT SET UP TO DO SUCH CASES.
C---------------------------------------------------------------START
      INTEGER ITL,ITR,ITD,ITF,ITB,ITT
      INTEGER INOZ
      DOUBLE PRECISION XMFI(NSP),XMFL(NSP),XMFR(NSP),
     1     XMFD(NSP),XMFF(NSP),XMFB(NSP),XMFT(NSP),
     2     TB,FR,ROCV,FRL,FRR,FRD,FRF,FRB,FRT
C---------------------------------------------------------------START
      EXTERNAL CHEMEV
      INTEGER LNREV
      PARAMETER (LNREV=3)
      DOUBLE PRECISION EQC(LNREV),PROG(LNREV),BAM(LNREV,LNREV),
     1     RHS(LNREV),SPDSV(NSP),SPDNEG(NSP),
     2 DPROG(LNREV),DSPD(LNREV),WORK(5*LNREV),SWORK(LNREV),
     3  RWORK(LNREV),CWORK(LNREV),BAMF(LNREV,LNREV)
      DOUBLE PRECISION ASV(LNREV),BSV(LNREV),CSV(LNREV),DSV(LNREV),
     1     ESV(LNREV),FBNANV(NSP,LNREV)
      INTEGER NLMV(LNREV),ANV(NSP,LNREV),BNV(NSP,LNREV),CNV(NSP,LNREV)
      INTEGER IWORK(LNREV),ISKIP(LNREV),IPCON(LNREV),IPIV(LNREV),
     1     IREV,NK
C---------------------------------------------------------------START
C
      DOUBLE PRECISION SDN(NSP)
      DOUBLE PRECISION UDEL,DVDX,DVDXMX,TKMAX,DELTA
      INTEGER IDEL
C---------------------------------------------------------------START
      double precision vdotAr,vdotH2,xiH2,amssAr,amssH2
      double precision penv,poweri,amssin,amass,
     1     Hdotin,Hdtout,swirnn,wpeak,ucl,tcl,ctot,
     2     rhof,sumhf,rocvi,areav,spdsnd,power,vratio,tratio,
     3     wpxp,rhofr,rhofc,priter,axial,radial,rordr,sratio
      double precision tl,pold,pnew,pcur,res
      DOUBLE PRECISION VCRIT
      integer iter,ipeak,ngo
c--- icin for nozzle radius, icout for torch outside of METCO-9MB
c      data icin,icout /11,30/
      DATA VCRIT /0.01D0/
C===============================================================START
      WRITE(NFOUT,1000)
 1000 FORMAT(5X,'-----------------------------------------------',/,
     1       5X,'START FROM BEGINNING')
      IPBC=0
C============== SET DERRIERE AND FRONT BOUNDARY CONDITIONS =====START
      IF(NY.GT.1) THEN
C---------------------------------------------------------------START
      DO 110 ICDF=1,NDF
c     PBCD(ICDF)=ZERO
c-----
      if(icdf.le.icout) then
        pbcd(icdf)=zero
      else
        pbcd(icdf)=one
      endif
c-----
      UVELD(ICDF)=ZERO
      VVELD(ICDF)=ZERO
      WVELD(ICDF)=ZERO
      PRESD(ICDF)=PAMB
c     TEMPD(ICDF)=TEMAMB
c-----  torch tip wall. but Why?
      if(icdf.gt.icin .and. icdf.le.icout) then
        tempd(icdf)=700.0d0+(temamb-700.0d0)
     1              *log(xc(icdf)/x(icin))/log(x(icout)/x(icin))
      else
        tempd(icdf)=temamb
      endif
c-----
      EDEND(ICDF)=ZERO
c     TEED(ICDF)=TEMAMB
      teed(icdf)=tempd(icdf)
      EDEED(ICDF)=ZERO
      TKED(ICDF)=SMALL
      IF(ITURB.EQ.0) TKED(ICDF)=ZERO
      IF(ITURB.EQ.2) EPSD(ICDF)=SMALL
c     PBCF(ICDF)=ZERO
      pbcf(icdf)=one
      UVELF(ICDF)=ZERO
      VVELF(ICDF)=ZERO
      WVELF(ICDF)=ZERO
      PRESF(ICDF)=PAMB
      TEMPF(ICDF)=TEMAMB
      EDENF(ICDF)=ZERO
      TEEF(ICDF)=TEMAMB
      EDEEF(ICDF)=ZERO
      TKEF(ICDF)=SMALL
      IF(ITURB.EQ.0) TKEF(ICDF)=ZERO
      IF(ITURB.EQ.2) EPSF(ICDF)=SMALL
      IF(PBCD(ICDF).EQ.ONE .OR. PBCF(ICDF).EQ.ONE) IPBC=1
  110 CONTINUE
C----------
      DO 140 ICDF=1,NDF
      ITD=INT(HUNDTH*TEMPD(ICDF))
      FRD=HUNDTH*TEMPD(ICDF)-DBLE(ITD)
      ITF=INT(HUNDTH*TEMPF(ICDF))
      FRF=HUNDTH*TEMPF(ICDF)-DBLE(ITF)
c----- air environment, 
c ----- should be in input mass fraction in air (O2+N2).
       xmfd(1)=0.00d0
       xmfd(2)=0.00d0
       xmfd(3)=0.00d0
       xmfd(4)=0.00d0
       xmfd(5)=0.00d0
       xmfd(6)=0.79d0
       xmfd(7)=0.00d0
       xmfd(8)=0.00d0
       xmfd(9)=0.00d0
       xmfd(10)=0.21d0
       xmfd(11)=0.00d0
       xmfd(12)=0.00d0
       xmfd(13)=0.00d0
       xmfd(14)=0.00d0
       xmfd(15)=0.00d0
c----- air environment
       xmff(1)=0.00d0
       xmff(2)=0.00d0
       xmff(3)=0.00d0
       xmff(4)=0.00d0
       xmff(5)=0.00d0
       xmff(6)=0.79d0
       xmff(7)=0.00d0
       xmff(8)=0.00d0
       xmff(9)=0.00d0
       xmff(10)=0.21d0
       xmff(11)=0.00d0
       xmff(12)=0.00d0
       xmff(13)=0.00d0
       xmfd(14)=0.00d0
       xmfd(15)=0.00d0
c-----
        DO 130 ISP=1,NSP
        SDND(ICDF,ISP)=XMFD(ISP)*MW(ISP)*PRESD(ICDF)/(RGAS*TEMPD(ICDF))
        SDNF(ICDF,ISP)=XMFF(ISP)*MW(ISP)*PRESF(ICDF)/(RGAS*TEMPF(ICDF))
c ----  is this internal energy?
        EDEND(ICDF)=EDEND(ICDF)+SDND(ICDF,ISP)
     1              *((ONE-FRD)*EK(ITD+1,ISP)+FRD*EK(ITD+2,ISP))
        EDENF(ICDF)=EDENF(ICDF)+SDNF(ICDF,ISP)
     1              *((ONE-FRF)*EK(ITF+1,ISP)+FRF*EK(ITF+2,ISP))
  130   CONTINUE
c ----- for non-thermal-equil.
      IF(NLTE.NE.0) THEN
        EDEED(ICDF)=SDND(ICDF,IELC)
     1    *((ONE-FRD)*EK(ITD+1,IELC)+FRD*EK(ITD+2,IELC))
        EDEEF(ICDF)=SDNF(ICDF,IELC)
     1    *((ONE-FRF)*EK(ITF+1,IELC)+FRF*EK(ITF+2,IELC))
      END IF
  140 CONTINUE
      END IF
C---------------------------------------- INPUT FOR CHEMEV -----START
C     DETERMINE ASV -- ESV, FBNANV,ANV,BNV,CNV,NLMV FOR CHEMEV
C---------------------------------------------------------------START
      DO 870 IREV=1,LNREV
C----------
      IF(irev.le.1) THEN
        ASV(IREV)=AKS(IREV)
        BSV(IREV)=BKS(IREV)
        CSV(IREV)=CKS(IREV)
        DSV(IREV)=DKS(IREV)
        ESV(IREV)=EKS(IREV)
      ELSE
        ASV(IREV)=AS(IREV-1)
        BSV(IREV)=BS(IREV-1)
        CSV(IREV)=CS(IREV-1)
        DSV(IREV)=DS(IREV-1)
        ESV(IREV)=ES(IREV-1)
      END IF
C----------
      NK=0
      DO 860 ISP=1,NSP
C----------
      IF(irev.le.1) THEN
        ANV(ISP,IREV)=AM(ISP,IREV)
        BNV(ISP,IREV)=BM(ISP,IREV)
      ELSE
        ANV(ISP,IREV)=AN(ISP,IREV-1)
        BNV(ISP,IREV)=BN(ISP,IREV-1)
         if(isp.eq.13) then
            ANV(ISP,IREV)=ANV(ISP,IREV)-1
            BNV(ISP,IREV)=BNV(ISP,IREV)-1
         endif
      END IF
C----------
      FBNANV(ISP,IREV)=DBLE(BNV(ISP,IREV)-ANV(ISP,IREV))
      IF(ANV(ISP,IREV).EQ.0 .AND. BNV(ISP,IREV).EQ.0) GO TO 860
      NK=NK+1
      CNV(NK,IREV)=ISP
  860 CONTINUE
      NLMV(IREV)=NK
  870 CONTINUE
C==================== conditions at the exit of the nozzle =====START
c--- THIS PART HAS BEEN CHANGED TO THE INPUT FILE
C     torch power = VOLTG * Amp * Efficiency * 1.0e7 (Joule to erg)
C---------------------------------------------------------------START
c*** if different pressure is desired, change in pamb in HOT
C      vdotAr=40.0d0
C      vdotH2=12.0d0
c reduce the flow rate see the convergence of the code
      xiH2=FLRT2/(FLRT1+FLRT2)
C       amp=600.0d0
C      VOLTG=40.0d0
C      eff=0.70d0
      penv=pamb/P1ATM
C---------------------------------------------------------------START
      PRINT*,        'TORCH OPERATING CONDITIONS ARE GIVEN AS'
      PRINT*,        CURCY,' Amp, ',VOLTG,' VOLTG, Eff = ',EFFCY
      PRINT*,        'Ar =',FLRT1,' standard l/min, ',penv,' ATM'
      PRINT*,        'H2 =',FLRT2,' standard l/min'
      PRINT*,        'SWIRL NUMBER = ',swirn
C--------------- mass flow rate g/s ----------------------------START
      poweri=1.0d7*VOLTG*CURCY*EFFCY
      vdotAr=FLRT1*1.0d3/60.0d0
      vdotH2=FLRT2*1.0d3/60.0d0
      amssAr=mw(1)*P1ATM/(Rgas*300.0d0)*VdotAr
      amssH2=mw(3)*P1ATM/(Rgas*300.0d0)*VdotH2
      amssin=amssAr+amssH2
C---------------------------------------------------------------START
       tb=hundth*temamb
       it=int(tb)
       fr=tb-dble(it)
      Hdotin=amssAr*((htform(1)+rgas*temamb)*rmw(1)
     1              +(one-fr)*ek(it+1,1)+fr*ek(it+2,1))
     2      +amssH2*((htform(3)+rgas*temamb)*rmw(3)
     3              +(one-fr)*ek(it+1,3)+fr*ek(it+2,3))
C---------------------------------------------------------------START
C      ipeak assumed to be Rin*0.8
C      SWIRL NUMBER =  9.999E-02 peak w= 44198.44959906817
C---------------------------------------------------------------START
      do 820 i=2,icin
      if(xc(i).gt.0.8d0*x(icin)) then
        ipeak=i-1
        go to 822
      endif
  820 continue
  822 continue
C------------------------------------------ ITERATION LOOP -----START
  900 continue
c      PRINT*,'Give powerv, powert, swirl number'
c      READ *,powerv,powert,swirn
c  ---- no read, through standard input by Y.P. WAN 7/2/97
C      powerv=20.0
C      powert=20.0
C      swirn=0.0
C---------------------------------------------------------------START
C     IF SWIRL PROFILE IS GIVEN, THEN SET SWIRN=0, AND GIVE SWIRL
C     PROFILE HERE
C---------------------------------------------------------------START
      if(swirn.eq.zero) swirnn=zero
      if(iswirl.eq.0 .and. swirn.ne.zero) call exita(1)
      if(iswirl.ne.0 .and. swirn.eq.zero) call exita(1)
      if(iswirl.eq.0) wpeak=zero
      if(iswirl.eq.1) wpeak=1000.0d0
C----------- why assume such a data ?   ------------------------START
      Ucl=178318.5307618468d0
      Tcl=13510.41156817594d0
  930 continue
C---------------------------------------------------------------START
      amass=zero
      Hdtout=zero
      do 750 i=2,icin
        tempd(i)=(tcl-700.0d0)*(one-(xc(i)/x(icin))**PWTEP)+700.0d0
        ctot=presd(i)/(Rgas*tempd(i))
C---------------------------------------------------------------START
C     get species densities using Newton's method and CHEMEV
C     initial guess needs to be consistent with all constraints.
C     (in array sdn).
C---------------------------------------------------------------START
c ---- this change of initial guessed values is done by Y.P. WAN 7/7/97
      if(i.eq.2) then
        sdn(1)= mw(1)*(one-xiH2)*ctot
        sdn(2)= 0.0d0
        sdn(3)= mw(3)*xiH2*ctot
        sdn(4)= 0.0d0
        sdn(5)= 0.0d0
        sdn(6)= 0.0d0
        sdn(7)= 0.0d0
        sdn(8)= 0.0d0
        sdn(9)= 0.0d0
        sdn(10)=0.0d0
        sdn(11)=0.0d0
        sdn(12)=0.0d0
        sdn(13)=0.0d0
        sdn(14)=0.0d0
        sdn(15)=0.0d0
      end if
C---------------------------------------------------------------START
      tl=tempd(i)
      pold=presd(i)
      iter=0
  700 continue
      iter=iter+1
      if(tl.gt.1250.0d0) then
      CALL CHEMEV(tl,sdn,NSP,LNREV,LNREV,
     1        ASV,BSV,CSV,DSV,ESV,ANV,BNV,CNV,NLMV,FBNANV,MW,RMW)
      endif
c      print*,i,iter,sdn(1),sdn(2),sdn(3),sdn(4),sdn(5)
c      print*,i,iter,tl
      ctot=zero
      do 720 isp=1,nsp
        ctot=ctot+sdn(isp)*rmw(isp)
  720 continue
      pcur=ctot*Rgas*tl
      res=pcur-presd(i)
      if(abs(res/presd(i)).lt.1.0d-6) go to 735
C----------------------------------------------- slope = 1 -----START
      pnew=pcur-res
      pnew=max(pnew,0.1d0*pcur)
      pnew=min(pnew,10.0d0*pcur)
C---------------------------------------------------------------START
        do 730 isp=1,nsp
          sdn(isp)=sdn(isp)*pnew/pcur
  730   continue
        pold=pcur
        pcur=pnew
      if(iter.gt.50) then
        write(*,*)'not converged, iter=',iter,', tl=',tl
        stop
      endif
        go to 700
C---------------------------------------------------------------START
  735 continue
C---------------------------------------------------------------START
        tb=hundth*tempd(i)
        it=int(tb)
        fr=tb-dble(it)
        edend(i)=zero
        rhof=zero
        sumhf=zero
        rocvi=zero
        do 740 isp=1,nsp
        sdnd(i,isp)=sdn(isp)
        edend(i)=edend(i)+sdnd(i,isp)*((1.0d0-fr)*ek(it+1,isp)
     1                                +fr*ek(it+2,isp))
        rhof=rhof+sdnd(i,isp)
        sumhf=sumhf+sdn(isp)*rmw(isp)*htform(isp)
        rocvi=rocvi+sdn(isp)*(ek(it+2,isp)-ek(it+1,isp))*hundth
  740   continue
        vveld(i)=ucl*(one-(xc(i)/x(icin))**PWVEL)
        areav=two*pi*xc(i)*dx(i)*vveld(i)
        amass=amass+rhof*areav
        Hdtout=Hdtout+areav*(edend(i)+presd(i)+sumhf)
c----- frozen speed of sound
        if(i.eq.2) spdsnd=hundth*sqrt((one+presd(i)/(rocvi*tempd(i)))
     1                                *presd(i)/rhof)
  750 continue
C---------------------------------------------------------------START
      power=Hdtout-Hdotin
c      vratio=amssin/amass
c      tratio=poweri/power
c ---- change the criteria for the convergence by Y.P. WAN 7/7/97
c      if(tratio.le.zero) then
c        write(*,*) 'ERROR!!!! --- negative tratio ....'
c        stop
c      endif
c      if(vratio.ge.0.9999d0 .and. vratio.le.1.0001d0 .and.
c     1   tratio.ge.0.9999d0 .and. tratio.le.1.0001d0) go to 960
c      ucl=vratio*ucl*0.25d0+0.75d0*ucl
c      tcl=tratio*tcl*0.10d0+0.90d0*tcl
      vratio=dabs((amssin-amass)/amssin)
      tratio=dabs((poweri-power)/poweri)
c -- the relative difference as criteria vcrit=1%
      if(vratio.le.vcrit.and.tratio.le.vcrit) go to 960
c      write(*,*)'****** ucl=',ucl,', tcl=',tcl,
c     & ', amssin=',amssin,',poweri=',poweri
c -- modify the values of ucl and tcl 
      ucl=amssin/amass*ucl*0.10d0+0.90d0*ucl
      tcl=poweri/power*tcl*0.10d0+0.90d0*tcl
c      write(*,*) '****** ucl=',ucl,', tcl=',tcl
c     & ',amass=',amass, ',power=',power
      go to 930
  960 continue
C--------- success on the iteration of u and T      -----------START
      vveld(1)=vveld(2)
      tempd(1)=tempd(2)
      edend(1)=edend(2)
      teed(1)=teed(2)
      edeed(1)=edeed(2)
      do 770 isp=1,nsp
      sdnd(1,isp)=sdnd(2,isp)
  770 continue
  
C--------------------- INFLOW PROFILE FOR SGS OR K-E MODEL -----START
C    MIXING LENGTH AT INLET FOR K-E MODEL IS SPECIFIED USING
C    THE FORMULA SUGGESTED BY LESCHZINER AND RODI FOR SWIRL JET
C    (AIAA 22, P1742, 1984).
C---------------------------------------------------------------START
      PRINT*,        'ucl', ucl,' cm/s, ',' tcl', tcl,' K'
      IF(ITURB.NE.0) THEN
        UDEL=TENTH*VVELD(2)
        DVDXMX=SMALL
        DO 785 I=2,ICIN
          IF(VVELD(I).GE.UDEL) IDEL=I
          DVDX=ABS((VVELD(I)-VVELD(I-1))/(XC(I)-XC(I-1)))
          DVDXMX=MAX(DVDXMX,DVDX)
  785   CONTINUE
        DELTA=TWO*(XC(IDEL)+(UDEL-VVELD(IDEL))*
     1        (XC(IDEL+1)-XC(IDEL))/(VVELD(IDEL+1)-VVELD(IDEL)))
        TKMAX=HUNDTH*VVELD(2)*VVELD(2)
C----- reduced turbulence intensity
c       tkmax=0.001d0*VVELD(2)*VVELD(2)
        DO 780 I=3,ICIN
          DVDX=ABS((VVELD(I)-VVELD(I-1))/(XC(I)-XC(I-1)))
          IF(ITURB.EQ.1) THEN
C------ USER MAY WANT TO USE DIFFERENT VALUE FOR SGS MODEL -----START
            TKED(I)=THRHAF*TKMAX*DVDX/DVDXMX
            TKED(I)=MAX(TKED(I),SMALL)
          END IF
          IF(ITURB.EQ.2) THEN
            TKED(I)=THRHAF*TKMAX*DVDX/DVDXMX
            TKED(I)=MAX(TKED(I),SMALL)
            EPSD(I)=TKED(I)**1.5*CMU**0.75/(0.075D0*DELTA)
          END IF
  780   CONTINUE
        TKED(2)=TKED(3)
        EPSD(2)=EPSD(3)
        TKED(1)=TKED(2)
        EPSD(1)=EPSD(2)
      END IF
C--------------------------------- inflow profile of swirl -----START
c     This part needs to be provided by the user
C---------------------------------------------------------------START
      if(iswirl.eq.0) go to 990
  800 continue
      do 825 i=2,ipeak
      wveld(i)=xc(i)*wpeak/xc(ipeak)
  825 continue
      wpxp=wpeak/(x(icin)-xc(ipeak))
      do 827 i=ipeak+1,icin
      wveld(i)=wpeak-wpxp*(xc(i)-xc(ipeak))
  827 continue
      wveld(1)=-wveld(2)
C---------------------------------------------------------------START
c     pressure at the inside wall of nozzle should be gussed
C---------------------------------------------------------------START
      rhofr=pamb*mw(1)/(Rgas*700.0d0)
      rhofc=zero
       do 805 isp=1,nsp
       rhofc=rhofc+sdnd(icin,isp)
  805  continue
      presd(icin)=pamb-wveld(icin)**2*dx(icin)
     1                 *(rhofc+rhofr)/(eight*x(icin))
C---------------------------------------------------------------START
      do 810 i=icin-1,2,-1
        rhofr=zero
        rhofc=zero
         do 807 isp=1,nsp
         rhofr=rhofr+sdnd(i+1,isp)
         rhofc=rhofc+sdnd(i,isp)
  807    continue
        presd(i)=presd(i+1)+twoo3*(rhofr*tked(i+1)-rhofc*tked(i))
     1           -dx(i)*(rhofr+rhofc)*(wveld(i+1)+wveld(i))**2
     2            /(eight*x(i))
  810 continue
c----- presd(1) contains pressure of previous iteration
      priter=presd(1)/presd(2)
      if(priter.gt.0.9999d0 .and. priter.le.1.0001d0) go to 830
        do 840 i=2,icin
        presd(i)=presd(i)*(0.6d0+0.4d0*priter)
  840   continue
        presd(1)=presd(2)
        go to 930
  830 continue
      presd(1)=presd(2)
C---------------------------------------------------------------START
      axial=zero
      radial=zero
      do 850 i=2,icin
        rhof=zero
        do 845 isp=1,nsp
        rhof=rhof+sdnd(i,isp)
  845   continue
        rordr=rhof*xc(i)*dx(i)
        axial=axial+rordr*(vveld(i)**2-half*wveld(i)**2)
        radial=radial+rordr*vveld(i)*wveld(i)*xc(i)
  850 continue
      swirnn=radial/(axial*x(icin))
      if(swirn.eq.zero) go to 990
        sratio=swirn/swirnn
        if(sratio.gt.0.9999d0 .and. sratio.le.1.0001d0) go to 990
       print*,'---',swirnn,sratio,wpeak
        wpeak=wpeak*sratio
       print*,'---',swirnn,sratio,wpeak
        go to 800
  990 continue
C---------------------------------------------------------------START
c----- Swank's experiment - swirl is clockwise viewed from torch
c     do 995 i=1,icin
c     wveld(i)=-wveld(i)
c 995 continue
C---------------------------------------------------------------START
      power=power*1.d-10
      PRINT*,        'MASS= ',amass,'g/s, POWER= ',power,'KW'
      PRINT*,        'Ucl=',ucl,' cm/s, powerv',PWVEL
      PRINT*,        'Tcl=',tcl,' K, powert',PWTEP
      PRINT*,        'Frozen Speed of Sound= ',spdsnd,' m/sec'
      PRINT*,        'SWIRL NUMBER = ',swirnn,' peak w=',wpeak
      PRINT*,        'Pcl = ',presd(2),' ,Pamb = ',penv,' ATM'
c      PRINT*,'Satisfied? If yes, give 0; give 1 otherwise.'
c      READ *,ngo
c  ---- by Y.P. WAN 7/2/97
      ngo=0
      if(ngo.ne.0) go to 900
C---------------------------------------------------------------START
      WRITE(NFOUT,*) '---------------------------------------'
      WRITE(NFOUT,*) 'TORCH OPERATING CONDITIONS ARE GIVEN AS'
      WRITE(NFOUT,*) CURCY,' Amp, ',VOLTG,' VOLTG, Eff = ',EFFCY
      WRITE(NFOUT,*) 'Ar =',vdotAr,' standard cm**3/sec'
      WRITE(NFOUT,*) 'H2 =',vdotH2,' standard cm**3/sec'
      WRITE(NFOUT,*) 'MASS= ',amass,'g/s, POWER= ',power,'KW'
      WRITE(NFOUT,*) 'Ucl=',ucl,' cm/s, powerv',PWVEL
      WRITE(NFOUT,*) 'Tcl=',tcl,' K, powert',PWTEP
      WRITE(NFOUT,*) 'Frozen Speed of Sound= ',spdsnd,' m/sec'
      WRITE(NFOUT,*) 'SWIRL NUMBER = ',swirn,' peak w=',wpeak
      WRITE(NFOUT,*) 'Pcl = ',presd(2),' ,Pamb = ',penv,' ATM'
      WRITE(NFOUT,*) '---------------------------------------'
C================== SET LEFT AND RIGHT BOUNDARY CONDITIONS =====START
      IF(NX.GT.1) THEN
      DO 10 ICLR=1,NLR
      PBCL(ICLR)=ZERO
      UVELL(ICLR)=ZERO
      VVELL(ICLR)=ZERO
      WVELL(ICLR)=ZERO
      PRESL(ICLR)=PAMB
      TEMPL(ICLR)=TEMAMB
      EDENL(ICLR)=ZERO
      TEEL(ICLR)=TEMAMB
      EDEEL(ICLR)=ZERO
      TKEL(ICLR)=SMALL
      IF(ITURB.EQ.0) TKEL(ICLR)=ZERO
      IF(ITURB.EQ.2) EPSL(ICLR)=SMALL
c     PBCR(ICLR)=ZERO
      pbcr(iclr)=one
      UVELR(ICLR)=ZERO
      VVELR(ICLR)=ZERO
      WVELR(ICLR)=ZERO
      PRESR(ICLR)=PAMB
      TEMPR(ICLR)=TEMAMB
      EDENR(ICLR)=ZERO
      TEER(ICLR)=TEMAMB
      EDEER(ICLR)=ZERO
      TKER(ICLR)=SMALL
      IF(ITURB.EQ.0) TKER(ICLR)=ZERO
      IF(ITURB.EQ.2) EPSR(ICLR)=SMALL
      IF(PBCR(ICLR).EQ.ONE .OR. PBCL(ICLR).EQ.ONE) IPBC=1
   10 CONTINUE
C---------------------------------------------------------------START
      DO 40 ICLR=1,NLR
      ITL=INT(HUNDTH*TEMPL(ICLR))
      FRL=HUNDTH*TEMPL(ICLR)-DBLE(ITL)
      ITR=INT(HUNDTH*TEMPR(ICLR))
      FRR=HUNDTH*TEMPR(ICLR)-DBLE(ITR)
c----- air environment
       xmfl(1)=0.00d0
       xmfl(2)=0.00d0
       xmfl(3)=0.00d0
       xmfl(4)=0.00d0
       xmfl(5)=0.00d0
       xmfl(6)=0.79d0
       xmfl(7)=0.00d0
       xmfl(8)=0.00d0
       xmfl(9)=0.00d0
       xmfl(10)=0.21d0
       xmfl(11)=0.00d0
       xmfl(12)=0.00d0
       xmfl(13)=0.00d0
       xmfl(14)=0.00d0
       xmfl(15)=0.00d0
c----- air environment
         xmfr(1)=0.00d0
         xmfr(2)=0.00d0
         xmfr(3)=0.00d0
         xmfr(4)=0.00d0
         xmfr(5)=0.00d0
         xmfr(6)=0.79d0
         xmfr(7)=0.00d0
         xmfr(8)=0.00d0
         xmfr(9)=0.00d0
         xmfr(10)=0.21d0
         xmfr(11)=0.00d0
         xmfr(12)=0.00d0
         xmfr(13)=0.00d0
         xmfr(14)=0.00d0
         xmfr(15)=0.00d0
c----- boundary values of rho and rhoe
        DO 30 ISP=1,NSP
        SDNL(ICLR,ISP)=XMFL(ISP)*MW(ISP)*PRESL(ICLR)/(RGAS*TEMPL(ICLR))
        SDNR(ICLR,ISP)=XMFR(ISP)*MW(ISP)*PRESR(ICLR)/(RGAS*TEMPR(ICLR))
        EDENL(ICLR)=EDENL(ICLR)+SDNL(ICLR,ISP)
     1              *((ONE-FRL)*EK(ITL+1,ISP)+FRL*EK(ITL+2,ISP))
        EDENR(ICLR)=EDENR(ICLR)+SDNR(ICLR,ISP)
     1              *((ONE-FRR)*EK(ITR+1,ISP)+FRR*EK(ITR+2,ISP))
   30   CONTINUE
      IF(NLTE.NE.0) THEN
        EDEEL(ICLR)=SDNL(ICLR,IELC)
     1    *((ONE-FRL)*EK(ITL+1,IELC)+FRL*EK(ITL+2,IELC))
        EDEER(ICLR)=SDNR(ICLR,IELC)
     1    *((ONE-FRR)*EK(ITR+1,IELC)+FRR*EK(ITR+2,IELC))
      END IF
   40 CONTINUE
      END IF
C================== SET TOP AND BOTTOM BOUNDARY CONDITIONS =====START
      IF(NZ.GT.1) THEN
      DO 210 ICBT=1,NBT
      PBCB(ICBT)=ZERO
      UVELB(ICBT)=ZERO
      VVELB(ICBT)=ZERO
      WVELB(ICBT)=ZERO
      PRESB(ICBT)=PAMB
      TEMPB(ICBT)=TEMAMB
      EDENB(ICBT)=ZERO
      TEEB(ICBT)=TEMAMB
      EDEEB(ICBT)=ZERO
      TKEB(ICBT)=SMALL
      IF(ITURB.EQ.0) TKEB(ICBT)=ZERO
      IF(ITURB.EQ.2) EPSB(ICBT)=SMALL
      PBCT(ICBT)=ZERO
      UVELT(ICBT)=ZERO
      VVELT(ICBT)=ZERO
      WVELT(ICBT)=ZERO
      PREST(ICBT)=PAMB
      TEMPT(ICBT)=TEMAMB
      EDENT(ICBT)=ZERO
      TEET(ICBT)=TEMAMB
      EDEET(ICBT)=ZERO
      TKET(ICBT)=SMALL
      IF(ITURB.EQ.0) TKET(ICBT)=ZERO
      IF(ITURB.EQ.2) EPST(ICBT)=SMALL
      IF(PBCB(ICBT).EQ.ONE .OR. PBCT(ICBT).EQ.ONE) IPBC=1
  210 CONTINUE
C----------
      DO 240 ICBT=1,NBT
      ITB=INT(HUNDTH*TEMPB(ICBT))
      FRB=HUNDTH*TEMPB(ICBT)-DBLE(ITB)
      ITT=INT(HUNDTH*TEMPT(ICBT))
      FRT=HUNDTH*TEMPT(ICBT)-DBLE(ITT)
c-----
         xmfb(1)=1.00d0
c-----
         xmft(1)=1.00d0
c-----
        DO 230 ISP=1,NSP
        SDNB(ICBT,ISP)=XMFB(ISP)*MW(ISP)*PRESB(ICBT)/(RGAS*TEMPB(ICBT))
        SDNT(ICBT,ISP)=XMFT(ISP)*MW(ISP)*PREST(ICBT)/(RGAS*TEMPT(ICBT))
        EDENB(ICBT)=EDENB(ICBT)+SDNB(ICBT,ISP)
     1              *((ONE-FRB)*EK(ITB+1,ISP)+FRB*EK(ITB+2,ISP))
        EDENT(ICBT)=EDENT(ICBT)+SDNT(ICBT,ISP)
     1              *((ONE-FRT)*EK(ITT+1,ISP)+FRT*EK(ITT+2,ISP))
  230   CONTINUE
      IF(NLTE.NE.0) THEN
        EDEEB(ICBT)=SDNB(ICBT,IELC)
     1    *((ONE-FRB)*EK(ITB+1,IELC)+FRB*EK(ITB+2,IELC))
        EDEET(ICBT)=SDNT(ICBT,IELC)
     1    *((ONE-FRT)*EK(ITT+1,IELC)+FRT*EK(ITT+2,IELC))
      END IF
  240 CONTINUE
      END IF
C========================================== INITIAL VALUES =====START
C     TE(IC) SHOULD BE INITICALIED WITH NON-ZERO VALUE FOR ALL CASES.
C     ROEE(IC) SHOULD BE INITIALIZED WITH SOME VALUE FOR NLTE=0
C---------------------------------------------------------------START
C     TRIPLE NESTED DO LOOP FOR INITIAL CONDITIONS - CONVINIENT
C---------------------------------------------------------------START
      DO 500 K=1,NZT
      ICK=(K-1)*NXYT
      DO 500 J=1,NYT
      ICJ=ICK+(J-1)*NXT
      DO 500 I=1,NXT
      IC=I+ICJ
      U(IC)=ZERO
      V(IC)=ZERO
c----- moving core initial condition
c     if(i.le.icin) then
c       v(ic)=vveld(i)*(yl-y(j))/yl
c     endif
c-----
      W(IC)=ZERO
      TEMP(IC)=TEMAMB
      TE(IC)=TEMAMB
      P(IC)=PAMB
      Q(IC)=RPGS2*P(IC)-PREF
      TKE(IC)=SMALL
      IF(ITURB.EQ.0) TKE(IC)=ZERO
      IF(ITURB.EQ.2) EPS(IC)=SMALL
  500 CONTINUE
      DO 510 ISP=1,LNRE
           DO 520 IC=1,NXYZT
           OMGDOT(IC,ISP)=ZERO
  520      CONTINUE
  510 CONTINUE
C---------------------------------------------------------------START
      DO 560 K=1,NZT
      ICK=(K-1)*NXYT
      DO 560 J=1,NYT
      ICJ=ICK+(J-1)*NXT
      DO 560 I=1,NXT
      IC=I+ICJ
c----- air environment
         xmfi(1)=0.00d0
         xmfi(2)=0.00d0
         xmfi(3)=0.00d0
         xmfi(4)=0.00d0
         xmfi(5)=0.00d0
         xmfi(6)=0.79d0
         xmfi(7)=0.00d0
         xmfi(8)=0.00d0
         xmfi(9)=0.00d0
         xmfi(10)=0.21d0
         xmfi(11)=0.00d0
         xmfi(12)=0.00d0
         xmfi(13)=0.00d0
         xmfi(14)=0.00d0
         xmfi(15)=0.00d0
c-----
      ROE(IC)=ZERO
      RO(IC)=ZERO
      ROEE(IC)=ZERO
      ROCV=ZERO
      DIV(IC)=ZERO
      IF(NLTE.EQ.0) CONDE(IC)=ZERO
      TB=HUNDTH*TEMP(IC)
      IT=INT(TB)
      FR=TB-DBLE(IT)
           DO 550 ISP=1,NSP
           SPD(IC,ISP)=XMFI(ISP)*MW(ISP)*P(IC)/(RGAS*TEMP(IC))
           RO(IC)=RO(IC)+SPD(IC,ISP)
           ROE(IC)=ROE(IC)+SPD(IC,ISP)*((ONE-FR)*EK(IT+1,ISP)
     1                                  +FR*EK(IT+2,ISP))
           ROCV=ROCV+HUNDTH*(EK(IT+2,ISP)-EK(IT+1,ISP))*SPD(IC,ISP)
  550      CONTINUE
      RON(IC)=RO(IC)
C---------------------------------------------------------------START
C     WE CAN USE GAMMA DEFINED HERE, WHEN NO ELECTRONS PRESENT
C     INITIALLY, OR LTE
C---------------------------------------------------------------START
      GAMMA(IC)=ONE+P(IC)/(ROCV*TEMP(IC))
      IF(NLTE.NE.0) THEN
      ROEE(IC)=SPD(IC,IELC)*((ONE-FR)*EK(IT+1,IELC)+FR*EK(IT+2,IELC))
      ELSE
        TE(IC)=TEMP(IC)
      END IF
  560 CONTINUE
C================ INITIAL VALUES OF PARTICLE MASS INJECTED =====START
      DO 600 INOZ=1,NUMNOZ
      TM1INJ(INOZ)=ZERO
  600 CONTINUE
C=========================== INITIALIZE CONTROL-PARAMETERS =====START
      NP=0
      NMJ=0
      TIME=ZERO
      TIMOVI=ZERO
C===============================================================START
      RETURN
      END
C
      SUBROUTINE RESETPWR(VOLTNEW)
C===============================================================START
C     THIS ROUTINE PROVIDES BOUNDARY CONDITIONS AND INITIAL VALUES
C     INITIALIZE EPS ONLY IF ITURB=2
C----------
C     CALLED BY LAVA
C===============================================================START
      INCLUDE 'COML.h'
C---------------------------------------------------------------START
C     TKE MUST BE SET TO BE ZERO FOR LAMINAR CASES
C     FOR BOTH INITIAL AND BOUNDARY CONDITIONS
C---------------------------------------------------------------START
C     IN CASE OF NLTE>=1 AND ELECTRON DENSITY IS NOT ZERO, SPECIAL
C     CODING IS REQUIRED TO DO THIS HERE AND BCCC. THIS SUBROUTINE IS
C     NOT SET UP TO DO SUCH CASES.
C---------------------------------------------------------------START
      INTEGER ITL,ITR,ITD,ITF,ITB,ITT
      INTEGER INOZ
      DOUBLE PRECISION XMFI(NSP),XMFL(NSP),XMFR(NSP),
     1     XMFD(NSP),XMFF(NSP),XMFB(NSP),XMFT(NSP),
     2     TB,FR,ROCV,FRL,FRR,FRD,FRF,FRB,FRT
C---------------------------------------------------------------START
      EXTERNAL CHEMEV
      INTEGER LNREV
      PARAMETER (LNREV=3)
      DOUBLE PRECISION EQC(LNREV),PROG(LNREV),BAM(LNREV,LNREV),
     1     RHS(LNREV),SPDSV(NSP),SPDNEG(NSP),
     2 DPROG(LNREV),DSPD(LNREV),WORK(5*LNREV),SWORK(LNREV),
     3  RWORK(LNREV),CWORK(LNREV),BAMF(LNREV,LNREV)
      DOUBLE PRECISION ASV(LNREV),BSV(LNREV),CSV(LNREV),DSV(LNREV),
     1     ESV(LNREV),FBNANV(NSP,LNREV)
      INTEGER NLMV(LNREV),ANV(NSP,LNREV),BNV(NSP,LNREV),CNV(NSP,LNREV)
      INTEGER IWORK(LNREV),ISKIP(LNREV),IPCON(LNREV),IPIV(LNREV),
     1     IREV,NK
C---------------------------------------------------------------START
C
      DOUBLE PRECISION SDN(NSP)
      DOUBLE PRECISION UDEL,DVDX,DVDXMX,TKMAX,DELTA
      INTEGER IDEL
C---------------------------------------------------------------START
      double precision vdotAr,vdotH2,xiH2,amssAr,amssH2
      double precision penv,poweri,amssin,amass,
     1     Hdotin,Hdtout,swirnn,wpeak,ucl,tcl,ctot,
     2     rhof,sumhf,rocvi,areav,spdsnd,power,vratio,tratio,
     3     wpxp,rhofr,rhofc,priter,axial,radial,rordr,sratio
      double precision tl,pold,pnew,pcur,res
      DOUBLE PRECISION VCRIT,VOLTNEW
      integer iter,ipeak,ngo
c--- icin for nozzle radius, icout for torch outside of METCO-9MB
c      data icin,icout /11,30/
      DATA VCRIT /0.01D0/
C===============================================================START
      WRITE(NFOUT,1000)TIME,VOLTNEW
 1000 FORMAT(5X,'-----------------------------------------------',/,
     1       5X,'NEW VOLTAGE AT T=',E12.3,'VOLT=',E12.3)
      IPBC=0
      VOLTG=VOLTNEW
C---------------------------------------- INPUT FOR CHEMEV -----START
C     DETERMINE ASV -- ESV, FBNANV,ANV,BNV,CNV,NLMV FOR CHEMEV
C---------------------------------------------------------------START
      DO 870 IREV=1,LNREV
C----------
      IF(irev.le.1) THEN
        ASV(IREV)=AKS(IREV)
        BSV(IREV)=BKS(IREV)
        CSV(IREV)=CKS(IREV)
        DSV(IREV)=DKS(IREV)
        ESV(IREV)=EKS(IREV)
      ELSE
        ASV(IREV)=AS(IREV-1)
        BSV(IREV)=BS(IREV-1)
        CSV(IREV)=CS(IREV-1)
        DSV(IREV)=DS(IREV-1)
        ESV(IREV)=ES(IREV-1)
      END IF
C----------
      NK=0
      DO 860 ISP=1,NSP
C----------
      IF(irev.le.1) THEN
        ANV(ISP,IREV)=AM(ISP,IREV)
        BNV(ISP,IREV)=BM(ISP,IREV)
      ELSE
        ANV(ISP,IREV)=AN(ISP,IREV-1)
        BNV(ISP,IREV)=BN(ISP,IREV-1)
         if(isp.eq.13) then
            ANV(ISP,IREV)=ANV(ISP,IREV)-1
            BNV(ISP,IREV)=BNV(ISP,IREV)-1
         endif
      END IF
C----------
      FBNANV(ISP,IREV)=DBLE(BNV(ISP,IREV)-ANV(ISP,IREV))
      IF(ANV(ISP,IREV).EQ.0 .AND. BNV(ISP,IREV).EQ.0) GO TO 860
      NK=NK+1
      CNV(NK,IREV)=ISP
  860 CONTINUE
      NLMV(IREV)=NK
  870 CONTINUE
C==================== conditions at the exit of the nozzle =====START
c--- THIS PART HAS BEEN CHANGED TO THE INPUT FILE
C     torch power = VOLTG * Amp * Efficiency * 1.0e7 (Joule to erg)
C---------------------------------------------------------------START
c*** if different pressure is desired, change in pamb in HOT
C      vdotAr=40.0d0
C      vdotH2=12.0d0
c reduce the flow rate see the convergence of the code
      xiH2=FLRT2/(FLRT1+FLRT2)
C       amp=600.0d0
C      VOLTG=40.0d0
C      eff=0.70d0
      penv=pamb/P1ATM
C---------------------------------------------------------------START
      PRINT*,        'TORCH OPERATING CONDITIONS ARE GIVEN AS'
      PRINT*,        CURCY,' Amp, ',VOLTG,' VOLTG, Eff = ',EFFCY
      PRINT*,        'Ar =',FLRT1,' standard l/min, ',penv,' ATM'
      PRINT*,        'H2 =',FLRT2,' standard l/min'
      PRINT*,        'SWIRL NUMBER = ',swirn
C--------------- mass flow rate g/s ----------------------------START
      poweri=1.0d7*VOLTG*CURCY*EFFCY
      vdotAr=FLRT1*1.0d3/60.0d0
      vdotH2=FLRT2*1.0d3/60.0d0
      amssAr=mw(1)*P1ATM/(Rgas*300.0d0)*VdotAr
      amssH2=mw(3)*P1ATM/(Rgas*300.0d0)*VdotH2
      amssin=amssAr+amssH2
C---------------------------------------------------------------START
       tb=hundth*temamb
       it=int(tb)
       fr=tb-dble(it)
      Hdotin=amssAr*((htform(1)+rgas*temamb)*rmw(1)
     1              +(one-fr)*ek(it+1,1)+fr*ek(it+2,1))
     2      +amssH2*((htform(3)+rgas*temamb)*rmw(3)
     3              +(one-fr)*ek(it+1,3)+fr*ek(it+2,3))
C---------------------------------------------------------------START
C      ipeak assumed to be Rin*0.8
C      SWIRL NUMBER =  9.999E-02 peak w= 44198.44959906817
C---------------------------------------------------------------START
      do 820 i=2,icin
      if(xc(i).gt.0.8d0*x(icin)) then
        ipeak=i-1
        go to 822
      endif
  820 continue
  822 continue
C------------------------------------------ ITERATION LOOP -----START
  900 continue
c      PRINT*,'Give powerv, powert, swirl number'
c      READ *,powerv,powert,swirn
c  ---- no read, through standard input by Y.P. WAN 7/2/97
C      powerv=20.0
C      powert=20.0
C      swirn=0.0
C---------------------------------------------------------------START
C     IF SWIRL PROFILE IS GIVEN, THEN SET SWIRN=0, AND GIVE SWIRL
C     PROFILE HERE
C---------------------------------------------------------------START
      if(swirn.eq.zero) swirnn=zero
      if(iswirl.eq.0 .and. swirn.ne.zero) call exita(1)
      if(iswirl.ne.0 .and. swirn.eq.zero) call exita(1)
      if(iswirl.eq.0) wpeak=zero
      if(iswirl.eq.1) wpeak=1000.0d0
C----------- why assume such a data ?   ------------------------START
      Ucl=178318.5307618468d0
      Tcl=13510.41156817594d0
  930 continue
C---------------------------------------------------------------START
      amass=zero
      Hdtout=zero
      do 750 i=2,icin
        tempd(i)=(tcl-700.0d0)*(one-(xc(i)/x(icin))**PWTEP)+700.0d0
        ctot=presd(i)/(Rgas*tempd(i))
C---------------------------------------------------------------START
C     get species densities using Newton's method and CHEMEV
C     initial guess needs to be consistent with all constraints.
C     (in array sdn).
C---------------------------------------------------------------START
c ---- this change of initial guessed values is done by Y.P. WAN 7/7/97
      if(i.eq.2) then
        sdn(1)= mw(1)*(one-xiH2)*ctot
        sdn(2)= 0.0d0
        sdn(3)= mw(3)*xiH2*ctot
        sdn(4)= 0.0d0
        sdn(5)= 0.0d0
        sdn(6)= 0.0d0
        sdn(7)= 0.0d0
        sdn(8)= 0.0d0
        sdn(9)= 0.0d0
        sdn(10)=0.0d0
        sdn(11)=0.0d0
        sdn(12)=0.0d0
        sdn(13)=0.0d0
        sdn(14)=0.0d0
        sdn(15)=0.0d0
      end if
C---------------------------------------------------------------START
      tl=tempd(i)
      pold=presd(i)
      iter=0
  700 continue
      iter=iter+1
      if(tl.gt.1250.0d0) then
      CALL CHEMEV(tl,sdn,NSP,LNREV,LNREV,
     1        ASV,BSV,CSV,DSV,ESV,ANV,BNV,CNV,NLMV,FBNANV,MW,RMW)
      endif
c      print*,i,iter,sdn(1),sdn(2),sdn(3),sdn(4),sdn(5)
c      print*,i,iter,tl
      ctot=zero
      do 720 isp=1,nsp
        ctot=ctot+sdn(isp)*rmw(isp)
  720 continue
      pcur=ctot*Rgas*tl
      res=pcur-presd(i)
      if(abs(res/presd(i)).lt.1.0d-6) go to 735
C----------------------------------------------- slope = 1 -----START
      pnew=pcur-res
      pnew=max(pnew,0.1d0*pcur)
      pnew=min(pnew,10.0d0*pcur)
C---------------------------------------------------------------START
        do 730 isp=1,nsp
          sdn(isp)=sdn(isp)*pnew/pcur
  730   continue
        pold=pcur
        pcur=pnew
      if(iter.gt.50) then
        write(*,*)'not converged, iter=',iter,', tl=',tl
        stop
      endif
        go to 700
C---------------------------------------------------------------START
  735 continue
C---------------------------------------------------------------START
        tb=hundth*tempd(i)
        it=int(tb)
        fr=tb-dble(it)
        edend(i)=zero
        rhof=zero
        sumhf=zero
        rocvi=zero
        do 740 isp=1,nsp
        sdnd(i,isp)=sdn(isp)
        edend(i)=edend(i)+sdnd(i,isp)*((1.0d0-fr)*ek(it+1,isp)
     1                                +fr*ek(it+2,isp))
        rhof=rhof+sdnd(i,isp)
        sumhf=sumhf+sdn(isp)*rmw(isp)*htform(isp)
        rocvi=rocvi+sdn(isp)*(ek(it+2,isp)-ek(it+1,isp))*hundth
  740   continue
        vveld(i)=ucl*(one-(xc(i)/x(icin))**PWVEL)
        areav=two*pi*xc(i)*dx(i)*vveld(i)
        amass=amass+rhof*areav
        Hdtout=Hdtout+areav*(edend(i)+presd(i)+sumhf)
c----- frozen speed of sound
        if(i.eq.2) spdsnd=hundth*sqrt((one+presd(i)/(rocvi*tempd(i)))
     1                                *presd(i)/rhof)
  750 continue
C---------------------------------------------------------------START
      power=Hdtout-Hdotin
c      vratio=amssin/amass
c      tratio=poweri/power
c ---- change the criteria for the convergence by Y.P. WAN 7/7/97
c      if(tratio.le.zero) then
c        write(*,*) 'ERROR!!!! --- negative tratio ....'
c        stop
c      endif
c      if(vratio.ge.0.9999d0 .and. vratio.le.1.0001d0 .and.
c     1   tratio.ge.0.9999d0 .and. tratio.le.1.0001d0) go to 960
c      ucl=vratio*ucl*0.25d0+0.75d0*ucl
c      tcl=tratio*tcl*0.10d0+0.90d0*tcl
      vratio=dabs((amssin-amass)/amssin)
      tratio=dabs((poweri-power)/poweri)
c -- the relative difference as criteria vcrit=1%
      if(vratio.le.vcrit.and.tratio.le.vcrit) go to 960
c      write(*,*)'****** ucl=',ucl,', tcl=',tcl,
c     & ', amssin=',amssin,',poweri=',poweri
c -- modify the values of ucl and tcl 
      ucl=amssin/amass*ucl*0.10d0+0.90d0*ucl
      tcl=poweri/power*tcl*0.10d0+0.90d0*tcl
c      write(*,*)'****** ucl=',ucl,', tcl=',tcl,
c     & ',amass=',amass, ',power=',power
      go to 930
  960 continue
C--------- success on the iteration of u and T      -----------START
      vveld(1)=vveld(2)
      tempd(1)=tempd(2)
      edend(1)=edend(2)
      teed(1)=teed(2)
      edeed(1)=edeed(2)
      do 770 isp=1,nsp
      sdnd(1,isp)=sdnd(2,isp)
  770 continue
C--------------------- INFLOW PROFILE FOR SGS OR K-E MODEL -----START
C    MIXING LENGTH AT INLET FOR K-E MODEL IS SPECIFIED USING
C    THE FORMULA SUGGESTED BY LESCHZINER AND RODI FOR SWIRL JET
C    (AIAA 22, P1742, 1984).
C---------------------------------------------------------------START
      IF(ITURB.NE.0) THEN
        UDEL=TENTH*VVELD(2)
        DVDXMX=SMALL
        DO 785 I=2,ICIN
          IF(VVELD(I).GE.UDEL) IDEL=I
          DVDX=ABS((VVELD(I)-VVELD(I-1))/(XC(I)-XC(I-1)))
          DVDXMX=MAX(DVDXMX,DVDX)
  785   CONTINUE
        DELTA=TWO*(XC(IDEL)+(UDEL-VVELD(IDEL))*
     1        (XC(IDEL+1)-XC(IDEL))/(VVELD(IDEL+1)-VVELD(IDEL)))
        TKMAX=HUNDTH*VVELD(2)*VVELD(2)
C----- reduced turbulence intensity
c       tkmax=0.001d0*VVELD(2)*VVELD(2)
        DO 780 I=3,ICIN
          DVDX=ABS((VVELD(I)-VVELD(I-1))/(XC(I)-XC(I-1)))
          IF(ITURB.EQ.1) THEN
C------ USER MAY WANT TO USE DIFFERENT VALUE FOR SGS MODEL -----START
            TKED(I)=THRHAF*TKMAX*DVDX/DVDXMX
            TKED(I)=MAX(TKED(I),SMALL)
          END IF
          IF(ITURB.EQ.2) THEN
            TKED(I)=THRHAF*TKMAX*DVDX/DVDXMX
            TKED(I)=MAX(TKED(I),SMALL)
            EPSD(I)=TKED(I)**1.5*CMU**0.75/(0.075D0*DELTA)
          END IF
  780   CONTINUE
        TKED(2)=TKED(3)
        EPSD(2)=EPSD(3)
        TKED(1)=TKED(2)
        EPSD(1)=EPSD(2)
      END IF
C--------------------------------- inflow profile of swirl -----START
c     This part needs to be provided by the user
C---------------------------------------------------------------START
      if(iswirl.eq.0) go to 990
  800 continue
      do 825 i=2,ipeak
      wveld(i)=xc(i)*wpeak/xc(ipeak)
  825 continue
      wpxp=wpeak/(x(icin)-xc(ipeak))
      do 827 i=ipeak+1,icin
      wveld(i)=wpeak-wpxp*(xc(i)-xc(ipeak))
  827 continue
      wveld(1)=-wveld(2)
C---------------------------------------------------------------START
c     pressure at the inside wall of nozzle should be gussed
C---------------------------------------------------------------START
      rhofr=pamb*mw(1)/(Rgas*700.0d0)
      rhofc=zero
       do 805 isp=1,nsp
       rhofc=rhofc+sdnd(icin,isp)
  805  continue
      presd(icin)=pamb-wveld(icin)**2*dx(icin)
     1                 *(rhofc+rhofr)/(eight*x(icin))
C---------------------------------------------------------------START
      do 810 i=icin-1,2,-1
        rhofr=zero
        rhofc=zero
         do 807 isp=1,nsp
         rhofr=rhofr+sdnd(i+1,isp)
         rhofc=rhofc+sdnd(i,isp)
  807    continue
        presd(i)=presd(i+1)+twoo3*(rhofr*tked(i+1)-rhofc*tked(i))
     1           -dx(i)*(rhofr+rhofc)*(wveld(i+1)+wveld(i))**2
     2            /(eight*x(i))
  810 continue
c----- presd(1) contains pressure of previous iteration
      priter=presd(1)/presd(2)
      if(priter.gt.0.9999d0 .and. priter.le.1.0001d0) go to 830
        do 840 i=2,icin
        presd(i)=presd(i)*(0.6d0+0.4d0*priter)
  840   continue
        presd(1)=presd(2)
        go to 930
  830 continue
      presd(1)=presd(2)
C---------------------------------------------------------------START
      axial=zero
      radial=zero
      do 850 i=2,icin
        rhof=zero
        do 845 isp=1,nsp
        rhof=rhof+sdnd(i,isp)
  845   continue
        rordr=rhof*xc(i)*dx(i)
        axial=axial+rordr*(vveld(i)**2-half*wveld(i)**2)
        radial=radial+rordr*vveld(i)*wveld(i)*xc(i)
  850 continue
      swirnn=radial/(axial*x(icin))
      if(swirn.eq.zero) go to 990
        sratio=swirn/swirnn
        if(sratio.gt.0.9999d0 .and. sratio.le.1.0001d0) go to 990
       print*,'---',swirnn,sratio,wpeak
        wpeak=wpeak*sratio
       print*,'---',swirnn,sratio,wpeak
        go to 800
  990 continue
C---------------------------------------------------------------START
c----- Swank's experiment - swirl is clockwise viewed from torch
c     do 995 i=1,icin
c     wveld(i)=-wveld(i)
c 995 continue
C---------------------------------------------------------------START
      power=power*1.d-10
      PRINT*,        'MASS= ',amass,'g/s, POWER= ',power,'KW'
      PRINT*,        'Ucl=',ucl,' cm/s, powerv',PWVEL
      PRINT*,        'Tcl=',tcl,' K, powert',PWTEP
      PRINT*,        'Frozen Speed of Sound= ',spdsnd,' m/sec'
      PRINT*,        'SWIRL NUMBER = ',swirnn,' peak w=',wpeak
      PRINT*,        'Pcl = ',presd(2),' ,Pamb = ',penv,' ATM'
c      PRINT*,'Satisfied? If yes, give 0; give 1 otherwise.'
      RETURN
      END
