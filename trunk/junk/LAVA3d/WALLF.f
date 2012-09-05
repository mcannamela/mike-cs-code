*DECK WALLF
      SUBROUTINE WALLF
C===============================================================WALLF
C     THIS ROUTINE UPDATES VELOCITIES AND OTHER SCALAR VARIABLES
C     BY TURBULENCE WALL FUNCTION
C----------
C     THIS ROUTINE UPDATES VELOCITIES AND ENERGY IMPLICITELY, AND
C     SPECIES DENSITIES EXPLICITELY.  SO SHOULD BE LOCATED BEFORE
C     CHMIMP.  BECAUSE BOUNDARY CONDITIONS FOR TKE AND EPS ARE SET
C     HERE (ROTKER AND ROEPSR), THIS SHOULD BE LOCATED AFTER PDVETC.
C     THIS COULD MAKE SOME CONFLICT IN NON-LTE CASE WHICH INVLOVES
C     IMPLICIT UPDATE OF ROER AND ROEER IN PDVETC.
C     MAY BE SOLVED BY UPDATE ROER AND ROEER IMPLICITELY HERE.
C     ALSO MIXING EXPLICIT-IMPICIT UPDATES MAY CAUSE INACCURATE
C     STEADY STATE SOLUTION AND TIME STEP RESTRICTION.
C----------
C     USER HAS TO SUPPLY THIS ROUTINE
C----------
C     CALLED BY LAVA
C===============================================================WALLF
      INCLUDE 'COML.h'
C----------
C
      DOUBLE PRECISION TWALL,SPDW(NSP),DIS,AREAM,ARV,
     1     MASSC,UCMX,UCMY,UPA,UPAR,VISL,REYNOL,FU,TAUW,
     2     DELTPP,PPAR,DELTAP,DELPX,DELPY,UAV,VAV,DPRPAR,
     3     DCOL(NSP),DSPDW(NSP),BDJ,SUMFL,SUMFQ,SUMQD,DIFW(NSP),
     4     SUMHD,CONDL,BHD,RCVRO,TFILM,TB,FR
      INTEGER IW
      double precision cAr,cH2,cN2,cO2,rctot,cwall,xAr,xH2,xN2,xO2
c--- icin for nozzle radius, icout for torch outside of METCO-9MB
c      integer icin,icout
c      data icin,icout /11,30/
C===============================================================WALLF
c----- along the torch face
      DO 100 i=icin+1,icout
      IW=nxt+i
c----------change
      TWALL=tempd(i)
       cAr=(spd(ic,1)+spd(ic,2))*rmw(1)
       cH2=(spd(iw,3)+spd(iw,4)+spd(iw,5))*rmw(3)
       cN2=(spd(iw,6)+spd(iw,7)+spd(iw,8)+spd(iw,9))*rmw(6)
       cO2=(spd(iw,10)+spd(iw,11)+spd(iw,12))*rmw(10)
       rctot=one/(cAr+cH2+cN2+cO2)
       xAr=cAr*rctot
       xH2=cH2*rctot
       xN2=cN2*rctot
       xO2=cO2*rctot
       cwall=p(iw)/(Rgas*twall)
       spdw(1)=xAr*cwall*mw(1)
       spdw(2)=zero
       spdw(3)=xH2*cwall*mw(3)
       spdw(4)=zero
       spdw(5)=zero
       spdw(6)=xN2*cwall*mw(6)
       spdw(7)=zero
       spdw(8)=zero
       spdw(9)=zero
       spdw(10)=xO2*cwall*mw(10)
       spdw(11)=zero
       spdw(12)=zero
       spdw(13)=zero
       spdw(14)=zero
       spdw(15)=zero
c----------
      dis=half*dy(iw)
      aream=rc(iw)*dx(iw)
      arv=one/dy(iw)
c----------
      massc=ron(iw)*dx(iw)*dy(iw)*rc(iw)
      upa=half*(un(iw)+un(iw-1))
      upar=abs(upa)
c----------no change
      VISL=VISC(IW)-VIST(IW)
      REYNOL=RON(IW)*DIS*UPAR/VISL
        IF(REYNOL.EQ.ZERO) THEN
          FU=ZERO
        ELSEIF(REYNOL.GT.ZERO .AND. REYNOL.LT.REYC) THEN
          FU=ONE/REYNOL
        ELSEIF(REYNOL.GE.REYC) THEN
          FU=ONE/(0.75D0+2.19D0*LOG(REYNOL))**2
        END IF
      TAUW=RON(IW)*UPAR*UPAR*FU
C------------------------------------FOR MOMENTUM (CHANGE) -----WALLF
      deltpp=-tauw*aream*dt(iw)*upa/(upar+small)
      ppar=massc*(upa+small)
      dprpar=one/(one-deltpp/ppar)
      uav=un(iw-1)+un(iw)
      if(uav.eq.zero) then
        deltap=zero
      else
        deltap=deltpp*dprpar/uav
      endif
c-----
      roru(iw-1)=roru(iw-1)
     1               +deltap*un(iw-1)/(dy(iw)*half*(dx(iw-1)+dx(iw)))
      roru(iw)=roru(iw)+deltap*un(iw)/(dy(iw)*half*(dx(iw+1)+dx(iw)))
C--------------------------- FOR K AND EPSILON (NO CHANGE) -----WALLF
      ROTKER(IW)=ROTKER(IW)+HALF*UPAR*UPAR*(ONE-DPRPAR*DPRPAR)
     1                      *RON(IW)*RC(IW)
      ROEPSR(IW)=0.4108D0*ROTKER(IW)
     1           *SQRT(ROTKER(IW)/(RON(IW)*RC(IW)))/DIS
C--------------------------------------------- FOR SPECIES -----WALLF
        DO 110 ISP=1,NSP
        DCOL(ISP)=DCOEF(IW,ISP)-RSCT*VIST(IW)/RON(IW)
        DSPDW(ISP)=SPDW(ISP)-SPDN(IW,ISP)
  110   CONTINUE
      BDJ=TAUW/(VISL*UPAR+SMALL)
      SUMFL=ZERO
      SUMFQ=ZERO
        DO 120 ISP=1,NSP
        IF(ISP.NE.IELC) THEN
          SUMFL=SUMFL+BDJ*DCOL(ISP)*DSPDW(ISP)
          SUMFQ=SUMFQ+MW(ISP)*QOM(ISP)*SPDN(IW,ISP)*DCOL(ISP)
        END IF
  120   CONTINUE
      SUMQD=ZERO
        DO 130 ISP=1,NSP
        IF(ISP.NE.IELC) THEN
          DIFW(ISP)=-BDJ*DCOL(ISP)*DSPDW(ISP)
     1            +SUMFL*SPDN(IW,ISP)/RON(IW)
     2      +BDJ*DSPDW(IELC)*(MW(ISP)*QOM(ISP)*SPDN(IW,ISP)*DCOL(ISP)
     3                           -SUMFQ*SPDN(IW,ISP)/RON(IW))
     4               /(QOM(IELC)*SPDN(IW,IELC)*MW(IELC)+SMALL)
          SUMQD=SUMQD+QOM(ISP)*DIFW(ISP)
        END IF
  130   CONTINUE
        DIFW(IELC)=-SUMQD/QOM(IELC)
        SUMHD=ZERO
        TFILM=HALF*(TEMP(IW)+TWALL)
        TB=HUNDTH*TFILM
        IT=INT(TB)
        IT=MIN(IT,NIT1)
        FR=TB-DBLE(IT)
        DO 140 ISP=1,NSP
        SPDR(IW,ISP)=SPDR(IW,ISP)-DT(IW)*DIFW(ISP)*ARV*RC(IW)
        SUMHD=SUMHD+DIFW(ISP)*((ONE-FR)*EK(IT+1,ISP)+FR*EK(IT+2,ISP)
     1        +RGAS*TFILM*RMW(ISP))
  140   CONTINUE
C---------------------------------------------- FOR ENERGY -----WALLF
      CONDL=COND(IW)-VIST(IW)*RPRT*GAMMA(IW)*P(IW)
     1               /((GAMMA(IW)-ONE)*RON(IW)*TEMP(IW))
      RCVRO=(GAMMA(IW)-ONE)*TEMP(IW)/P(IW)
      BHD=ARV*DT(IW)*RC(IW)*CONDL*TAUW/(VISL*UPAR+SMALL)
      ROER(IW)=(ROER(IW)-ARV*DT(IW)*RC(IW)*SUMHD
     1          +BHD*RCVRO*ROEN(IW)-BHD*(TEMP(IW)-TWALL))
     2        /(ONE+BHD*RCVRO)
  100 CONTINUE
C===============================================================WALLF
      RETURN
      END
