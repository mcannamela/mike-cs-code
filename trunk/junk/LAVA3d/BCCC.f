*DECK BCCC
      SUBROUTINE BCCC
C================================================================BCCC
C     THIS ROUTINE SETS BOUNDARY CONDITIONS ON CELL-CENTERED
C     QUANTITIES
C----------
C     CALLED BY LAVA
C----------
C     INPUT:
C     OUTPUT:
C================================================================BCCC
      INCLUDE 'COML.h'
C----------
      DOUBLE PRECISION SWIRL,SPDI,ROEI
C
c      integer icin,icout
c--- icin for nozzle radius, icout for torch outside of METCO-9MB
c ---- this data should be put into input data file
c      data icin,icout /11,30/
C================================================================BCCC
C
C     AT PRESENT THESE ARE SET HALF A CELL OUTSIDE THE BOUNDARY
C
C================================================================BCCC
      SWIRL=DBLE(ISWIRL)
C=================== SET LEFT AND RIGHT BOUNDARY CONDITIONS =====BCCC
      IF(NX.EQ.1) GO TO 100
      DO 10 ICLR=1,NLR
      ICL=1+(ICLR-1)*NXT
      ICR=ICL+NX+1
      RO(ICL)=ZERO
      RO(ICR)=ZERO
c     P(ICL)=PRESL(ICLR)
c     ROE(ICL)=EDENL(ICLR)
c     TEMP(ICL)=TEMPL(ICLR)
c     ROEE(ICL)=EDEEL(ICLR)
c     TE(ICL)=TEEL(ICLR)
c     TKE(ICL)=TKEL(ICLR)
c     IF(ITURB.EQ.2) EPS(ICL)=EPSL(ICLR)
      W(ICL)=SWIRL*WVELL(ICLR)+(ONE-SWIRL)*W(ICL)
c----- sym. axis
      p(icl)=p(icl+1)
      q(icl)=q(icl+1)
      roe(icl)=roe(icl+1)
      temp(icl)=temp(icl+1)
      roee(icl)=roee(icl+1)
      te(icl)=te(icl+1)
      tke(icl)=tke(icl+1)
      if(iturb.eq.2) eps(icl)=eps(icl+1)
C----------------------------------------------------------------BCCC
      P(ICR)=PRESR(ICLR)
      if(u(icr-1).lt.zero) then
      ROE(ICR)=EDENR(ICLR)
      TEMP(ICR)=TEMPR(ICLR)
c     ROEE(ICR)=EDEER(ICLR)
c     TE(ICR)=TEER(ICLR)
        roee(icr)=roee(icr-1)
        te(icr)=te(icr-1)
      TKE(ICR)=TKER(ICLR)
      IF(ITURB.EQ.2) EPS(ICR)=EPSR(ICLR)
      W(ICR)=SWIRL*WVELR(ICLR)+(ONE-SWIRL)*W(ICR)
      else
        roe(icr)=roe(icr-1)
        temp(icr)=temp(icr-1)
        roee(icr)=roee(icr-1)
        te(icr)=te(icr-1)
        tke(icr)=tke(icr-1)
        if(iturb.eq.2) eps(icr)=eps(icr-1)
        w(icr)=w(icr-1)
      endif
   10 CONTINUE
C----------------------------------------------------------------BCCC
      DO 50 ISP=1,NSP
           DO 30 ICLR=1,NLR
           ICL=1+(ICLR-1)*NXT
           ICR=ICL+NX+1
c          SPD(ICL,ISP)=SDNL(ICLR,ISP)
c----- sym. axis
       spd(icl,isp)=spd(icl+1,isp)
c-----open boundary
       if(u(icr-1).ge.zero) then
         spd(icr,isp)=spd(icr-1,isp)
       else
           SPD(ICR,ISP)=SDNR(ICLR,ISP)
       endif
           RO(ICL)=RO(ICL)+SPD(ICL,ISP)
           RO(ICR)=RO(ICR)+SPD(ICR,ISP)
   30      CONTINUE
   50 CONTINUE
C---- BERNOULLI BOUNDARY CONDITION FOR OPEN INFLOW BOUNDARY -----BCCC
      DO 70 ICLR=1,NLR
        ICL=1+(ICLR-1)*NXT
        ICR=ICL+NX+1
        IF(U(ICL).GT.ZERO) P(ICL)=PRESL(ICLR)-PBCL(ICLR)*HALF*
     1                                        RO(ICL)*U(ICL)*U(ICL)
        IF(U(ICR-1).LT.ZERO) P(ICR)=PRESR(ICLR)-PBCR(ICLR)*HALF*
     1                                      RO(ICR)*U(ICR-1)*U(ICR-1)
   70 CONTINUE
C=============== SET DERRIERE AND FRONT BOUNDARY CONDITIONS =====BCCC
  100 IF(NY.EQ.1) GO TO 200
      ICDF=0
      DO 120 K=1,NZT
      ICK=(K-1)*NXYT
        DO 110 I=1,NXT
        ICD=I+ICK
        ICF=ICD+NXYT-NXT
        ICDF=ICDF+1
        RO(ICD)=ZERO
        RO(ICF)=ZERO
        P(ICD)=PRESD(ICDF)
c----- inflow boundary
      if(i.le.icin) then
        ROE(ICD)=EDEND(ICDF)
        TEMP(ICD)=TEMPD(ICDF)
        ROEE(ICD)=EDEED(ICDF)
        TE(ICD)=TEED(ICDF)
        TKE(ICD)=TKED(ICDF)
        IF(ITURB.EQ.2) EPS(ICD)=EPSD(ICDF)
        W(ICD)=SWIRL*WVELD(ICDF)+(ONE-SWIRL)*W(ICD)
c----- torch face, BC done in WALLF.f
      else if(i.gt.icin .and. i.le.icout) then
          roe(icd)=roe(icd+nxt)
          temp(icd)=temp(icd+nxt)
          roee(icd)=roee(icd+nxt)
          te(icd)=te(icd+nxt)
          tke(icd)=tke(icd+nxt)
          if(iturb.eq.2) eps(icd)=eps(icd+nxt)
c----- free-slip for swirl -- no WALLF for swirl yet
          w(icd)=w(icd+nxt)
c----- open boundary
      else if(i.gt.icout) then
        if(v(icd).gt.zero) then
          roe(icd)=edend(icdf)
          temp(icd)=tempd(icdf)
          roee(icd)=roee(icd+nxt)
          te(icd)=te(icd+nxt)
          tke(icd)=tked(icdf)
          if(iturb.eq.2) eps(icd)=epsd(icdf)
          w(icd)=swirl*wveld(icdf)+(one-swirl)*w(icd)
        else
          roe(icd)=roe(icd+nxt)
          temp(icd)=temp(icd+nxt)
          roee(icd)=roee(icd+nxt)
          te(icd)=te(icd+nxt)
          tke(icd)=tke(icd+nxt)
          if(iturb.eq.2) eps(icd)=eps(icd+nxt)
          w(icd)=w(icd+nxt)
        end if
      end if
C----------------------------------------------------------------BCCC
c----- open boundary
        P(ICF)=PRESF(ICDF)
      if(v(icf-nxt).ge.zero) then
        roe(icf)=roe(icf-nxt)
        temp(icf)=temp(icf-nxt)
        roee(icf)=roee(icf-nxt)
        te(icf)=te(icf-nxt)
        tke(icf)=tke(icf-nxt)
        if(iturb.eq.2) eps(icf)=eps(icf-nxt)
        w(icf)=w(icf-nxt)
      else
        ROE(ICF)=EDENF(ICDF)
        TEMP(ICF)=TEMPF(ICDF)
c       ROEE(ICF)=EDEEF(ICDF)
c       TE(ICF)=TEEF(ICDF)
        roee(icf)=roee(icf-nxt)
        te(icf)=te(icf-nxt)
        TKE(ICF)=TKEF(ICDF)
        IF(ITURB.EQ.2) EPS(ICF)=EPSF(ICDF)
        W(ICF)=SWIRL*WVELF(ICDF)+(ONE-SWIRL)*W(ICF)
      end if
  110      CONTINUE
  120 CONTINUE
C----------------------------------------------------------------BCCC
      DO 150 ISP=1,NSP
      ICDF=0
        DO 140 K=1,NZT
        ICK=(K-1)*NXYT
          DO 130 I=1,NXT
          ICD=I+ICK
          ICF=ICD+NXYT-NXT
          ICDF=ICDF+1
c         SPD(ICD,ISP)=SDND(ICDF,ISP)
c         SPD(ICF,ISP)=SDNF(ICDF,ISP)
c----- inflow boundary
      if(i.le.icin) then
c------ spd at inflow is set
          spd(icd,isp)=sdnd(icdf,isp)
c----- torch face, BC done in WALLF.f
      elseif(i.gt.icin .and. i.le.icout) then
          spd(icd,isp)=spd(icd+nxt,isp)
c----- open boundary
      elseif(i.gt.icout) then
        if(v(icd).gt.zero) then
          spd(icd,isp)=sdnd(icdf,isp)
        else
          spd(icd,isp)=spd(icd+nxt,isp)
        endif
      endif
c----- open boundary
      if(v(icf-nxt).ge.zero) then
        spd(icf,isp)=spd(icf-nxt,isp)
      else
        spd(icf,isp)=sdnf(icdf,isp)
      endif
          RO(ICD)=RO(ICD)+SPD(ICD,ISP)
          RO(ICF)=RO(ICF)+SPD(ICF,ISP)
  130     CONTINUE
  140   CONTINUE
  150 CONTINUE
C---- BERNOULLI BOUNDARY CONDITION FOR OPEN INFLOW BOUNDARY -----BCCC
      ICDF=0
      DO 180 K=1,NZT
      ICK=(K-1)*NXYT
        DO 170 I=1,NXT
        ICD=I+ICK
        ICF=ICD+NXYT-NXT
        ICDF=ICDF+1
        IF(V(ICD).GT.ZERO)
     1   P(ICD)=PRESD(ICDF)-PBCD(ICDF)*HALF*RO(ICD)*V(ICD)*V(ICD)
        IF(V(ICF-NXT).LT.ZERO) P(ICF)=PRESF(ICDF)
     1                 -PBCF(ICDF)*HALF*RO(ICF)*V(ICF-NXT)*V(ICF-NXT)
  170   CONTINUE
  180 CONTINUE
C=================== SET TOP AND BOTTOM BOUNDARY CONDITIONS =====BCCC
  200 IF(NZ.EQ.1) GO TO 300
      DO 210 ICBT=1,NBT
      ICB=ICBT
      ICT=ICBT+NXYZT-NXYT
      RO(ICB)=ZERO
      RO(ICT)=ZERO
      P(ICB)=PRESB(ICBT)
      ROE(ICB)=EDENB(ICBT)
      TEMP(ICB)=TEMPB(ICBT)
      ROEE(ICB)=EDEEB(ICBT)
      TE(ICB)=TEEB(ICBT)
      TKE(ICB)=TKEB(ICBT)
      IF(ITURB.EQ.2) EPS(ICB)=EPSB(ICBT)
C----------------------------------------------------------------BCCC
      P(ICT)=PREST(ICBT)
      ROE(ICT)=EDENT(ICBT)
      TEMP(ICT)=TEMPT(ICBT)
      ROEE(ICT)=EDEET(ICBT)
      TE(ICT)=TEET(ICBT)
      TKE(ICT)=TKET(ICBT)
      IF(ITURB.EQ.2) EPS(ICT)=EPST(ICBT)
  210 CONTINUE
C----------------------------------------------------------------BCCC
      DO 250 ISP=1,NSP
        DO 230 ICBT=1,NBT
        ICB=ICBT
        ICT=ICBT+NXYZT-NXYT
        SPD(ICB,ISP)=SDNB(ICBT,ISP)
        SPD(ICT,ISP)=SDNT(ICBT,ISP)
        RO(ICB)=RO(ICB)+SPD(ICB,ISP)
        RO(ICT)=RO(ICT)+SPD(ICT,ISP)
  230   CONTINUE
  250 CONTINUE
C---- BERNOULLI BOUNDARY CONDITION FOR OPEN INFLOW BOUNDARY -----BCCC
      DO 270 ICBT=1,NBT
        ICB=ICBT
        ICT=ICBT+NXYZT-NXYT
        IF(W(ICB).GT.ZERO)
     1   P(ICB)=PRESB(ICBT)-PBCB(ICBT)*HALF*RO(ICB)*W(ICB)*W(ICB)
        IF(W(ICT-NXYT).LT.ZERO) P(ICT)=PREST(ICBT)
     1               -PBCT(ICBT)*HALF*RO(ICT)*W(ICT-NXYT)*W(ICT-NXYT)
  270 CONTINUE
C================================================================BCCC
  300 CONTINUE
C========================================== BLOCKAGE LOGIC ======BCCC
C     PUT BLOCKAGE LOGIC HERE. BOUNDARY CONDITIONS ALONG THIS WALL
C          SHOULD BE PROVIDED AT THE CORRESPONDING BOUNDARY CONDITION
C          LOGIC.
C----------------------------------------------------------------BCCC
C----- ADIABATIC ALONG THE BLOCKAGE
C     SPDI=PAMB*MW(1)/(RGAS*TEMAMB)
C     ROEI=SPDI*EK(4,1)
c     do 400 j=1,26
c     icj=(j-1)*nxt
c     do 400 i=icin+1,icout
c     ic=i+icj
C     IF(RC(IC).LE.SMALL) THEN
c       ROE(IC)=ROEI
c       TEMP(IC)=TEMAMB
c       ROEE(IC)=SMALL
c       TE(IC)=TEMAMB
c       P(IC)=PAMB
c       Q(IC)=RPGS2*P(IC)-PREF
c       TKE(IC)=SMALL
c       IF(ITURB.EQ.2) EPS(IC)=SMALL
c-----
c       roe(ic)=edend(i)
c       temp(ic)=tempd(i)
c       roee(ic)=edeed(i)
c       te(ic)=teed(i)
c       p(ic)=presd(i)
c       q(ic)=rpgs2*p(ic)-pref
c       tke(ic)=tked(i)
c       eps(ic)=epsd(i)
c       w(icd)=wveld(i)
c-----
c       RO(IC)=ZERO
c       DO 410 ISP=1,NSP
cc        SPD(IC,ISP)=SPDI
c         spd(ic,isp)=sdnd(i,isp)
c         RO(IC)=RO(IC)+SPD(IC,ISP)
c 410   CONTINUE
C     END IF
c 400 CONTINUE
C================================================================BCCC
      RETURN
      END
