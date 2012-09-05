*DECK BCMOM
      SUBROUTINE BCMOM
C===============================================================BCMOM
C     THIS ROUTINE SETS VELOCITY BOUNDARY CONDITIONS
C----------
C     CALLED BY LAVA
C----------
C     INPUT:
C     OUTPUT:
C===============================================================BCMOM
      INCLUDE 'COML.h'
C----------
      DOUBLE PRECISION SWIRL
c      integer icin,icout,ictip,icdj,icdd
c--- icin for nozzle radius, icout for torch outside of METCO-9MB
c      data icin,icout /11,30/
C
C===============================================================BCMOM
      SWIRL=DBLE(ISWIRL)
C================== SET LEFT AND RIGHT BOUNDARY CONDITIONS =====BCMOM
      IF(NX.EQ.1) GO TO 100
      DO 10 ICLR=1,NLR
      ICL=1+(ICLR-1)*NXT
      ICR=ICL+NX
      U(ICL)=PBCL(ICLR)*U(ICL)+(ONE-PBCL(ICLR))*UVELL(ICLR)
c     V(ICL)=VVELL(ICLR)
      W(ICL)=SWIRL*W(ICL)+(ONE-SWIRL)*WVELL(ICLR)
c----- axis
      v(icl)=v(icl+1)
C---------------------------------------------------------------BCMOM
      U(ICR)=PBCR(ICLR)*U(ICR)+(ONE-PBCR(ICLR))*UVELR(ICLR)
      U(ICR+1)=U(ICR)
c     V(ICR+1)=VVELR(ICLR)
      W(ICR+1)=SWIRL*W(ICR+1)+(ONE-SWIRL)*WVELR(ICLR)
c----- entrainment boundary
      if(u(icr).gt.zero) then
        v(icr+1)=v(icr)
      else
        v(icr+1)=vvelr(iclr)
      endif
   10 CONTINUE
C============== SET DERRIERE AND FRONT BOUNDARY CONDITIONS =====BCMOM
  100 IF(NY.EQ.1) GO TO 200
      ICDF=0
      DO 120 K=1,NZT
      ICK=(K-1)*NXYT
           DO 110 I=1,NXT
           ICD=I+ICK
           ICF=ICD+NXYT-NXT-NXT
           ICDF=ICDF+1
           V(ICD)=PBCD(ICDF)*V(ICD)+(ONE-PBCD(ICDF))*VVELD(ICDF)
c          U(ICD)=UVELD(ICDF)
           W(ICD)=SWIRL*W(ICD)+(ONE-SWIRL)*WVELD(ICDF)
c----- inflow
      if(i.le.icin) then
        u(icd)=uveld(icdf)
      elseif(i.gt.icin .and. i.le.icout) then
c----- bc for the wall is done at WALLF
        u(icd)=u(icd+nxt)
      elseif(i.gt.icout) then
        if(v(icd).le.zero) then
          u(icd)=u(icd+nxt)
        else
          u(icd)=uveld(icdf)
        endif
      endif
C---------------------------------------------------------------BCMOM
           V(ICF)=PBCF(ICDF)*V(ICF)+(ONE-PBCF(ICDF))*VVELF(ICDF)
           V(ICF+NXT)=V(ICF)
c          U(ICF+NXT)=UVELF(ICDF)
           W(ICF+NXT)=SWIRL*W(ICF+NXT)+(ONE-SWIRL)*WVELF(ICDF)
c----- outflow - open boundary
      if(v(icf).gt.zero) then
        u(icf+nxt)=u(icf)
c       v(icf-nxt)=v(icf-2*nxt)
c       v(icf)=v(icf-nxt)
      else
        u(icf+nxt)=uvelf(icdf)
      endif
  110      CONTINUE
  120 CONTINUE
C================== SET TOP AND BOTTOM BOUNDARY CONDITIONS =====BCMOM
  200 IF(NZ.EQ.1) GO TO 300
      DO 210 ICBT=1,NBT
      ICB=ICBT
      ICT=ICBT+NXYZT-NXYT-NXYT
      W(ICB)=PBCB(ICBT)*W(ICB)+(ONE-PBCB(ICBT))*WVELB(ICBT)
      U(ICB)=UVELB(ICBT)
      V(ICB)=VVELB(ICBT)
      W(ICT)=PBCT(ICBT)*W(ICT)+(ONE-PBCT(ICBT))*WVELT(ICBT)
      W(ICT+NXYT)=W(ICT)
      U(ICT+NXYT)=UVELT(ICBT)
      V(ICT+NXYT)=VVELT(ICBT)
  210 CONTINUE
C===============================================================BCMOM
  300 CONTINUE
C========================================== BLOCKAGE LOGIC =====BCMOM
C     PUT BLOCKAGE LOGIC HERE. BOUNDARY CONDITIONS ALONG THIS WALL
C          SHOULD BE PROVIDED AT THE CORRESPONDING BOUNDARY CONDITION
C          LOGIC.
C---------------------------------------------------------------BCMOM
C----- FREE SLIP ALONG THE BLOCKAGE
c     do 400 j=1,26
c     ICJ=(J-1)*NXT
c     do 400 i=icin,icout
c     IC=I+ICJ
c     IF(RR(IC).LE.SMALL) U(IC)=ZERO
c     IF(RF(IC).LE.SMALL) V(IC)=ZERO
cc    if(rt(ic).le.small) w(ic)=zero
c 400 CONTINUE
C===============================================================BCMOM
      RETURN
      END
