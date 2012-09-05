*DECK SETQ
      SUBROUTINE SETQ
C================================================================SETQ
C     THIS ROUTINE PROVIDES PRESSURE "Q"
C     USING VARIOUS ACCELERATION TECHNIQUE
C----------
C     CALLED BY LAVA
C----------
C     INPUT: DIV, DIVE, P, GAMMA
C     OUTPUT: Q
C     PARAMETER: DAMP, MACCEL,THETA,RELAX,PQDIF
C================================================================SETQ
      INCLUDE 'COML.h'
C----------
      DOUBLE PRECISION DELDIV,ROD,COMB,BEETA,PSUM,VSUM,PSUBE
C
C========================================= NO ACCELERATION ======SETQ
      IF(MACCEL.EQ.0) THEN
        DO 10 IC=1,NXYZT
        DELDIV=DIV(IC)-DIVE(IC)
        Q(IC)=RPGS2*P(IC)*(ONE-DT(IC)*DAMP*DELDIV*GAMMA(IC))
   10   CONTINUE
C============== ACCELERATION METHOD 1 FOR COMPRESSIBLE FLOW =====SETQ
      ELSEIF(MACCEL.EQ.1) THEN
C----------------------------------------------------------------SETQ
C     IF WE DO NOT HAVE ANY PRESSURE BOUNDARY, PREF SHOULD BE VOLUME
C     AVERAGED PRESSURE OF THE DOMAIN. FOR CYLINDRICAL COORDINATES,
C     CIRCUMFERENTIAL VOLUME SHOULD BE USED.
C----------------------------------------------------------------SETQ
      IF(IPBC.EQ.0) THEN
        PSUM=ZERO
        VSUM=ZERO
        DO 130 K=K1,K2
        ICK=(K-1)*NXYT
          DO 120 J=J1,J2
          ICJ=(J-1)*NXT+ICK
            DO 110 I=I1,I2
            IC=I+ICJ
            PSUM=PSUM+P(IC)*VOL(IC)
            VSUM=VSUM+VOL(IC)
  110       CONTINUE
  120     CONTINUE
  130   CONTINUE
        PREF=PSUM/VSUM
      END IF
C----------------------------------------------------------------SETQ
      IF(NLTE.EQ.0) THEN
        DO 140 IC=IC1,IC2
        COMB=RO(IC)/(DT(IC)*(RDX(IC)*RDX(IC)+RDY(IC)*RDY(IC)
     1                      +RDZ(IC)*RDZ(IC)))
        BEETA=THETA*COMB/(DT(IC)*GAMMA(IC)*P(IC))
        BEETA=MIN(BEETA,ONE)
        ROD=DAMP*COMB
        Q(IC)=QN(IC)+BEETA*(P(IC)-PN(IC)-RELAX*(QN(IC)-PN(IC)+PREF))
     1              -ROD*(DIV(IC)-DIVOLD(IC))*DT(IC)/DTOLD(IC)
  140   CONTINUE
C---------------------------------------------- FOR NON-LTE -----SETQ
      ELSE
        DO 150 IC=IC1,IC2
        COMB=RO(IC)/(DT(IC)*(RDX(IC)*RDX(IC)+RDY(IC)*RDY(IC)
     1                      +RDZ(IC)*RDZ(IC)))
        PSUBE=SPDN(IC,IELC)*RGAS*RMW(IELC)*TE(IC)
        BEETA=THETA*COMB
     1        /(DT(IC)*(GAMMA(IC)*P(IC)+(FIVO3-GAMMA(IC))*PSUBE))
        BEETA=MIN(BEETA,ONE)
        ROD=DAMP*COMB
        Q(IC)=QN(IC)+BEETA*(P(IC)-PN(IC)-RELAX*(QN(IC)-PN(IC)+PREF))
     1              -ROD*(DIV(IC)-DIVOLD(IC))*DT(IC)/DTOLD(IC)
  150   CONTINUE
      END IF
C============ ACCELERATION METHOD 2 FOR INCOMPRESSIBLE FLOW =====SETQ
      ELSEIF(MACCEL.EQ.2) THEN
        DO 210 IC=IC1,IC2
        COMB=RO(IC)/(DT(IC)*(RDX(IC)*RDX(IC)+RDY(IC)*RDY(IC)
     1                      +RDZ(IC)*RDZ(IC)))
        ROD=DAMP*COMB
        Q(IC)=QN(IC)-THETA*COMB*(DIV(IC)-DIVE(IC))
     1              -ROD*(DIV(IC)-DIVOLD(IC))*DT(IC)/DTOLD(IC)
  210   CONTINUE
      END IF
C=========================================== AT BOUNDARIES ======SETQ
      IF(NX.NE.1) THEN
        DO 510 ICLR=1,NLR
        ICL=1+(ICLR-1)*NXT
        ICR=ICL+NX+1
        Q(ICL)=RPGS2*P(ICL)-PREF
        Q(ICR)=RPGS2*P(ICR)-PREF
  510   CONTINUE
      END IF
C----------
  600 IF(NY.NE.1) THEN
        DO 620 K=1,NZT
        ICK=(K-1)*NXYT
          DO 610 I=1,NXT
          ICD=I+ICK
          ICF=ICD+NXYT-NXT
          Q(ICD)=RPGS2*P(ICD)-PREF
          Q(ICF)=RPGS2*P(ICF)-PREF
  610     CONTINUE
  620   CONTINUE
      END IF
C----------
      IF(NZ.NE.1) THEN
        DO 710 ICBT=1,NBT
        ICB=ICBT
        ICT=ICBT+NXYZT-NXYT
        Q(ICB)=RPGS2*P(ICB)-PREF
        Q(ICT)=RPGS2*P(ICT)-PREF
  710   CONTINUE
      END IF
C================================================================SETQ
      RETURN
      END
