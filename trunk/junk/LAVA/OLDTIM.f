*DECK OLDTIM
      SUBROUTINE OLDTIM
C==============================================================OLDTIM
C     THIS ROUTINE SAVES OLD TIME VALUES,
C     CALCULATES DIV(IC) AND DIVOLD(IC)
C----------
C     CALLED BY LAVA
C----------
C     INPUT: SPD,UN,VN,WN,RR,RF,RT,RC,U,V,W,ROE,P,Q,
C            DX,DY,DZ,RDX,RDY,RDZ
C     OUTPUT: SPDN,RON,DIV,DIVOLD,DIVE,UN,VN,WN,ROEN,PN,QN
C==============================================================OLDTIM
      INCLUDE 'COML.h'
C==============================================================OLDTIM
      DO 100 IC=1,NXYZT
      RON(IC)=ZERO
      DIVOLD(IC)=DIV(IC)
      DIVE(IC)=ZERO
  100 CONTINUE
C----------- SAVE OLD SPECIES DENSITIES AND TOTAL DENSITY -----OLDTIM
      DO 200 ISP=1,NSP
           DO 220 IC=1,NXYZT
           SPDN(IC,ISP)=SPD(IC,ISP)
           RON(IC)=RON(IC)+SPD(IC,ISP)
  220      CONTINUE
  200 CONTINUE
C------------------------------------ SET OLD TIME ARRAYS -----OLDTIM
C     WSDT IS SET ZERO FOR FUTURE USE IN SUBROUTINE PINT
C--------------------------------------------------------------OLDTIM
      DO 300 IC=1,NXYZT
      UN(IC)=U(IC)
      VN(IC)=V(IC)
      WN(IC)=W(IC)
      ROEN(IC)=ROE(IC)
      ROEEN(IC)=ROEE(IC)
      PN(IC)=P(IC)
      QN(IC)=Q(IC)
      TKEN(IC)=TKE(IC)
      EPSN(IC)=EPS(IC)
      WSDT(IC)=ZERO
      DIV(IC)=ZERO
      UDELR(IC)=ZERO
  300 CONTINUE
C========================== CALCULATE DIVERGENCE AND UDELR ====OLDTIM
C     TRIPLE NESTED DO LOOP IS NECESSAARY TO AVOID BOUNDARY CELLS
C     DIVERGENCE SHOULD NOT BE CACULATED FOR BOUNDARY CELLS
C     UDELR AT BOUNDARY CELL IS SET TO BE THE SAME AS THE ONE INSIDE
C--------------------------------------------------------------OLDTIM
      IF(NX.NE.1) THEN
        DO 400 K=K1,K2
        ICK=(K-1)*NXYT
        DO 400 J=J1,J2
        ICJ=ICK+(J-1)*NXT
        DO 400 I=I1,I2
        IC=ICJ+I
        DIV(IC)=DIV(IC)+(RR(IC)*U(IC)-RR(IC-1)*U(IC-1))*RDX(IC)
        UDELR(IC)=UDELR(IC)
     1            +HALF*(U(IC)+U(IC-1))*(RR(IC)-RR(IC-1))*RDX(IC)
  400   CONTINUE
        DO 410 ICLR=1,NLR
        ICL=1+(ICLR-1)*NXT
        ICR=ICL+NX+1
        UDELR(ICL)=UDELR(ICL+1)
        UDELR(ICR)=UDELR(ICR-1)
  410   CONTINUE
      END IF
C----------
      IF(NY.NE.1) THEN
        DO 500 K=K1,K2
        ICK=(K-1)*NXYT
        DO 500 J=J1,J2
        ICJ=ICK+(J-1)*NXT
        DO 500 I=I1,I2
        IC=ICJ+I
        DIV(IC)=DIV(IC)+(RF(IC)*V(IC)-RF(IC-NXT)*V(IC-NXT))*RDY(IC)
        UDELR(IC)=UDELR(IC)
     1            +HALF*(V(IC)+V(IC-NXT))*(RF(IC)-RF(IC-NXT))*RDY(IC)
  500   CONTINUE
        DO 520 K=1,NZT
        ICK=(K-1)*NXYT
          DO 510 I=1,NXT
          ICD=I+ICK
          ICF=ICD+NXYT-NXT
          UDELR(ICD)=UDELR(ICD+NXT)
          UDELR(ICF)=UDELR(ICF-NXT)
  510     CONTINUE
  520   CONTINUE
      END IF
C----------
      IF(NZ.NE.1) THEN
        DO 600 K=K1,K2
        ICK=(K-1)*NXYT
        DO 600 J=J1,J2
        ICJ=ICK+(J-1)*NXT
        DO 600 I=I1,I2
        IC=ICJ+I
        DIV(IC)=DIV(IC)+(RT(IC)*W(IC)-RT(IC-NXYT)*W(IC-NXYT))*RDZ(IC)
        UDELR(IC)=UDELR(IC)
     1          +HALF*(W(IC)+W(IC-NXYT))*(RT(IC)-RT(IC-NXYT))*RDZ(IC)
  600   CONTINUE
        DO 610 ICBT=1,NBT
        ICB=ICBT
        ICT=ICBT+NXYZT-NXYT
        UDELR(ICB)=UDELR(ICB+NXYT)
        UDELR(ICT)=UDELR(ICT-NXYT)
  610   CONTINUE
      END IF
C----------
      DO 700 IC=IC1,IC2
      DIV(IC)=DIV(IC)*RRC(IC)
      UDELR(IC)=UDELR(IC)*RRC(IC)
  700 CONTINUE
C==============================================================OLDTIM
      RETURN
      END
