*DECK TIMSTP
      SUBROUTINE TIMSTP
C==============================================================TIMSTP
C     THIS ROUTINE CONTROLS THE TIME-STEP AUTOMATICALLY
C     INPUT PARAMETERS ARE LDT AND SAFEDT
C     LDT=0 : DT = CONST
C     LDT=1 : DT = F(TIME)
C     LDT=2 : DT = F(TIME,X,Y,Z), but not available.
C----------
C     CALLED BY LAVA
C==============================================================TIMSTP
      INCLUDE 'COML.h'
C----------
      DOUBLE PRECISION AA1,DTCMIN,DTDMIN,DTSMIN,
     1     UC,VC,WC,DTC1,DTC2,DTC,ALF1,ALF2,ALF2E,PSUBE,ALF3,
     2     DCOMAX,DCOE1,DTURB,ALF,DTD,DTS,DTCON
      integer indexdt
C==============================================================TIMSTP
      IF(LDT.EQ.0) RETURN
C----------
      DO 100 IC=1,NXYZT
      DTOLD(IC)=DT(IC)
  100 CONTINUE
C----------
      DTCMIN=LARGE
      DTDMIN=LARGE
      DTSMIN=LARGE
C--------------------------------------------------------------TIMSTP
      DO 500 K=K1,K2
      ICK=(K-1)*NXYT
      DO 500 J=J1,J2
      ICJ=ICK+(J-1)*NXT
      DO 500 I=I1,I2
      IC=ICJ+I
C----------------------------------------- BLOCKAGE LOGIC -----TIMSTP
      IF(RC(IC).LE.SMALL) GO TO 500
      AA1=ONE/(RDX(IC)*RDX(IC)+RDY(IC)*RDY(IC)+RDZ(IC)*RDZ(IC))
C--------------------------------------- CONVECTION LIMIT -----TIMSTP
      IF(NX.GT.1) THEN
           UC=HALF*ABS(U(IC)+U(IC-1))
      ELSE
           UC=ZERO
      END IF
      IF(NY.GT.1) THEN
           VC=HALF*ABS(V(IC)+V(IC-NXT))
      ELSE
           VC=ZERO
      END IF
      IF(NZ.GT.1) THEN
           WC=HALF*ABS(W(IC)+W(IC-NXYT))
      ELSE
           WC=ZERO
      END IF
      DTC1=UC*UC*RDX(IC)*RDX(IC)+VC*VC*RDY(IC)*RDY(IC)
     1    +WC*WC*RDZ(IC)*RDZ(IC)
      DTC2=UC*RDX(IC)+VC*RDY(IC)+WC*RDZ(IC)
      DTC=MAX(SMALL,SQRT(DTC1))/MAX(SMALL,DTC2*DTC2)
      DTCMIN=MIN(DTCMIN,DTC)
      DT(IC)=DTC
C---------------------------------------- DIFFUSION LIMIT -----TIMSTP
      ALF1=(TWO+VISRAT)*VISC(IC)/RO(IC)
      IF(ISWIRL.EQ.1) ALF1=EIGHT*VISC(IC)/(THREE*RO(IC))
      IF(NLTE.EQ.0) THEN
        ALF2=COND(IC)*TEMP(IC)*(GAMMA(IC)-ONE)/P(IC)
      ELSE
        PSUBE=SPD(IC,IELC)*RGAS*RMW(IELC)*TE(IC)
        ALF2=COND(IC)*TEMP(IC)*(GAMMA(IC)-ONE)/(P(IC)-PSUBE)
        ALF2E=CONDE(IC)*TE(IC)*TWOO3/MAX(PSUBE,SMALL)
        ALF2=MAX(ALF2,ALF2E)
      END IF
      IF(NODIFF.NE.1) THEN
        DCOMAX=SMALL
        DO 200 ISP=1,NSP
        IF(NLTE.EQ.0) THEN
          DCOMAX=MAX(TWO*DCOEF(IC,NSP),DCOMAX)
        ELSE
          DCOMAX=MAX(DCOEF(IC,NSP)*(ONE+TE(IC)/TEMP(IC)),DCOMAX)
        END IF
  200   CONTINUE
        ALF3=DCOMAX
      ELSE
        ALF3=ZERO
      END IF
      ALF=MAX(ALF1,ALF2,ALF3)
      DTD=HALF*AA1/ALF
      DTDMIN=MIN(DTDMIN,DTD)
      DT(IC)=MIN(DT(IC),DTD)
C-------------------------------------- SOUND SPEED LIMIT -----TIMSTP
      IF(MACCEL.EQ.0) THEN
        IF(NLTE.EQ.0) THEN
          DTS=PGS*SQRT(AA1*RO(IC)/(GAMMA(IC)*P(IC)))
        ELSE
          PSUBE=SPD(IC,IELC)*RGAS*RMW(IELC)*TE(IC)
          DTS=PGS*SQRT(AA1*RO(IC)
     1                 /(GAMMA(IC)*P(IC)+(FIVO3-GAMMA(IC))*PSUBE))
        END IF
        DTSMIN=MIN(DTSMIN,DTS)
        DT(IC)=MIN(DT(IC),DTS)
      END IF
C---------------------------------- SPATIALLY VARIABLE DT -----TIMSTP
      DT(IC)=MIN(DT(IC)*SAFEDT,DTGROW*DTOLD(IC),DTMAX)
  500 CONTINUE
C==============================================================TIMSTP
c ----  change the value of safedt, safdtd to 10 times of its old values
c       to enforce a quick convergence, by Y.P. WAN
      IF(LDT.EQ.1) THEN
      DTCON=MIN(SAFEDT*DTCMIN,SAFDTD*DTDMIN,SAFEDT*DTSMIN,
     1          DTGROW*DTOLD(IC1),DTMAX)
      if(dtcon.eq.SAFEDT*DTCMIN) then
        indexdt=1
      else if(dtcon.eq.SAFDTD*DTDMIN) then
        indexdt=2
      else if(dtcon.eq.SAFEDT*DTSMIN) then
        indexdt=3
      else if(dtcon.eq.DTGROW*DTOLD(IC1)) then
        indexdt=4
      else if(dtcon.eq.DTMAX) then
        indexdt=5
      end if
c by Y.P. WAN -- see the stability limitation 
c      write(*,601)ncyc,indexdt,dtcon,dtcmin,dtdmin,dtsmin,
c     1            dtold(ic1),dtmax
c 601  format(1x,2i5,1x,7e10.3)
      DO 600 IC=1,NXYZT
      DT(IC)=DTCON
  600 CONTINUE
      END IF
C==============================================================TIMSTP
      RETURN
      END
