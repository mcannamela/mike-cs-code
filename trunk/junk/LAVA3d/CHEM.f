*DECK CHEM
      SUBROUTINE CHEM
C----------------------------------------------------------------CHEM
C      CALCULATES THE CHANGE IN SPECIES DENSITIES AND INTERNAL ENERGY
C      DUE TO KINETIC CHEMICAL REACTIONS
C----------
C      CALLED BY LAVA
C----------------------------------------------------------------CHEM
       INCLUDE 'COML.h'
C----------
      DOUBLE PRECISION TIC,RTIC,RP,PP,ROM,KB,KF,
     1     TEB,TEF,EKB,EKF,OMEG,RMIN,
     2     FLAM,FLBM,CTOP,CBOT,DECHEM,DECHK,DOMEGA,TCHEM,SS,PHI,
     3     ROMN,SUM,ETA,TOP,BOT,TEIC,RTEIC,DECHME,ROM3RD,
     4     TA,TEA,TALOG,TEALOG,EQCC,EQCCE
      INTEGER IR,IREF,NE,KK,L3RD
C----------
      DIMENSION DOMEGA(LNRK)
C----------------------------------------------------------------CHEM
      TCHEM=TCHEMI
C----------------------------------------------------------------CHEM
C     TRIPLE NESTED DO LOOP IS NECESSAARY TO AVOID BOUNDARY CELLS
C----------------------------------------------------------------CHEM
      DO 100 K=K1,K2
      ICK=(K-1)*NXYT
      DO 100 J=J1,J2
      ICJ=ICK+(J-1)*NXT
      DO 100 I=I1,I2
      IC=ICJ+I
C------------------------------------------- BLOCKAGE LOGIC -----CHEM
      IF(RC(IC).LE.SMALL) GO TO 100
      DO 70 IR=1,NRK
C----------------------------------------------------------------CHEM
C     FOWARD REACTION COEF. AND EQULIIBRIUM CONSTATNTS ARE BASED ON
C     EITHER T OR TE.  OTHERWISE (E.G. TWO-TEMPERATURE SAHA), SPECIAL
C     CODING IS NEEDED.
C----------------------------------------------------------------CHEM
      TIC=(ONE-TFLAGK(IR))*TEMP(IC)+TFLAGK(IR)*TE(IC)
      IF(TIC.LT.TCUTL(IR) .OR. TIC.GT.TCUTH(IR)) GO TO 70
      TA=THOUTH*TIC
      TALOG=LOG(TA)
      RTIC=ONE/TIC
      RP=ONE
      PP=ONE
      NE=NELEM(IR)
      DO 20 KK=1,NE
      ISP=CM(KK,IR)
      ROM=SPDN(IC,ISP)*RMW(ISP)
      IF(AM(ISP,IR).NE.0) THEN
        IF(ROM.LE.ZERO) RP=ZERO
        IF(ROM.GT.ZERO) RP=RP*ROM**AE(ISP,IR)
      END IF
      IF(BM(ISP,IR).NE.0) THEN
        IF(ROM.LE.ZERO) PP=ZERO
        IF(ROM.GT.ZERO) PP=PP*ROM**BE(ISP,IR)
      END IF
   20 CONTINUE
C-------------------------------------- THIRD BODY REACTION -----CHEM
C     EXPONENT FOR THE THIRD BODY IS ASSUMED TO BE ONE.  OTHERWISE,
C     USER SHOULD PROVIDE PROPER FORTRAN CODING.
C     FOR EXAMPLE, RP=RP*ROM3RD SHOULD BE RP=RP*ROM3RD**EXPONENT.
C     (BE CAREFUL ABOUT THE CHAPERONE EFFICIENCY.)
C----------------------------------------------------------------CHEM
      L3RD=0
      ROM3RD=ZERO
      DO 10 ISP=1,NSP
       IF(CHAEFF(ISP,IR).GT.SMALL) THEN
        L3RD=1
        ROM3RD=ROM3RD+CHAEFF(ISP,IR)*SPDN(IC,ISP)*RMW(ISP)
       END IF
   10 CONTINUE
      IF(L3RD.EQ.1) THEN
        RP=RP*ROM3RD
        PP=PP*ROM3RD
      END IF
C----------------------------------------------------------------CHEM
      KB=ZERO
      KF=ZERO
      TEB=ONE
      TEF=ONE
      EKB=ONE
      EKF=ONE
C---------- FORWARD REACTION COEFFICIENT
      IF(EF(IR).NE.ZERO) EKF=EXP(-EF(IR)*RTIC)
      IF(ZETAF(IR).NE.ZERO) TEF=TIC**ZETAF(IR)
      KF=CF(IR)*TEF*EKF
C---------- BACKWARD REACTION COEFFICIENT
      KB=KF/EXP(AKS(IR)*TALOG+BKS(IR)/TA
     1          +CKS(IR)+TA*(DKS(IR)+EKS(IR)*TA))
C----------------------------------------------------------------CHEM
C     IF ANY RATE COEFFICIENTS CANNOT BE PUT IN STANDARD
C     FORM, CODE THEM BY HAND AND PUT THEM HERE
C
C     FIND THE REFERENCE SPECIES (THE ONE IN GREATEST DANGER
C     OF BEING DRIVEN NEGATIVE)
C----------------------------------------------------------------CHEM
      OMEG=KF*RP-KB*PP
      RMIN=ZERO
      IF(OMEG.EQ.ZERO) GO TO 70
      DO 50 KK=1,NE
      ISP=CM(KK,IR)
      IF(SPDR(IC,ISP).LE.ZERO) GO TO 50
      ROM=OMEG*FBMAM(ISP,IR)*MW(ISP)*RC(IC)/SPDR(IC,ISP)
      IF(ROM.LT.RMIN) THEN
        IREF=ISP
        RMIN=ROM
      END IF
   50 CONTINUE
C----------------------------------------------------------------CHEM
      ROM=SPDR(IC,IREF)*RRC(IC)*RMW(IREF)
      ROMN=SPDN(IC,IREF)*RMW(IREF)
      SS=SIGN(ONE,OMEG)
      SUM=KF*RP+KB*PP
      PHI=HALF*(OMEG+SS*SUM)
      ETA=HALF*(OMEG-SS*SUM)
      TOP=PHI*ROM+ETA*ROMN
      BOT=ROMN-DT(IC)*FBMAM(IREF,IR)*PHI
      DOMEGA(IR)=DT(IC)*TOP/BOT
C----------------------------------------------------------------CHEM
      DO 60 ISP=1,NSP
      SPDR(IC,ISP)=SPDR(IC,ISP)+RC(IC)*MW(ISP)*FBMAM(ISP,IR)*DOMEGA(IR)
   60 CONTINUE
      DECHEM=QR(IR)*DOMEGA(IR)*RC(IC)
      DECHK=ABS(DECHEM/ROER(IC))
      ROER(IC)=ROER(IC)+DECHEM*(ONE-FRK(IR))
      IF(NLTE.NE.0) THEN
        DECHME=FBMAM(IELC,IR)*THRHAF*RGAS*TE(IC)*DOMEGA(IR)*RC(IC)
        ROEER(IC)=ROEER(IC)+DECHEM*FEK(IR)+DECHME
      END IF
      TCHEM=MAX(TCHEM,DECHK)
   70 CONTINUE
  100 CONTINUE
C----------------------------------------------------------------CHEM
      RETURN
      END
