*DECK RDIVID
      SUBROUTINE RDIVID
C==============================================================RDIVID
C     CALLED BY LAVA
C==============================================================RDIVID
      INCLUDE 'COML.h'
C
C==============================================================RDIVID
      IF(INCOMP.EQ.1) GO TO 490
      DO 200 IC=1,NXYZT
      RO(IC)=ZERO
      ROE(IC)=ROER(IC)*RRC(IC)
  200 CONTINUE
C==============================================================RDIVID
      IF(NLTE.NE.0) THEN
        DO 210 IC=1,NXYZT
        ROEE(IC)=ROEER(IC)*RRC(IC)
  210   CONTINUE
      END IF
C==============================================================RDIVID
      DO 400 ISP=1,NSP
           DO 300 IC=1,NXYZT
           SPD(IC,ISP)=SPDR(IC,ISP)*RRC(IC)
           RO(IC)=RO(IC)+SPD(IC,ISP)
  300      CONTINUE
  400 CONTINUE
  490 CONTINUE
C==============================================================RDIVID
      IF(ISWIRL.EQ.0) GO TO 590
      DO 500 IC=1,NXYZT
      W(IC)=RORW(IC)*RRC(IC)/RO(IC)
  500 CONTINUE
  590 CONTINUE
C==============================================================RDIVID
C     EPS IS NOT ALLOWED TO BE NEGATIVE.
C--------------------------------------------------------------RDIVID
      IF(ITURB.EQ.0) GO TO 690
      DO 600 IC=1,NXYZT
      TKE(IC)=ROTKER(IC)*RRC(IC)/RO(IC)
  600 CONTINUE
      IF(ITURB.EQ.1) GO TO 690
      DO 610 IC=1,NXYZT
      EPS(IC)=ROEPSR(IC)*RRC(IC)/RO(IC)
      IF(EPS(IC).LE.ZERO) THEN
        EPS(IC)=CMU*TKE(IC)*TKE(IC)*RO(IC)/VIST(IC)
      END IF
  610 CONTINUE
  690 CONTINUE
C==============================================================RDIVID
      RETURN
      END
