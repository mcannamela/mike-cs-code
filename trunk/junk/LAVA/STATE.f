*DECK STATE
      SUBROUTINE STATE
C---------------------------------------------------------------STATE
C     CALCULATE STATE VARIABLES: TEMPERATURE, PRESSURE AND GAMMA
C----------
C     CALLED BY LAVA
C---------------------------------------------------------------STATE
       INCLUDE 'COML.h'
C----------
      DOUBLE PRECISION ELO,EHI,FR,CVRO
      DOUBLE PRECISION ROEP,SPDP,TEBEF
C
C---------------------------------------------------------------STATE
C     TRIPLE NESTED DO LOOP IS MORE EFFICIENT FOR TABLE LOOK UP
C---------------------------------------------------------------STATE
      DO 70 K=K1,K2
      ICK=(K-1)*NXYT
      DO 70 J=J1,J2
      ICJ=ICK+(J-1)*NXT
      DO 70 I=I1,I2
      IC=ICJ+I
C---------------------------------------------------------------STATE
C     ROEP = HEAVY PARTICLE ENERGY FOR NLTE>=1
C     ROEP = TOTAL ENERGY FOR NLTE=0
C---------------------------------------------------------------STATE
      IF(NLTE.EQ.0) THEN
        ROEP=MAX(ROE(IC),SMALL)
      ELSE
C----------------------------- OBTAIN ELECTRON TEMPERATURE -----STATE
        TEBEF=TE(IC)
        TE(IC)=TWOO3*MW(IELC)*ROEE(IC)/(SPD(IC,IELC)*RGAS+SMALL)
        ROEP=MAX(ROE(IC)-ROEE(IC),SMALL)
        IF(SPD(IC,IELC).LE.SMALL) TE(IC)=TEBEF
      END IF
C----------------------- OBTAIN HEAVY PARTICLE TEMPERATURE -----STATE
      IT=INT(HUNDTH*TEMP(IC))
      IT=MIN(IT,NIT1)
   10 ELO=ZERO
      EHI=ZERO
      DO 20 ISP=1,NSP
      SPDP=MAX(SPD(IC,ISP),SMALL*MW(ISP))
      ELO=ELO+SPDP*EK(IT+1,ISP)
      EHI=EHI+SPDP*EK(IT+2,ISP)
   20 CONTINUE
      IF(NLTE.NE.0) THEN
       ELO=ELO-MAX(SPD(IC,IELC),SMALL*MW(IELC))*EK(IT+1,IELC)
       EHI=EHI-MAX(SPD(IC,IELC),SMALL*MW(IELC))*EK(IT+2,IELC)
      END IF
      IF(ROEP.LE.EHI) GO TO 30
      IF(IT.EQ.NIT1) GO TO 40
      IT=IT+1
      GO TO 10
   30 IF(ROEP.GE.ELO) GO TO 40
      IF(IT.EQ.0) GO TO 40
      IT=IT-1
      GO TO 10
   40 FR=(ROEP-ELO)/(EHI-ELO)
      TEMP(IC)=(DBLE(IT)+FR)*HUNDRD
C------------------------------- OBTAIN PRESSURE AND GAMMA -----STATE
      CVRO=(EHI-ELO)*HUNDTH
      P(IC)=ZERO
      DO 50 ISP=1,NSP
      SPDP=MAX(SPD(IC,ISP),SMALL*MW(ISP))
      P(IC)=P(IC)+SPDP*RMW(ISP)*DBLE(IGAS(ISP))
   50 CONTINUE
      IF(NLTE.NE.0) THEN
       P(IC)=P(IC)-MAX(SPD(IC,IELC),SMALL*MW(IELC))*RMW(IELC)
      END IF
C---------------------------------------------------------------STATE
C     P IS THE TOTAL PRESSURE.
C     GAMMA IS THE SPECIFIC HEAT RATIO FOR HEAVY PARTICLES FOR NLTE>0
C---------------------------------------------------------------STATE
      P(IC)=P(IC)*RGAS*TEMP(IC)
      GAMMA(IC)=ONE+P(IC)/(CVRO*TEMP(IC))
C---------------------------------------------------------------STATE
C     ADD P AND PSUBE FROM UPDATED ELECTRON QUANTITIES.
C     SET TE = T FOR NLTE = 0
C---------------------------------------------------------------STATE
      IF(NLTE.EQ.0) THEN
        TE(IC)=TEMP(IC)
      ELSE
        P(IC)=P(IC)+SPD(IC,IELC)*RGAS*RMW(IELC)*TE(IC)
      END IF
   70 CONTINUE
C---------------------------------------------------------------STATE
      RETURN
      END
