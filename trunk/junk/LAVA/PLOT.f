*DECK PLOT
      SUBROUTINE PLOT(IPLOT)
C================================================================PLOT
C     THIS ROUTINE GENERATES MOVIE FILES FOR LAVA FOR EVERY DTJUMP
C     TIME
C----------
C     IPLOT = 0 :   GENERATE FIXED GRID
C           = 1 :   PRINT FOR MOVIE FILES (EVERY DTJUMP TIME)
C----------
C     CALLED BY LAVA
C================================================================PLOT
      INCLUDE 'COML.h'
C----------
      DOUBLE PRECISION DGOM(NPAR)
      DOUBLE PRECISION UC,VC,WC,DTH,ANGLE,COSTH,SINTH,SWIRL,UX,UZ,
     1     EPDGOM
      INTEGER IPLOT,KREC,NDTH
C----------
      CHARACTER*4 JXXX,JYYY,JZZZ,JUUU,JVVV,JWWW,JTTT,JPPP,PART
      CHARACTER*3 PXX,PYY,PZZ,PUU,PVV,PWW,PTT,PMM,PSS
      CHARACTER*1 PTYPE
      CHARACTER*9 JXNMJ,JYNMJ,JZNMJ,JUNMJ,JVNMJ,JWNMJ,JTNMJ,JPNMJ
      CHARACTER*9 PASCI,PXNMJ,PYNMJ,PZNMJ,PUNMJ,PVNMJ,PWNMJ,PTNMJ,
     1     PMNMJ,PSNMJ
      DATA JXXX,JYYY,JZZZ,JUUU,JVVV,JWWW,JTTT,JPPP
     1    /'jxxx','jyyy','jzzz','juuu','jvvv','jwww','jttt','jppp'/
      DATA PART,PXX,PYY,PZZ,PUU,PVV,PWW,PTT,PMM,PSS
     1 /'part','pxx','pyy','pzz','puu','pvv','pww','ptt','pmm','pss'/
C----------------------------------------------------------------PLOT
C     ANGLE = 90 DEGREE
C----------------------------------------------------------------PLOT
      DATA NDTH /4/
C
C================================================== FORMATS =====PLOT
 9010 FORMAT(A4,I5.5)
 9020 FORMAT(12(1X,E12.5))
 9030 FORMAT(A3,A1,I5.5)
C================================================================PLOT
      DTH=TWO*PI/DBLE(NDTH)
C=============================== INITIAL STAGE - GRID FILES =====PLOT
      IF(IPLOT.EQ.0) THEN
        WRITE(JXNMJ,9010) JXXX,NMJ
        WRITE(JYNMJ,9010) JYYY,NMJ
        WRITE(JZNMJ,9010) JZZZ,NMJ
C----------------------------------------------------------------PLOT
        OPEN(UNIT=20,FILE=JXNMJ,ACCESS='DIRECT',RECL=8)
        OPEN(UNIT=21,FILE=JYNMJ,ACCESS='DIRECT',RECL=8)
        IF(I3DP2D.EQ.1 .OR. NZ.GT.1)
     1    OPEN(UNIT=22,FILE=JZNMJ,ACCESS='DIRECT',RECL=8)
C----------------------------------------------------------------PLOT
      IF(I3DP2D.EQ.0) THEN
        KREC=0
        DO 130 K=K1,K2
          DO 120 J=1,J2
            DO 110 I=1,I2
            KREC=KREC+1
            WRITE(20,REC=KREC) X(I)
            WRITE(21,REC=KREC) Y(J)
            IF(NZ.GT.1) WRITE(22,REC=KREC) HALF*(Z(K)+Z(K-1))
  110       CONTINUE
  120     CONTINUE
  130   CONTINUE
C----------------------------------------------------------------PLOT
C     PSEUDO 3-D PARTICLE DISPLAYED IN CYLINDRICAL COORDINATES
C----------------------------------------------------------------PLOT
      ELSE
        KREC=0
        DO 230 K=1,NDTH
          ANGLE=DTH*DBLE(K-1)
          COSTH=COS(ANGLE)
          SINTH=SIN(ANGLE)
          DO 220 J=1,J2
            DO 210 I=1,I2
            KREC=KREC+1
            WRITE(20,REC=KREC) X(I)*COSTH
            WRITE(21,REC=KREC) Y(J)
            WRITE(22,REC=KREC) X(I)*SINTH
  210       CONTINUE
  220     CONTINUE
  230   CONTINUE
      END IF
C----------------------------------------------------------------PLOT
        CLOSE(UNIT=20)
        CLOSE(UNIT=21)
        IF(I3DP2D.EQ.1 .OR. NZ.GT.1) CLOSE(UNIT=22)
      END IF
C================================ GENERATE FLOW MOVIE FILES =====PLOT
      WRITE(JUNMJ,9010) JUUU,NMJ
      WRITE(JVNMJ,9010) JVVV,NMJ
      WRITE(JWNMJ,9010) JWWW,NMJ
      WRITE(JTNMJ,9010) JTTT,NMJ
      WRITE(JPNMJ,9010) JPPP,NMJ
C----------------------------------------------------------------PLOT
      OPEN(UNIT=23,FILE=JUNMJ,ACCESS='DIRECT',RECL=8)
      OPEN(UNIT=24,FILE=JVNMJ,ACCESS='DIRECT',RECL=8)
      IF(I3DP2D.EQ.1 .OR. NZ.GT.1)
     1  OPEN(UNIT=25,FILE=JWNMJ,ACCESS='DIRECT',RECL=8)
      OPEN(UNIT=26,FILE=JTNMJ,ACCESS='DIRECT',RECL=8)
      OPEN(UNIT=27,FILE=JPNMJ,ACCESS='DIRECT',RECL=8)
C----------------------------------------------------------------PLOT
      IF(I3DP2D.EQ.0) THEN
        KREC=0
        DO 330 K=K1,K2
        ICK=(K-1)*NXYT
          DO 320 J=J1,J2
          ICJ=ICK+(J-1)*NXT
            DO 310 I=I1,I2
            IC=ICJ+I
            KREC=KREC+1
            UC=HALF*(U(IC)+U(IC-1))
            VC=HALF*(V(IC)+V(IC-NXT))
            WRITE(23,REC=KREC) UC
            WRITE(24,REC=KREC) VC
            IF(NZ.GT.1) THEN
              WC=HALF*(W(IC)+W(IC-NXYT))
              WRITE(25,REC=KREC) WC
            END IF
            WRITE(26,REC=KREC) TEMP(IC)
            WRITE(27,REC=KREC) P(IC)
  310       CONTINUE
  320     CONTINUE
  330   CONTINUE
C-------------------------------------- PSEUDO 3-D PARTICLE -----PLOT
      ELSE
        SWIRL=DBLE(ISWIRL)
        KREC=0
        DO 430 K=1,NDTH
          ANGLE=DTH*DBLE(K-1)
          COSTH=COS(ANGLE)
          SINTH=SIN(ANGLE)
          DO 420 J=J1,J2
          ICJ=(J-1)*NXT
            DO 410 I=I1,I2
            IC=ICJ+I
            KREC=KREC+1
            UC=HALF*(U(IC)+U(IC-1))
            VC=HALF*(V(IC)+V(IC-NXT))
            WC=SWIRL*W(IC)
            UX=UC*COSTH-WC*SINTH
            UZ=UC*SINTH+WC*COSTH
            WRITE(23,REC=KREC) UX
            WRITE(24,REC=KREC) VC
            WRITE(25,REC=KREC) UZ
            WRITE(26,REC=KREC) TEMP(IC)
            WRITE(27,REC=KREC) P(IC)
  410       CONTINUE
  420     CONTINUE
  430   CONTINUE
      END IF
C----------------------------------------------------------------PLOT
      CLOSE(UNIT=23)
      CLOSE(UNIT=24)
      IF(I3DP2D.EQ.1 .OR. NZ.GT.1) CLOSE(UNIT=25)
      CLOSE(UNIT=26)
      CLOSE(UNIT=27)
C============================ GENERATE PARTICLE MOVIE FILES =====PLOT
      IF(IPTCL.EQ.1) THEN
        WRITE(PASCI,9010) PART,NMJ
        OPEN(UNIT=30,FILE=PASCI)
        WRITE(30,*) NP
        DO 500 IP=1,NP
        EPDGOM=MAX(EP(IP),EPM(IP))
        EPDGOM=MIN(EPDGOM,EPL(IP))
        DGOM(IP)=(EPDGOM-EPM(IP))/(EPL(IP)-EPM(IP))
        IF(I3DP2D.EQ.0 .AND. NZ.EQ.1) THEN
          WRITE(30,9020) XP(IP),YP(IP),UP(IP),VP(IP),
     1                   TP(IP),DGOM(IP),RADP(IP),PARTN(IP)
        ELSE
          WRITE(30,9020) XP(IP),YP(IP),ZP(IP),
     1         UP(IP),VP(IP),WP(IP),TP(IP),DGOM(IP),RADP(IP),PARTN(IP)
        END IF
        WRITE(30,9020) (PMSP(IP,IPSP),IPSP=1,NPSP)
  500   CONTINUE
        CLOSE(UNIT=30)
C----------------------------------------------------------------PLOT
        DO 530 IPSP=1,NPSP
          PTYPE=CHAR(96+IPSP)
          WRITE(PXNMJ,9030) PXX,PTYPE,NMJ
          WRITE(PYNMJ,9030) PYY,PTYPE,NMJ
          WRITE(PZNMJ,9030) PZZ,PTYPE,NMJ
          WRITE(PUNMJ,9030) PUU,PTYPE,NMJ
          WRITE(PVNMJ,9030) PVV,PTYPE,NMJ
          WRITE(PWNMJ,9030) PWW,PTYPE,NMJ
          WRITE(PTNMJ,9030) PTT,PTYPE,NMJ
          WRITE(PMNMJ,9030) PMM,PTYPE,NMJ
          WRITE(PSNMJ,9030) PSS,PTYPE,NMJ
C----------------------------------------------------------------PLOT
          OPEN(UNIT=31,FILE=PXNMJ,ACCESS='DIRECT',RECL=8)
          OPEN(UNIT=32,FILE=PUNMJ,ACCESS='DIRECT',RECL=8)
          IF(I3DP2D.EQ.1 .OR. NZ.GT.1) THEN
            OPEN(UNIT=33,FILE=PZNMJ,ACCESS='DIRECT',RECL=8)
            OPEN(UNIT=34,FILE=PWNMJ,ACCESS='DIRECT',RECL=8)
          END IF
          OPEN(UNIT=35,FILE=PYNMJ,ACCESS='DIRECT',RECL=8)
          OPEN(UNIT=36,FILE=PVNMJ,ACCESS='DIRECT',RECL=8)
          OPEN(UNIT=37,FILE=PTNMJ,ACCESS='DIRECT',RECL=8)
          OPEN(UNIT=38,FILE=PMNMJ,ACCESS='DIRECT',RECL=8)
          OPEN(UNIT=39,FILE=PSNMJ,ACCESS='DIRECT',RECL=8)
C----------------------------------------------------------------PLOT
C     THIS LOGIC ONLY WORKS FOR PURE SUBSTANCES, NOT MIXTURES.
C----------------------------------------------------------------PLOT
          KREC=0
          DO 520 IP=1,NP
            IF(PMSP(IP,IPSP)/PMASS(IP).GE.HALF) THEN
              KREC=KREC+1
              WRITE(31,REC=KREC) XP(IP)
              WRITE(32,REC=KREC) UP(IP)
              IF(I3DP2D.EQ.1 .OR. NZ.GT.1) THEN
                WRITE(33,REC=KREC) ZP(IP)
                WRITE(34,REC=KREC) WP(IP)
              END IF
              WRITE(35,REC=KREC) YP(IP)
              WRITE(36,REC=KREC) VP(IP)
              WRITE(37,REC=KREC) TP(IP)
              WRITE(38,REC=KREC) DGOM(IP)
              WRITE(39,REC=KREC) RADP(IP)
            END IF
  520     CONTINUE
C----------------------------------------------------------------PLOT
          CLOSE(UNIT=31)
          CLOSE(UNIT=32)
          IF(I3DP2D.EQ.1 .OR. NZ.GT.1) THEN
            CLOSE(UNIT=33)
            CLOSE(UNIT=34)
          END IF
          CLOSE(UNIT=35)
          CLOSE(UNIT=36)
          CLOSE(UNIT=37)
          CLOSE(UNIT=38)
          CLOSE(UNIT=39)
  530   CONTINUE
C----------------------------------------------------------------PLOT
      END IF
C================================================================PLOT
      NMJ=NMJ+1
C================================================================PLOT
      RETURN
      END

