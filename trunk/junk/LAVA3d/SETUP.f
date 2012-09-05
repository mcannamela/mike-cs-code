*DECK SETUP
      SUBROUTINE SETUP
C===============================================================SETUP
C     THIS ROUTINE PROVIDES GRID INFORMATIONS
C     TWO GRID OPTIONS --- UNIFORM GRID AND LINEARLY INCREASING GRID
C     USER SPECIFIED GRID IS PROVIDED BY USERX, USERY, AND USERZ
C----------
C     CALLED BY LAVA
C     CALLS USERX, USERY, USERZ, USERSL
C----------
C     INPUT: NX, NY, NZ, IGRIDX, IGRIDY, IGRIDZ, from DATA BLOCK
C            XL, YL, ZL, X(1), Y(1), Z(1)
C     OUTPUT: I1, I2, J1, J2, K1, K2, NXT, NYT, NZT, NXYT, NXYZT,
C             IC1, IC2, IC1U, IC1V, IC1W,
C             DX, DY, DZ, XC, RDX, RDY, RDZ,
C             INITIAL VALUES, DIVN
C===============================================================SETUP
      INCLUDE 'COML.h'
C----------
      EXTERNAL USERX,USERY,USERZ,USERSL
      INTEGER NOZ,INOZ,JNOZ,KNOZ
      DOUBLE PRECISION DXUNI,DYUNI,DZUNI,DXINC,DYINC,DZINC,DELY,DELZ,
     1     RINJ
C
C===============================================================SETUP
      I1=2
      IF(NX.EQ.1) THEN
           I1=1
           IGRIDX=0
      END IF
      J1=2
      IF(NY.EQ.1) THEN
           J1=1
           IGRIDY=0
      END IF
      K1=2
      IF(NZ.EQ.1) THEN
           K1=1
           IGRIDZ=0
      END IF
C----------
      I2=I1+NX-1
      J2=J1+NY-1
      K2=K1+NZ-1
C----------
C     NXT=2*I2-NX
C     NYT=2*J2-NY
C     NZT=2*K2-NZ
C----------
      NXYT=NXT*NYT
      NXYZT=NXYT*NZT
      IC1=I1+(J1-1)*NXT+(K1-1)*NXYT
      IC2=I2+(J2-1)*NXT+(K2-1)*NXYT
      IC1U=1+(J1-1)*NXT+(K1-1)*NXYT
      IC1V=I1+(K1-1)*NXYT
      IC1W=I1+(J1-1)*NXT
      NIT1=NIT-2
C----------
      NXS=(NXT-NX)/2
      NXL=(NXT+NX)/2
      NYS=(NYT-NY)/2
      NYL=(NYT+NY)/2
      NZS=(NZT-NZ)/2
      NZL=(NZT+NZ)/2
C================ UNIFORM GRID OR LINEARLY INCREASING GRID =====SETUP
      DXUNI=XL/DBLE(NX)
      DYUNI=YL/DBLE(NY)
      DZUNI=ZL/DBLE(NZ)
C--------------------------------------------- X-DIRECTION -----SETUP
C -------- Eqi-distance grid
      IF(IGRIDX.EQ.0) THEN
           X(0)=X(1)-DXUNI
           DO 110 I=1,NXT
           X(I)=X(I-1)+DXUNI
  110      CONTINUE
           IF(NX.EQ.1) THEN
                X(0)=ZERO
                X(1)=ONE
           END IF
      END IF
C---------- USER SPECIFIED GRID
      IF(IGRIDX.EQ.1) THEN
           CALL USERX
      END IF
C--------------------------------------------- Y-DIRECTION -----SETUP
      IF(IGRIDY.EQ.0) THEN
           Y(0)=Y(1)-DYUNI
           DO 210 J=1,NYT
           Y(J)=Y(J-1)+DYUNI
  210      CONTINUE
           IF(NY.EQ.1) THEN
                Y(0)=ZERO
                Y(1)=ONE
           END IF
      END IF
C---------- USER SPECIFIED GRID
      IF(IGRIDY.EQ.1) THEN
           CALL USERY
      END IF
C--------------------------------------------- Z-DIRECTION -----SETUP
      IF(IGRIDZ.EQ.0) THEN
           Z(0)=Z(1)-DZUNI
           DO 310 K=1,NZT
           Z(K)=Z(K-1)+DZUNI
  310      CONTINUE
           IF(NZ.EQ.1) THEN
                Z(0)=ZERO
                Z(1)=ONE
           END IF
      END IF
C---------- USER SPECIFIED GRID
      IF(IGRIDZ.EQ.1) THEN
           CALL USERZ
      END IF
C================================ STORE DX, DY, DZ, AND XC =====SETUP
      DO 400 K=1,NZT
      DELZ=Z(K)-Z(K-1)
      ICK=(K-1)*NXYT
        DO 410 J=1,NYT
        DELY=Y(J)-Y(J-1)
        IC=ICK+(J-1)*NXT+1
          DO 420 I=1,NXT
          DX(IC)=X(I)-X(I-1)
          XC(IC)=X(I-1)+HALF*DX(IC)
          DY(IC)=DELY
          DZ(IC)=DELZ
          IF(NX.EQ.1) THEN
            RDX(IC)=ZERO
          ELSE
            RDX(IC)=ONE/DX(IC)
          END IF
          IF(NY.EQ.1) THEN
            RDY(IC)=ZERO
          ELSE
            RDY(IC)=ONE/DY(IC)
          END IF
          IF(NZ.EQ.1) THEN
            RDZ(IC)=ZERO
          ELSE
            RDZ(IC)=ONE/DZ(IC)
          END IF
          IC=IC+1
  420     CONTINUE
  410   CONTINUE
  400 CONTINUE
C-------------------- VOLUME OF CELLS, INITIALIZE DT ARRAY -----SETUP
C     AREA ARRAY CONTAINES AREA OR VOLUME OF CELLS TO BE USED FOR
C       INTERPOLATION OF VISCOSITIES IN VISCX AND VISCY.
C       AREA ARRAY IS THE SAME AS VOL ARRAY EXCEPT FOR THE CASE OF
C       CYLINDRICAL COORDINATES.
C     VOL IN AXISYMMETRIC CASE IS TOTAL CIRCUMFERENTIAL VOLUME, NOT
C       PER UNIT ANGLE (RADIAN). NOTE THAT ABS(XC) IS USED.
C     HRDXCR, HRDYCF, AND HRDZCT ARRAYS ARE INITIALIZED HERE.
C---------------------------------------------------------------SETUP
      DO 500 IC=1,NXYZT
      AREA(IC)=DX(IC)*DY(IC)*DZ(IC)
      VOL(IC)=AREA(IC)
      IF(ICYL.EQ.1) VOL(IC)=PI2*ABS(XC(IC))*DX(IC)*DY(IC)
      DT(IC)=DTI
      DTOLD(IC)=DTI
  500 CONTINUE
      DO 510 IC=1,IC2
      HRDXCR(IC)=ONE/(DX(IC)+DX(IC+1))
      HRDYCF(IC)=ONE/(DY(IC)+DY(IC+NXT))
      HRDZCT(IC)=ONE/(DZ(IC)+DZ(IC+NXYT))
  510 CONTINUE
C========================== SET SUB-GRID LENGTH SCALE SGSL =====SETUP
      IF(ITURB.EQ.1) THEN
           CALL USERSL
      END IF
C============================= ICNOZ AT LOCATION OF NOZZLE =====SETUP
      IF(IPTCL.EQ.1 .OR. ICARRY.EQ.1) THEN
      DO 790 NOZ=1,NUMNOZ
      IF(NZ.NE.1) THEN
        DO 700 K=2,NZT-1
        IF(ZINJ(NOZ).GE.Z(K-1) .AND. ZINJ(NOZ).LT.Z(K)) KNOZ=K
  700   CONTINUE
      ELSE
        KNOZ=1
      END IF
        DO 710 J=2,NYT-1
        IF(YINJ(NOZ).GE.Y(J-1) .AND. YINJ(NOZ).LT.Y(J)) JNOZ=J
  710   CONTINUE
        RINJ=SQRT(XINJ(NOZ)**2+ZINJ(NOZ)**2)
        DO 720 I=2,NXT-1
        IF(RINJ.GE.X(I-1) .AND. RINJ.LT.X(I)) INOZ=I
  720   CONTINUE
      ICNOZ(NOZ)=INOZ+(JNOZ-1)*NXT+(KNOZ-1)*NXYT
  790 CONTINUE
      END IF
C===============================================================SETUP
      RETURN
      END
