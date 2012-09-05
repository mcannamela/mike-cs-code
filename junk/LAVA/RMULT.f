*DECK RMULT
      SUBROUTINE RMULT
C===============================================================RMULT
C     THIS ROUTINE DEFINES COMPOSITE VARIABLES CONTAINING R
C     AND DENSITY
C----------
C     CALLED BY LAVA
C---------------------------------------------------------------RMULT
C     INPUT:SPDN,RON,ROEN,UN,VN,WN,DX,DY,DZ,RC
C     OUTPUT: SPDR,ROER,RORU,RORV,RORW
C===============================================================RMULT
      INCLUDE 'COML.h'
C----------
      DOUBLE PRECISION DXC,DXR,DYC,DYF,DZC,DZT,ROF
C
C=============================== DEFINE (SPECIES RO)*R =========RMULT
      IF(INCOMP.EQ.1) GO TO 100
      DO 20 ISP=1,NSP
           DO 10 IC=1,NXYZT
           SPDR(IC,ISP)=SPDN(IC,ISP)*RC(IC)
   10      CONTINUE
   20 CONTINUE
C============================== DEFINE RO*R*E AND ROE*R*EE =====RMULT
      DO 200 IC=1,NXYZT
      ROER(IC)=ROEN(IC)*RC(IC)
  200 CONTINUE
      IF(NLTE.NE.0) THEN
        DO 300 IC=1,NXYZT
        ROEER(IC)=ROEEN(IC)*RC(IC)
  300   CONTINUE
      END IF
  100 CONTINUE
C================== DEFINE MOMENTUM DENSITY RO*R*VELOCITY ======RMULT
      IF(NX.GT.1) THEN
        DO 440 IC=IC1U,IC2
        DXC=DX(IC)
        DXR=DX(IC+1)
        ROF=(DXR*RON(IC)+DXC*RON(IC+1))/(DXC+DXR)
        RORU(IC)=ROF*RR(IC)*UN(IC)
  440   CONTINUE
      END IF
C----------
      IF(NY.GT.1) THEN
        DO 460 IC=IC1V,IC2
        DYC=DY(IC)
        DYF=DY(IC+NXT)
        ROF=(DYF*RON(IC)+DYC*RON(IC+NXT))/(DYC+DYF)
        RORV(IC)=ROF*RF(IC)*VN(IC)
  460   CONTINUE
      END IF
C----------
      IF(NZ.GT.1) THEN
        DO 480 IC=IC1W,IC2
        DZC=DZ(IC)
        DZT=DZ(IC+NXYT)
        ROF=(DZT*RON(IC)+DZC*RON(IC+NXYT))/(DZC+DZT)
        RORW(IC)=ROF*RT(IC)*WN(IC)
  480   CONTINUE
      END IF
C================================= DEFINE RO*R*W FOR SWIRL =====RMULT
      IF(ISWIRL.EQ.1) THEN
        DO 500 IC=1,NXYZT
        RORW(IC)=RON(IC)*RC(IC)*WN(IC)
  500   CONTINUE
      END IF
C=============================== DEFINE RO*R*TKE, RO*R*EPS =====RMULT
      IF(ITURB.NE.0) THEN
        DO 600 IC=1,NXYZT
        ROTKER(IC)=TKEN(IC)*RON(IC)*RC(IC)
  600   CONTINUE
      END IF
      IF(ITURB.EQ.2) THEN
        DO 620 IC=1,NXYZT
        ROEPSR(IC)=EPSN(IC)*RON(IC)*RC(IC)
  620   CONTINUE
      END IF
C===============================================================RMULT
      RETURN
      END
