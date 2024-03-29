*DECK VISCX
      SUBROUTINE VISCX
C===============================================================VISCX
C     THIS ROUTINE UPDATES RORU, RORV, AND RORW
C     DUE TO VISCOUS STRESSES IN THE X-DIRECTION
C----------
C     CALLED BY LAVA
C===============================================================VISCX
      INCLUDE 'COML.h'
C----------
      DOUBLE PRECISION AMU,AMU1,DUDX,
     1     DTDX,DTDY,DTDZ,
     2     RRFR,RRFF,RRF2,TRF,HDUDY,HDVDX,
     3     RRTR,RRTT,RRT2,TRT,HDUDZ,HDWDX,
     4     URX,TMU,TMU1,TMU11,TMU2,DISST
C
C==================================== SRXX AND DISSIPATION =====VISCX
      IF(ITURB.EQ.0) THEN
        DO 10 IC=IC1,IC2
        DUDX=(UN(IC)-UN(IC-1))*RDX(IC)
        SRXX(IC)=RC(IC)*VISC(IC)*(TWO*DUDX+VISRAT*DIV(IC))
        ROER(IC)=ROER(IC)+DT(IC)*SRXX(IC)*DUDX
   10   CONTINUE
      ELSEIF(ITURB.EQ.1) THEN
        DO 11 IC=IC1,IC2
        DUDX=(UN(IC)-UN(IC-1))*RDX(IC)
        SRXX(IC)=RC(IC)*VISC(IC)*(TWO*DUDX+VISRAT*DIV(IC))
        ROTKER(IC)=ROTKER(IC)+DT(IC)*SRXX(IC)*DUDX
   11   CONTINUE
      ELSEIF(ITURB.EQ.2) THEN
        DO 12 IC=IC1,IC2
        DUDX=(UN(IC)-UN(IC-1))*RDX(IC)
        SRXX(IC)=RC(IC)*VISC(IC)*(TWO*DUDX+VISRAT*DIV(IC))
        ROTKER(IC)=ROTKER(IC)+DT(IC)*SRXX(IC)*DUDX
      ROEPSR(IC)=ROEPSR(IC)+DT(IC)*SRXX(IC)*DUDX*CE1*EPSN(IC)/TKEN(IC)
   12   CONTINUE
      END IF
C---------------------------------------------------------------VISCX
      DO 20 ICLR=1,NLR
      ICL=1+(ICLR-1)*NXT
      ICR=ICL+NX+1
      SRXX(ICL)=SRXX(ICL+1)
      SRXX(ICR)=SRXX(ICR-1)
   20 CONTINUE
      DO 40 IC=IC1U,IC2
      RORU(IC)=RORU(IC)+DT(IC)*(SRXX(IC+1)-SRXX(IC))*TWO*HRDXCR(IC)
   40 CONTINUE
C===============================================================VISCX
      IF(NY.EQ.1) GO TO 200
      DO 110 IC=IC1U-NXT,IC2
      RRFR=(DY(IC)*RR(IC+NXT)+DY(IC+NXT)*RR(IC))*HRDYCF(IC)
      RRFF=(DX(IC)*RF(IC+1)+DX(IC+1)*RF(IC))*HRDXCR(IC)
      RRF2=RRFR+RRFF
      AMU1=AREA(IC)/VISC(IC)+AREA(IC+NXT+1)/VISC(IC+NXT+1)
     1    +AREA(IC+NXT)/VISC(IC+NXT)+AREA(IC+1)/VISC(IC+1)
      AMU=(AREA(IC)+AREA(IC+1)+AREA(IC+NXT)+AREA(IC+NXT+1))/AMU1
      HDUDY=(UN(IC+NXT)-UN(IC))*HRDYCF(IC)
      HDVDX=(VN(IC+1)-VN(IC))*HRDXCR(IC)
      SRXY(IC)=RRF2*AMU*(HDUDY+HDVDX)
      DISS(IC)=HALF*DT(IC)*SRXY(IC)*(HDUDY+HDVDX)
  110 CONTINUE
      DO 140 IC=IC1V,IC2
      RORV(IC)=RORV(IC)+DT(IC)*(SRXY(IC)-SRXY(IC-1))*RDX(IC)
  140 CONTINUE
C--------------------------------- DISSIPATION DUE TO SRXY -----VISCX
      IF(ITURB.EQ.0) THEN
        DO 150 IC=IC1,IC2
        DISST=DISS(IC)+DISS(IC-1)+DISS(IC-NXT)+DISS(IC-NXT-1)
        ROER(IC)=ROER(IC)+DISST
  150   CONTINUE
      ELSEIF(ITURB.EQ.1) THEN
        DO 160 IC=IC1,IC2
        DISST=DISS(IC)+DISS(IC-1)+DISS(IC-NXT)+DISS(IC-NXT-1)
        ROTKER(IC)=ROTKER(IC)+DISST
  160   CONTINUE
      ELSEIF(ITURB.EQ.2) THEN
        DO 170 IC=IC1,IC2
        DISST=DISS(IC)+DISS(IC-1)+DISS(IC-NXT)+DISS(IC-NXT-1)
        ROTKER(IC)=ROTKER(IC)+DISST
        ROEPSR(IC)=ROEPSR(IC)+DISST*CE1*EPSN(IC)/TKEN(IC)
  170   CONTINUE
      END IF
  200 CONTINUE
C===============================================================VISCX
      IF(NZ.EQ.1) GO TO 300
      DO 210 IC=IC1U-NXYT,IC2
      RRTR=(DZ(IC)*RR(IC+NXYT)+DZ(IC+NXYT)*RR(IC))*HRDZCT(IC)
      RRTT=(DX(IC)*RT(IC+1)+DX(IC+1)*RT(IC))*HRDXCR(IC)
      RRT2=RRTR+RRTT
      AMU1=AREA(IC)/VISC(IC)+AREA(IC+NXYT+1)/VISC(IC+NXYT+1)
     1    +AREA(IC+NXYT)/VISC(IC+NXYT)+AREA(IC+1)/VISC(IC+1)
      AMU=(AREA(IC)+AREA(IC+1)+AREA(IC+NXYT)+AREA(IC+NXYT+1))/AMU1
      HDUDZ=(UN(IC+NXYT)-UN(IC))*HRDZCT(IC)
      HDWDX=(WN(IC+1)-WN(IC))*HRDXCR(IC)
      SRXZ(IC)=RRT2*AMU*(HDUDZ+HDWDX)
      DISS(IC)=HALF*DT(IC)*SRXZ(IC)*(HDUDZ+HDWDX)
  210 CONTINUE
      DO 240 IC=IC1W,IC2
      RORW(IC)=RORW(IC)+DT(IC)*(SRXZ(IC)-SRXZ(IC-1))*RDX(IC)
  240 CONTINUE
C--------------------------------- DISSIPATION DUE TO SRXZ -----VISCX
      IF(ITURB.EQ.0) THEN
        DO 250 IC=IC1,IC2
        DISST=DISS(IC)+DISS(IC-1)+DISS(IC-NXYT)+DISS(IC-NXYT-1)
        ROER(IC)=ROER(IC)+DISST
  250   CONTINUE
      ELSEIF(ITURB.EQ.1) THEN
        DO 260 IC=IC1,IC2
        DISST=DISS(IC)+DISS(IC-1)+DISS(IC-NXYT)+DISS(IC-NXYT-1)
        ROTKER(IC)=ROTKER(IC)+DISST
  260   CONTINUE
      ELSEIF(ITURB.EQ.2) THEN
        DO 270 IC=IC1,IC2
        DISST=DISS(IC)+DISS(IC-1)+DISS(IC-NXYT)+DISS(IC-NXYT-1)
        ROTKER(IC)=ROTKER(IC)+DISST
        ROEPSR(IC)=ROEPSR(IC)+DISST*CE1*EPSN(IC)/TKEN(IC)
  270   CONTINUE
      END IF
  300 CONTINUE
C===============================================================VISCX
      RETURN
      END
