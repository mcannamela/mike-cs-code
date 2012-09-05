*DECK VISCZ
      SUBROUTINE VISCZ
C===============================================================VISCZ
C     THIS ROUTINE UPDATES RORU, RORV, AND RORW
C     DUE TO VISCOUS STRESSES IN THE Z-DIRECTION
C----------
C     CALLED BY LAVA
C===============================================================VISCZ
      INCLUDE 'COML.h'
C----------
      DOUBLE PRECISION DWDZ,
     1     AMU,TMU,TMU1,DISST
C
C==================================== SRZZ AND DISSIPATION =====VISCZ
      IF(ITURB.EQ.0) THEN
        DO 10 IC=IC1,IC2
        DWDZ=(WN(IC)-WN(IC-NXYT))*RDZ(IC)
        SRZZ(IC)=RC(IC)*VISC(IC)*(TWO*DWDZ+VISRAT*DIV(IC))
        ROER(IC)=ROER(IC)+DT(IC)*SRZZ(IC)*DWDZ
   10   CONTINUE
      ELSEIF(ITURB.EQ.1) THEN
        DO 11 IC=IC1,IC2
        DWDZ=(WN(IC)-WN(IC-NXYT))*RDZ(IC)
        SRZZ(IC)=RC(IC)*VISC(IC)*(TWO*DWDZ+VISRAT*DIV(IC))
        ROTKER(IC)=ROTKER(IC)+DT(IC)*SRZZ(IC)*DWDZ
   11   CONTINUE
      ELSEIF(ITURB.EQ.2) THEN
        DO 12 IC=IC1,IC2
        DWDZ=(WN(IC)-WN(IC-NXYT))*RDZ(IC)
        SRZZ(IC)=RC(IC)*VISC(IC)*(TWO*DWDZ+VISRAT*DIV(IC))
        ROTKER(IC)=ROTKER(IC)+DT(IC)*SRZZ(IC)*DWDZ
      ROEPSR(IC)=ROEPSR(IC)+DT(IC)*SRZZ(IC)*DWDZ*CE1*EPSN(IC)/TKEN(IC)
   12   CONTINUE
      END IF
C---------------------------------------------------------------VISCZ
      DO 20 ICBT=1,NBT
      ICB=ICBT
      ICT=ICBT+NXYZT-NXYT
      SRZZ(ICB)=SRZZ(ICB+NXYT)
      SRZZ(ICT)=SRZZ(ICT-NXYT)
   20 CONTINUE
      DO 40 IC=IC1W,IC2
      RORW(IC)=RORW(IC)+DT(IC)*(SRZZ(IC+NXYT)-SRZZ(IC))*TWO*HRDZCT(IC)
   40 CONTINUE
C===============================================================VISCZ
      IF(NX.NE.1) THEN
        DO 140 IC=IC1U,IC2
        RORU(IC)=RORU(IC)+DT(IC)*(SRXZ(IC)-SRXZ(IC-NXYT))*RDZ(IC)
  140   CONTINUE
      END IF
C---------------------------------------------------------------VISCZ
      IF(NY.NE.1) THEN
        DO 240 IC=IC1V,IC2
        RORV(IC)=RORV(IC)+DT(IC)*(SRYZ(IC)-SRYZ(IC-NXYT))*RDZ(IC)
  240   CONTINUE
      END IF
C===============================================================VISCZ
      RETURN
      END
