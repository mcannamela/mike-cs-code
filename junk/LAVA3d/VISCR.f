*DECK VISCR
      SUBROUTINE VISCR
C===============================================================VISCR
C     THIS ROUTINE UPDATES RORU, RORV, AND RORW
C     DUE TO VISCOUS STRESSES DUE TO R FACTOR IN ALL DIRECTION
C----------
C     CALLED BY LAVA
C===============================================================VISCR
      INCLUDE 'COML.h'
C----------
      DOUBLE PRECISION DICYL,AMUR,AMUF,AMUT,UDELRR,UDELRF,UDELRT,
     1     DIVR,DIVF,DIVT,DRCDX,DRCDY,DRCDZ,TMU,URX,DISST,XR
C
C===============================================================VISCR
      DICYL=DBLE(ICYL)
C--------------------------------------------- X-DIRECTION -----VISCR
      IF(NX.NE.1) THEN
        DO 100 IC=IC1U,IC2
        AMUR=(DX(IC+1)*VISC(IC)+DX(IC)*VISC(IC+1))*HRDXCR(IC)
        UDELRR=(DX(IC+1)*UDELR(IC)+DX(IC)*UDELR(IC+1))*HRDXCR(IC)
        DIVR=(DX(IC+1)*DIV(IC)+DX(IC)*DIV(IC+1))*HRDXCR(IC)
        DRCDX=(RC(IC+1)-RC(IC))*TWO*HRDXCR(IC)
        XR=MAX(SMALL,XC(IC)+HALF*DX(IC))
        TMU=AMUR*(TWO*DICYL/XR*(UN(IC)*DRCDX
     1              +RR(IC)*(UDELRR-TWO*UN(IC)/XR))
     2           -(TWO*UDELRR+VISRAT*DIVR)*DRCDX)
        RORU(IC)=RORU(IC)+DT(IC)*TMU
  100   CONTINUE
      END IF
C--------------------------------------------- Y-DIRECTION -----VISCR
      IF(NY.NE.1) THEN
        DO 200 IC=IC1V,IC2
        AMUF=(DY(IC+NXT)*VISC(IC)+DY(IC)*VISC(IC+NXT))*HRDYCF(IC)
        UDELRF=(DY(IC+NXT)*UDELR(IC)+DY(IC)*UDELR(IC+NXT))*HRDYCF(IC)
        DIVF=(DY(IC+NXT)*DIV(IC)+DY(IC)*DIV(IC+NXT))*HRDYCF(IC)
        DRCDY=(RC(IC+NXT)-RC(IC))*TWO*HRDYCF(IC)
        TMU=(DICYL*HRDYCF(IC)*((UN(IC)+UN(IC-1))*DY(IC+NXT)
     1                      +(UN(IC+NXT)+UN(IC+NXT-1))*DY(IC))/XC(IC)
     2       -TWO*UDELRF-VISRAT*DIVF)*AMUF*DRCDY
        RORV(IC)=RORV(IC)+DT(IC)*TMU
  200 CONTINUE
      END IF
C--------------------------------------------- Z-DIRECTION -----VISCR
      IF(NZ.NE.1) THEN
        DO 300 IC=IC1W,IC2
        AMUT=(DZ(IC+NXYT)*VISC(IC)+DZ(IC)*VISC(IC+NXYT))*HRDZCT(IC)
      UDELRT=(DZ(IC+NXYT)*UDELR(IC)+DZ(IC)*UDELR(IC+NXYT))*HRDZCT(IC)
        DIVT=(DZ(IC+NXYT)*DIV(IC)+DZ(IC)*DIV(IC+NXYT))*HRDZCT(IC)
        DRCDZ=(RC(IC+NXYT)-RC(IC))*TWO*HRDZCT(IC)
        RORW(IC)=RORW(IC)-DT(IC)*(TWO*UDELRT+VISRAT*DIVT)*AMUT*DRCDZ
  300   CONTINUE
      END IF
C============================ DISSIPATION DUE TO R FACTORS =====VISCR
      IF(ITURB.EQ.0) THEN
        DO 500 IC=IC1,IC2
        URX=HALF*(UN(IC-1)+UN(IC))/XC(IC)
        DISST=DT(IC)*VISC(IC)*RC(IC)*(DICYL*FOUR*URX*(URX-UDELR(IC))
     1                    +UDELR(IC)*(TWO*UDELR(IC)+VISRAT*DIV(IC)))
        ROER(IC)=ROER(IC)+DISST
  500   CONTINUE
      ELSEIF(ITURB.EQ.1) THEN
        DO 510 IC=IC1,IC2
        URX=HALF*(UN(IC-1)+UN(IC))/XC(IC)
        DISST=DT(IC)*VISC(IC)*RC(IC)*(DICYL*FOUR*URX*(URX-UDELR(IC))
     1                    +UDELR(IC)*(TWO*UDELR(IC)+VISRAT*DIV(IC)))
        ROTKER(IC)=ROTKER(IC)+DISST
  510   CONTINUE
      ELSEIF(ITURB.EQ.2) THEN
        DO 520 IC=IC1,IC2
        URX=HALF*(UN(IC-1)+UN(IC))/XC(IC)
        DISST=DT(IC)*VISC(IC)*RC(IC)*(DICYL*FOUR*URX*(URX-UDELR(IC))
     1                    +UDELR(IC)*(TWO*UDELR(IC)+VISRAT*DIV(IC)))
        ROTKER(IC)=ROTKER(IC)+DISST
        ROEPSR(IC)=ROEPSR(IC)+DISST*CE1*EPSN(IC)/TKEN(IC)
  520   CONTINUE
      END IF
C===============================================================VISCR
      RETURN
      END
