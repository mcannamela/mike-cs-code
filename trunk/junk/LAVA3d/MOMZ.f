*DECK MOMZ
      SUBROUTINE MOMZ
C================================================================MOMZ
C     THIS ROUTINE UPDATES RORW DUE TO PRESSURE (Q) GRADIENTS,
C     AND CONVERTS BACK FROM RORW TO W.
C----------
C     CALLED BY LAVA
C----------
C     INPUT: DZ,Q,RORW,ROR
C     OUTPUT: W
C================================================================MOMZ
      INCLUDE 'COML.h'
C----------
      DOUBLE PRECISION DQDZ,RORT
C
C================================================================MOMZ
      DO 10 IC=IC1W,IC2
      DQDZ=TWO*HRDZCT(IC)*(Q(IC+NXYT)-Q(IC)
     1          +TWOO3*(RO(IC+NXYT)*TKE(IC+NXYT)-RO(IC)*TKE(IC)))
      RORT=HRDZCT(IC)*(DZ(IC+NXYT)*RO(IC)+DZ(IC)*RO(IC+NXYT))*RT(IC)
      W(IC)=(RORW(IC)-DT(IC)*RT(IC)*DQDZ)/RORT
   10 CONTINUE
C================================================================MOMZ
      RETURN
      END
