*DECK MOMX
      SUBROUTINE MOMX
C================================================================MOMX
C     THIS ROUTINE UPDATES RORU DUE TO PRESSURE (Q) GRADIENTS,
C     AND CONVERTS BACK FROM RORU TO U.
C----------
C     CALLED BY LAVA
C----------
C     INPUT: DX,Q,RORU,ROR
C     OUTPUT: U
C================================================================MOMX
      INCLUDE 'COML.h'
C----------
      DOUBLE PRECISION DQDX,RORR
C
C================================================================MOMX
      DO 10 IC=IC1U,IC2
      DQDX=TWO*HRDXCR(IC)*(Q(IC+1)-Q(IC)
     1              +TWOO3*(RO(IC+1)*TKE(IC+1)-RO(IC)*TKE(IC)))
      RORR=HRDXCR(IC)*(DX(IC+1)*RO(IC)+DX(IC)*RO(IC+1))*RR(IC)
      RORR=MAX(RORR,SMALL)
      U(IC)=(RORU(IC)-DT(IC)*RR(IC)*DQDX)/RORR
   10 CONTINUE
C================================================================MOMX
      RETURN
      END
