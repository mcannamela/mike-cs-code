*DECK MOMY
      SUBROUTINE MOMY
C================================================================MOMY
C     THIS ROUTINE UPDATES RORV DUE TO PRESSURE (Q) GRADIENTS,
C     AND CONVERTS BACK FROM RORV TO V.
C----------
C     CALLED BY LAVA
C----------
C     INPUT: DY,Q,RORV,ROR
C     OUTPUT: V
C================================================================MOMY
      INCLUDE 'COML.h'
C----------
      DOUBLE PRECISION DQDY,RORF
C
C================================================================MOMY
      DO 10 IC=IC1V,IC2
      DQDY=TWO*HRDYCF(IC)*(Q(IC+NXT)-Q(IC)
     1            +TWOO3*(RO(IC+NXT)*TKE(IC+NXT)-RO(IC)*TKE(IC)))
      RORF=HRDYCF(IC)*(DY(IC+NXT)*RO(IC)+DY(IC)*RO(IC+NXT))*RF(IC)
      V(IC)=(RORV(IC)-DT(IC)*RF(IC)*DQDY)/RORF
   10 CONTINUE
C================================================================MOMY
      RETURN
      END
