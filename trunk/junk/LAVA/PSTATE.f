*DECK PSTATE
      SUBROUTINE PSTATE
C==============================================================PSTATE
C     EQUATION OF STATE FOR PARTICLE
C     USER HAS TO PROVIDE THIS ROUTINE
C     THIS ROUTINE CONVERTS EP TO TP, AND PROVIDE DTDE
C----------
C     CALLED BY LAVA
C==============================================================PSTATE
      INCLUDE 'COML.h'
C==============================================================PSTATE
C     LINEAR RELATION BETWEEN EP AND TP IS ASSUMED. MORE COMPLEX
C     RELATIONS CAN BE PUT IN AS DESIRED (SUBROUTINE RINPUT SHOULD BE
C     MODIFIED ACCORDINGLY).
C--------------------------------------------------------------PSTATE
      DO 100 IP=1,NP
      IF(EP(IP).LT.EPM(IP)) THEN
        DTPDEP(IP)=RPSPHS(IP)
        TP(IP)=RPSPHS(IP)*EP(IP)
C--------------------------------------------------------------PSTATE
      ELSEIF(EP(IP).GE.EPM(IP) .AND. EP(IP).LE.EPL(IP)) THEN
        DTPDEP(IP)=(TML(IP)-TMM(IP))/(EPL(IP)-EPM(IP))
        TP(IP)=TMM(IP)+DTPDEP(IP)*(EP(IP)-EPM(IP))
C--------------------------------------------------------------PSTATE
      ELSEIF(EP(IP).GT.EPL(IP)) THEN
        DTPDEP(IP)=RPSPHL(IP)
        TP(IP)=RPSPHL(IP)*(EP(IP)-EPL(IP))+TML(IP)
      END IF
  100 CONTINUE
C==============================================================PSTATE
      RETURN
      END
