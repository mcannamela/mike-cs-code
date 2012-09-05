*DECK PMOVE
      SUBROUTINE PMOVE
C===============================================================PMOVE
C     THIS ROUTINE PERFORMS THE TRANSPORT OF SPRAY PARTICLES
C
C
C----------
C     CALLED BY LAVA
C===============================================================PMOVE
      INCLUDE 'COML.h'
C================================= TRANSPORT THE PARTICLES =====PMOVE
      DO 100 IP=1,NP
      XP(IP)=XP(IP)+DT(IC1)*UP(IP)
  100 CONTINUE
C----------
      DO 200 IP=1,NP
      YP(IP)=YP(IP)+DT(IC1)*VP(IP)
  200 CONTINUE
C----------
      IF(I3DP2D.EQ.1 .OR. NZ.GT.1) THEN
        DO 300 IP=1,NP
        ZP(IP)=ZP(IP)+DT(IC1)*WP(IP)
  300   CONTINUE
      END IF
C===============================================================PMOVE
      RETURN
      END
