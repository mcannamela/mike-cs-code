*DECK EXITA
      SUBROUTINE EXITA (IEXIT)
C===============================================================EXITA
C     ABNORMAL STOPS
C----------
      EXTERNAL PRINT,DUMP
      INTEGER IEXIT,IDUMP
C===============================================================EXITA
      IF(IEXIT.EQ.1) THEN
       WRITE(*,*) '--- TERMINATING DUE TO PARAMETER MISMATCH ---'
      ELSEIF(IEXIT.EQ.2) THEN
       WRITE(*,*) '--- LAVA TERMINATING FROM CHECK ---'
c       CALL PRINT(3)
      ELSEIF(IEXIT.EQ.11) THEN
       WRITE(*,*) '--- LAVA TERMINATING FROM INJECT ---'
c       CALL PRINT(3)
      END IF
C===============================================================EXITA
      STOP
      END
