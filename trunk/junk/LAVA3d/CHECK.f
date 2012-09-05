*DECK CHECK
      SUBROUTINE CHECK
C===============================================================CHECK
C     CHECK THERMODYNAMIC VARIABLES. (T AND P)
C----------
C     CALLED BY LAVA
C===============================================================CHECK
      INCLUDE 'COML.h'
C----------
      INTEGER ICS
      EXTERNAL EXITA
C
C===============================================================CHECK
      ICS=0
      DO 10 IC=IC1,IC2
      IF(RC(IC).GT.SMALL) THEN
      IF(TEMP(IC).LT.ONE .OR.P(IC).LE.ZERO .OR.P(IC).GT.LARGE) ICS=IC
      END IF
   10 CONTINUE
      IF(NLTE.NE.0) THEN
        DO 20 IC=IC1,IC2
        IF(RC(IC).GT.SMALL .AND. TE(IC).LT.ZERO) ICS=IC
   20   CONTINUE
      END IF
C---------------------------------------------------------------CHECK
      IF(ICS.NE.0) THEN
      WRITE(NFOUT,900) TIME,NCYC,ICS,TEMP(ICS),TE(ICS),P(ICS),Q(ICS)
      WRITE(*,900)     TIME,NCYC,ICS,TEMP(ICS),TE(ICS),P(ICS),Q(ICS)
        CALL EXITA(2)
      END IF
C===============================================================CHECK
      RETURN
C
  900 FORMAT(' OVERFLOW AT TIME=',1PD12.5,' CYC',I6,' IC=',I6/
     1 ' TEMP=',D12.5,'TE=',D12.5,' P=',D12.5, ' Q=',D12.5)
      END
