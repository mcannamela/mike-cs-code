*DECK LAVA
      PROGRAM LAVA
C================================================================LAVA
C
C     MAIN PROGRAM
C
C----------------------------------------------------------------LAVA
C     WHEN USING BLOCKAGE LOGIC, SET RC = 0, RRC = 1, THEN USING RC
C     BYPASS CHEMISTRY ROUTINES, CHECK, STATE, TIMSTP.
C----------------------------------------------------------------LAVA
      INCLUDE 'COML.h'
      INCLUDE 'COMMPTC.h'
C
      INTEGER NTIMP, myIdx
      DOUBLE PRECISION VOLTMAX,VOLTMIN,TPERIOD,DTNVOL,TNTIM,
     &                 VOLTNEW,TIMNVOL
C================================================================LAVA
      EXTERNAL BCCC,BCMOM,CARRY,CCFLUX,CCFLUY,CCFLUZ,CHECK,
     1     CHEM,CHEMEQ,CHMIMP,
     2     DIFFX,DIFFY,DIFFZ,DUMP,EXITA,INJECT,
     3     MOMFLX,MOMFLY,MOMFLZ,MOMSRC,MOMX,MOMY,MOMZ,OLDTIM,PDVETC,
     4     PFIND,PINT,PLOT,PMOVE,PRINT,PSTATE,RDIVID,REPACK,
     5     RSTART,RINPUT,RMULT,SETR,SETQ,SETUP,START,STATE,
     6     TIMSTP,WALLF,XPORT,VISCR,VISCX,VISCY,VISCZ
C ----- VARIABLES FOR FLUCTUATION CALC. BY Y.P. WAN 1-11-98
C    
      DATA VOLTMAX,VOLTMIN,TPERIOD,DTNVOL /80.D0,40.D0,2.D-4,1.D-5/
C=============================================== DISCLAIMER =====LAVA
c     WRITE(*,1000) 
C================================================================LAVA
      CALL RINPUT
      IF(NCYC.EQ.0) THEN
           CALL SETUP
C  -----  SETUP THE COMPUTATIONAL GEOMETRY FOR PARTICLE HEATING BY Y.P. WAN
C                                                           12-10-97
           IF(ICONDP.EQ.1) THEN
             CALL SETUPPTC
           END IF
           CALL START
           CALL SETR
           CALL XPORT
           CALL BCCC
c -   RLW 7-8-98; Turned off call to PLOT since now using TECPLOT
c          CALL PLOT(0)
           CALL TECPLOT
           CALL PRINT(0)
c           CALL PLOT(0)
      ELSE
           CALL RSTART
           DO 100 IC=1,NXYZT
           DT(IC)=DTOLD(IC)
  100      CONTINUE
           CALL PRINT(0)
      END IF
C ----- FOR POST-PROCESSING DATA, by Y.P. WAN
      IF(IPOST.EQ.1) THEN
        CALL PRINT(11)
        CALL HDFOUT
        STOP
      END IF
C ----- SET THE TIME FOR NEXT VOLTAGE SAMPLING
      TIMNVOL=DTNVOL+TIME
C========================================= MAIN LOOP BEGINS =====LAVA
   10 CONTINUE
   
      NCYC=NCYC+1
      TIME=TIME+DT(IC1)
      CALL OLDTIM
c      DO 11 myIdx=1,NX*NY*NZ
c	     IF(V(myIdx).GT.10.0) THEN
c		    WRITE(*,*) V(myIdx)/10.0
c	   	 END IF
c   11 CONTINUE
C ----- INTRODUCING TORCH FLUCTUATION DUE TO VOLTAGE FLUCTUATION
C                          BY Y.P. WAN  1-11-97
      IF(IFLUCT.EQ.1.AND.TIME.GE.TIMNVOL-HALF*DT(IC1).AND.
     &  TIME.LT.TIMNVOL+HALF*DT(IC1)) THEN
        NTIMP=INT(TIME/TPERIOD)
        TNTIM=NTIMP*TPERIOD
C --    ASSUME SAW TOOTH FLUCTUATION
c        IF(TIME.GE.TNTIM.AND.TIME.LT.(TNTIM+0.5*TPERIOD)) THEN
c          VOLTNEW=(VOLTMAX-VOLTMIN)*2.0*(TIME-TNTIM)/TPERIOD+VOLTMIN
c        ELSE IF(TIME.GE.(TNTIM+0.5*TPERIOD).
c     &          AND.TIME.LT.(TNTIM+TPERIOD)) THEN
c          VOLTNEW=-(VOLTMAX-VOLTMIN)*2.0*(TIME-TPERIOD-TNTIM)/
c     &            TPERIOD+VOLTMIN
c        END IF
C --   ASSUME SINUSOIDAL FLUCTUATION
        VOLTNEW=SIN((TIME-TNTIM)*PI2/TPERIOD)*(VOLTMAX-VOLTMIN)*0.5D0
     &   +(VOLTMIN+VOLTMAX)*0.5D0
C ---  SET THE NEW BOUNDARY CONDITIONS AT SOME SPECIFIED TIME
        WRITE(*,*)'NEW VOLTAGE AT TIME=',TIME,VOLTNEW
        CALL RESETPWR(VOLTNEW)
        TIMNVOL=TIMNVOL+DTNVOL
      END IF
      CALL RMULT
C=========================================== PARTICLE LOGIC =====LAVA
      IF(TIME.GE.TIMPTL .AND. IPTCL.EQ.1) THEN
           IF(NUMNOZ.GE.1) CALL INJECT
C----------------------------------------------------------------LAVA
C     PARTICLE GENERATION LOGIC (BREAKUP, NUCLEATION) AND
C     PARTICLE THERMOPHYSICAL PROPERTIES FOR NEW PARTICLES
C     WILL BE PUT IN LATER, HERE
C----------------------------------------------------------------LAVA
           CALL PMOVE
C----------------------------------------------------------------LAVA
C     FOR PARTICLE STATISTICS MAKE OUTPUT FILE BEFORE REPACK
C----------------------------------------------------------------LAVA
           CALL PRINT(4)
C----------------------------------------------------------------LAVA
           CALL REPACK
C ---- Deleting particles fully evaporated, by Y.P. WAN
           IF(IEVAP.EQ.1) CALL REPACKEVAP
           CALL PFIND
C---------------------------- FOR SIGNIFICANT PARTICLE VOID -----LAVA
C          CALL SETR
C---------------------- PARTICLE HEAT AND MOMENTUM TRANSFER -----LAVA
           CALL PINT
C----------------------------------------------------------------LAVA
C     WHEN THERE IS MASS TRANSFER BETWEEN PARTICLE AND GAS, CHANGES
C     IN PARTICLE THERMOPHYSICAL PROPERTIES, PARTIAL MASS, ETC. NEEDS
C     TO BE REEVALUATED HERE OR IN PINT.
C----------------------------------------------------------------LAVA
C ----   FOR NEW HEATING MODEL, THE TEMP IS CALC. DIRECTLY, NO STATE
C        CORELATION ARE NEEDED. BY Y.P. WAN 12-10-97
           IF(ICONDP.NE.1) THEN
             CALL PSTATE
           END IF
      END IF
c ---- If no coupling between plasma and particle is considered, skip
c      the following proccedure
c      go to 9999
C==================================== FLUXES IN X-DIRECTION =====LAVA
      IF(NX.GT.1) THEN
        CALL CCFLUX
        CALL MOMFLX
        CALL DIFFX
        CALL VISCX
      END IF
C------------------------------------ FLUXES IN Y-DIRECTION -----LAVA
      IF(NY.GT.1) THEN
        CALL CCFLUY
        CALL MOMFLY
        CALL DIFFY
        CALL VISCY
      END IF
C------------------------------------ FLUXES IN Z-DIRECTION -----LAVA
      IF(NZ.GT.1) THEN
        CALL CCFLUZ
        CALL MOMFLZ
        CALL DIFFZ
        CALL VISCZ
      END IF
C---------------------------- VISCOUS TERMS DUE TO R FACTOR -----LAVA
      CALL VISCR
C==================== CARRIER GAS INJECTION THROUGH NOZZLES =====LAVA
c     IF(TIME.GE.TIMPTL .AND. ICARRY.EQ.1) THEN
c          CALL CARRY
c     END IF
C======================================= KINETIC CHEMISTRY ======LAVA
      IF(NRK.GT.0) CALL CHEM
C================== PDV WORK, RADIATION, WALL FUNCTION ETC. =====LAVA
      CALL PDVETC
      IF(ITURB.GT.0) CALL WALLF
C========================== FAST AND EQUILIBRIUM CHEMISTRY ======LAVA
      IF(NRE.GT.0) THEN
        CALL CHMIMP
C----- EQUILIBRIUM CHEMISTRY
CC      CALL CHEMEQ
      END IF
C======================================== SET P, T, Q, ETC. =====LAVA
      CALL RDIVID
      IF(INCOMP.EQ.0) CALL STATE
C----------------------------------------------------------------LAVA
C     EVEN THOUGH BCCC MAY NEED OUTPUT FROM XPORT, WE WILL USE OLD
C     TIME VALUE OF PROPERTIES.
C----------------------------------------------------------------LAVA
      CALL BCCC
      CALL CHECK
      CALL SETQ
C--------------------------------------- PRESSURE GRADIENTS -----LAVA
      IF(NX.GT.1) CALL MOMX
      IF(NY.GT.1) CALL MOMY
      IF(NZ.GT.1) CALL MOMZ
C----------------------------------------- MOMENTUM SOURCES -----LAVA
      CALL MOMSRC
C================================================================LAVA
      CALL BCMOM
      CALL XPORT
c ---- for skipping of the gas phase recalculation
c 9999 continue
      CALL TIMSTP
C================================================================LAVA
      IF(MOD(NCYC,JUMP).EQ.0) THEN
        CALL PRINT(1)
        CALL HDFOUT
      END IF
C----------------------------------------------------------------LAVA
      IF(TIME.GE.TIMOVI-HALF*DT(IC1) .AND.
     1   TIME.LT.TIMOVI+HALF*DT(IC1)) THEN
c
c      RLW removed next statement since this is done in TECPLOT.f
c       TIMOVI=TIMOVI+DTJUMP
        CALL PRINT(2)
c        CALL PLOT(1)
c      Moved tecplot call so plot file is written every DTJUMP interval
         CALL TECPLOT
c 
        IF(IDUMP.EQ.1) THEN
          DO 200 IC=1,NXYZT
          DTOLD(IC)=DT(IC)
  200     CONTINUE
          CALL DUMP
          WRITE(*,*) "RESTART FILE CREATED"
        END IF
      END IF
C================================================================LAVA
      IF(TIME.LT.TMAX.AND.NCYC.LT.NMAX) GO TO 10
      IF(IDUMP.EQ.1) THEN
        DO 300 IC=1,NXYZT
        DTOLD(IC)=DT(IC)
  300   CONTINUE
        CALL DUMP
        WRITE(*,*) "RESTART FILE CREATED"
      END IF
      CALL PRINT(-1)
      CALL PRINT(11)
C ----  output for TecPlot and HDF data format by Y.P. WAN
c     CALL TECPLOT
c      CALL HDFOUT
C================================================================LAVA
      STOP
C================================================================LAVA
 1000 FORMAT (5X,'THIS MATERIAL RESULTED FROM WORK DEVELOPED UNDER'
     1 /5X,'COVERNMENT CONTRACT NO. DE-AC07-94ID13223 AT THE IDAHO'
     2 /5X,'NATIONAL ENGINEERING LABORATORY AND IS SUBJECT TO A BROAD'
     3 /5X,'GOVERNMENT LICENSE.  NEITHER THE UNITED STATES NOR THE'
     4 /5X,'UNITED STATES DEPARTMENT OF ENERGY, NOR ANY OF THEIR'
     5 /5X,'EMPLOYEES MAKE ANY WARRANTY, EXPRESS OR IMPLIED,'
     6 /5X,'INCLUDING THE WARRANTIES OF MERCHANTABILITY OR FITNESS'
     7 /5X,'FOR A PRTICULAR PURPOSE, OR ASSUMES ANY LEGAL LIABILITY'
     8 /5X,'OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR'
     9 /5X,'USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, OR'
     9 /5X,'PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT'
     1 /5X,'INFRINGE PRIVATELY-OWNED RIGHTS.')
C================================================================LAVA
      END

