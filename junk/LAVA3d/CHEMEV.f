*DECK CHEMEV
      SUBROUTINE CHEMEV(TIC,SPD,NSP,LNRE,NRE,
     1                AS,BS,CS,DS,ES,AN,BN,CN,NLM,FBNAN,MW,RMW)
C==============================================================CHEMEV
C     THIS SUBROUTINE CALCULATES THE CHANGE IN SPECIES DENSITIES AND
C     INTERNAL ENERGY DUE TO IMPLICIT FAST KINETIC AND EQUILIBRIUM
C     CHEMICAL REACTIONS, FOR A GENERAL REACTION SET.
C----------
C     CHEMEV CALLS THE FOLLOWING LAPACK SUBROUTINES: DGESVX AND
C     DGELSS WHICH ARE DESCRIBED BELOW WHERE THE CALL IS MADE.
C==============================================================CHEMEV
      IMPLICIT NONE
      INTEGER NSP,LNRE,NRE,NSP1,LNRE1
      PARAMETER (NSP1=15,LNRE1=3)
      INTEGER IRE,NE,KK,KREF,AREF,BREF,M,MA,IREA,NREA,ISKIP(LNRE1),
     1     IPCON(LNRE1),ICONV,N,NDIAG,NCHEM,
     2     IWORK(LNRE1),IPIV(LNRE1),INFO,IRANK
      CHARACTER*1 EQUED
      EXTERNAL DGESVX,DGELSS
C----------
      DOUBLE PRECISION TIC,RTIC,TA,TALOG,EQC(LNRE1),
     1     EQCC,REQCC,
     2     BAM(LNRE1,LNRE1),RHS(LNRE1),FS,RDIAG,
     3     PROG(LNRE1),DPROG(LNRE1),PRMAX,PRMIN,REF,RP,RPP,PP,
     4     TEST,ROM,AMBROM,PWR,DIR,DIP,PPP,VS,VSPM,DVDOSS,DVDOST,
     5     ASS,BSS,AST,BST,US,USPM,
     6     DUDOSS,DUDOST,FACT,DSPD(NSP1),DROE,FACTOR,AVG,SCALE,
     7     WORK(5*LNRE1),SWORK(LNRE1),RCOND,FERR(1),BERR(1),
     8     RWORK(LNRE1),CWORK(LNRE1),BAMF(LNRE1,LNRE1),MCHACC
C--------------------------------------------------------------CHEMEV
      DOUBLE PRECISION AS(LNRE1),BS(LNRE1),CS(LNRE1),
     1     DS(LNRE1),ES(LNRE1),
     1     FBNAN(NSP1,LNRE1),MW(NSP1),RMW(NSP1)
     2    ,SPD(NSP1)
      INTEGER NLM(LNRE1),AN(NSP1,LNRE1),
     1        BN(NSP1,LNRE1),CN(NSP1,LNRE1)
      INTEGER ISP
C--------------------------------------------------------------CHEMEV
      DOUBLE PRECISION EPSCHM,OMGCHM,OMGCM2
      DOUBLE PRECISION ZERO,ONE,HUNDTH,THOUTH,SMALL,LARGE,PNTNIN
C--------------------------------------------------------------CHEMEV
      DATA ZERO,ONE,HUNDTH,THOUTH /0.0D0,1.0D0,1.0D-2,1.0D-3/
      DATA PNTNIN /0.9D0/
      DATA SMALL,LARGE /1.0D-80,1.0D+100/
      DATA EPSCHM,OMGCHM,OMGCM2 /0.02D0,1.0D0,0.98D0/
C--------------------------------------------------------------CHEMEV
      DATA NCHEM,NDIAG /50,10/
C---------- MCHACC=100*MACHINE ACCURACY
      DATA MCHACC/1.0D-13/
C==============================================================CHEMEV
c  ---- check if the dimension is correct,  by Y.P. WAN
      IF(NSP.NE.NSP1.OR.LNRE.NE.LNRE1) THEN
        WRITE(*,*)'****** ERROR DIMENSION NSP,LNRE=',
     1             NSP,NSP1,LNRE,LNRE1
        STOP
      END IF
      RTIC=ONE/TIC
      TA=THOUTH*TIC
      TALOG=LOG(TA)
C--------------------------------------------------------------CHEMEV
      NREA=0
      DO 10 IRE=1,NRE
          ISKIP(IRE)=0
          NREA=NREA+1
          EQC(IRE)=EXP(AS(IRE)*TALOG+BS(IRE)/TA
     &                 +CS(IRE)+TA*(DS(IRE)+ES(IRE)*TA))
C--------------------------------------------------------------CHEMEV
          PROG(NREA)=ZERO
          DPROG(NREA)=ZERO
          IPCON(IRE)=1
   10 CONTINUE
C--------------------------------------------------------------CHEMEV
      N=0
      FACTOR=ONE
C-------------------------------------- THE ITERATION LOOP ----CHEMEV
   20 N=N+1
      ICONV=0
      IREA=0
      DO 100 IRE=1,NRE
       IF(ISKIP(IRE).EQ.0) THEN
        IREA=IREA+1
        PROG(IREA)=PROG(IREA)+FACTOR*DPROG(IREA)
        EQCC=EQC(IRE)
        REQCC=ONE/EQCC
C--------------------------------------------------------------CHEMEV
        PRMAX=ZERO
        PRMIN=ZERO
        REF=ZERO
        RP=ONE
        PP=ONE
        ASS=ZERO
        BSS=ZERO
        NE=NLM(IRE)
C--------- FIND REFERENCE SPECIES AND CALCULATE PP AND RP -----CHEMEV
        DO 40 KK=1,NE
         ISP=CN(KK,IRE)
         TEST=SMALL*MW(ISP)
         IF(N.EQ.1 .AND. SPD(ISP).LE.TEST) SPD(ISP)=TEST
         ROM=SPD(ISP)*RMW(ISP)
         AMBROM=-FBNAN(ISP,IRE)/ROM
         TEST=ABS(AMBROM)
         IF(TEST.GT.REF) KREF=ISP
         REF=MAX(REF,TEST)
         PRMAX=MAX(PRMAX,AMBROM)
         PRMIN=MIN(PRMIN,AMBROM)
         IF(AN(ISP,IRE).GT.0) RP=RP*ROM**AN(ISP,IRE)
         IF(BN(ISP,IRE).GT.0) PP=PP*ROM**BN(ISP,IRE)
         FACT=FBNAN(ISP,IRE)*MW(ISP)/SPD(ISP)
         ASS=ASS+AN(ISP,IRE)*FACT
         BSS=BSS+BN(ISP,IRE)*FACT
   40   CONTINUE
        PRMAX=ONE/PRMAX
        PRMIN=ONE/PRMIN
        AREF=AN(KREF,IRE)
        BREF=BN(KREF,IRE)
C--------------------------------------------------------------CHEMEV
        PWR=ONE
        SCALE=MAX(RP,REQCC*PP)
        TEST=REQCC*PP-RP
        IF(ABS(TEST).GT.EPSCHM*SCALE) ICONV=1
C--------------------------------------------------------------CHEMEV
        DIR=DBLE(AREF*(AREF-1))*RP
        DIP=DBLE(BREF*(BREF-1))*REQCC*PP
        IF(DIP.GT.DIR) THEN
         US=EQCC*RP
         PPP=PP
         USPM=ONE
         IF(IPCON(IRE).EQ.1 .AND. BREF.GT.1) THEN
          PWR=ONE/DBLE(BREF)
          PPP=PP**PWR
          USPM=US**(PWR-ONE)
         END IF
         FS=PPP-USPM*US
         DUDOSS=EQCC*RP*ASS
         RDIAG=ONE/(PWR*(PPP*BSS-USPM*DUDOSS))
        ELSE
         VS=REQCC*PP
         RPP=RP
         VSPM=ONE
         IF(IPCON(IRE).EQ.1 .AND. AREF.GT.1) THEN
          PWR=ONE/DBLE(AREF)
          RPP=RP**PWR
          VSPM=VS**(PWR-ONE)
         END IF
         FS=RPP-VSPM*VS
         DVDOSS=REQCC*PP*BSS
         RDIAG=ONE/(PWR*(RPP*ASS-VSPM*DVDOSS))
        END IF
C--------------------------------------------------------------CHEMEV
        IF(N.LE.NDIAG) THEN
C----------------- ONE-STEP GAUSS-SEIDEL-NEWTON ITERATION -----CHEMEV
         DPROG(IREA)=-OMGCHM*RDIAG*FS
         DPROG(IREA)=MIN(DPROG(IREA),PNTNIN*PRMAX)
         DPROG(IREA)=MAX(DPROG(IREA),PNTNIN*PRMIN)
C--------------------------------------------------------------CHEMEV
CDIR$ IVDEP
         DO 60 KK=1,NE
          ISP=CN(KK,IRE)
   60     SPD(ISP)=SPD(ISP)+MW(ISP)*FBNAN(ISP,IRE)*DPROG(IREA)
        ELSE
C----------------------------------- FOR NEWTON ITERATION -----CHEMEV
         BAM(IREA,IREA)=ONE
         RHS(IREA)=-RDIAG*FS
         MA=0
         DO 90 M=1,NRE
         IF(ISKIP(M).EQ.0) THEN
          MA=MA+1
          IF(M.NE.IRE) THEN
           AST=ZERO
           BST=ZERO
           DO 80 KK=1,NE
            ISP=CN(KK,IRE)
            FACT=FBNAN(ISP,M)*MW(ISP)/SPD(ISP)
            AST=AST+AN(ISP,IRE)*FACT
   80       BST=BST+BN(ISP,IRE)*FACT
C--------------------------------------------------------------CHEMEV
           IF(DIP.GT.DIR) THEN
            DUDOST=EQCC*RP*AST
            BAM(IREA,MA)=PWR*(PPP*BST-USPM*DUDOST)*RDIAG
           ELSE
            DVDOST=REQCC*PP*BST
            BAM(IREA,MA)=PWR*(RPP*AST-VSPM*DVDOST)*RDIAG
           END IF
C--------------------------------------------------------------CHEMEV
          END IF
         END IF
   90    CONTINUE
        END IF
      END IF
  100 CONTINUE
C==============================================================CHEMEV
      IF(N.LE.NDIAG) GO TO 180
C-------------------------------------------- LAPACK CALL -----CHEMEV
C     LAPACK SUBROUTINES NEEDED ARE;
C     DGESVX, DGERFS, DGETRS, DGECON, DLATRS, DLACON, DGETRF, DLASWP,
C     DGETF2, DLAQGE, DGEEQU, DGELSS, DGELQF, DGELQ2, DLACPY, DRSCL,
C     DLABAD, DBDSQR, DLAS2, DLASV2, DLASR, DLARTG, DORGBR, DORGLQ,
C     DORGL2, DORGQR, DORG2R, DORMBR, DORMLQ, DORML2, DORMQR, DORM2R,
C     DGEBRD, DGEBD2, DLABRD, DGEQRF, DLARFB, DLARFT, DGEQR2, DLARF,
C     DLARFG, DLASET, DLASCL, DLAMC1, DLAMC2, DLAMC4, DLAMC5, DLASSQ,
C     DSYMM, DSYRK, DSYR2K, DGEMM, DGBMV, DSYMV, DSBMV, DSPMV, DTBMV,
C     DTPMV, DTBSV, DTPSV, DSYR, DSPR, DSYR2, DSPR2, DGEMV, DGER, DTRMM,
C     DTRMV, DTRSM, XERBLA, DTRSV, DCOPY, DSCAL, DSWAP, DROT, DAXPY,
C--------------------------------------------------------------CHEMEV
C     SUBROUTINE DGESVX SOLVES THE LINEAR SYSTEM A*X=B, USING LU
C     DECOMPOSITION.  HERE, A IS THE MATRIX BAM(L,M), X IS THE
C     VECTOR DPROG(L), AND B IS THE VECTOR RHS(L).
C-------------------------------------------- LAPACK CALL -----CHEMEV
C     SUBROUTINE DGELSS (LAPACK) SOLVES THE SAME SYSTEM A*X=B BY 
C     SINGULAR VALUE DECOMPOSITION, WHEN A IS SINGULAR 
C     (UNDERDETERMINED SYSTEM).  THIS ROUTINE IS NECESSARY FOR A
C     LINEARLY DEPENDENT SYSTEM OF KINETIC CHEMCIAL REACTIONS IN
C     EQUILIBRIUM LIMIT.  THIS ROUTINE RETURNS X IN THE SAME STORAGE
C     USED FOR B.
C--------------------------------------------------------------CHEMEV
      CALL DGESVX('N','N',NREA,1,BAM,LNRE,BAMF,LNRE,IPIV,EQUED,
     1            RWORK,CWORK,RHS,LNRE,DPROG,LNRE,RCOND,FERR,BERR,
     2            WORK,IWORK,INFO)
      IF(RCOND.LE.MCHACC) THEN
C--------------------------------------------------------------CHEMEV
C     THIS WRITE STATEMENT SHOULD BE DISABLED, IF SINGULAR MATRIX IS
C     EXPECTED -- LINEARLY DEPENDENT EQUILIBRIUM REACTIONS
C--------------------------------------------------------------CHEMEV
        WRITE(*,910) N
        CAll DGELSS(NREA,NREA,1,BAM,LNRE,RHS,LNRE,SWORK,MCHACC,IRANK,
     1              WORK,5*LNRE,INFO)
        DO 110 IREA=1,NREA
         DPROG(IREA)=RHS(IREA)
 110    CONTINUE
        IF(INFO.NE.0) THEN
          WRITE(*,*) '-- INFO=',INFO
          GO TO 200
        END IF
      END IF
C--------------------------------------------------------------CHEMEV
      DO 120 ISP=1,NSP
  120 DSPD(ISP)=ZERO
      IREA=0
      DO 140 IRE=1,NRE
       IF(ISKIP(IRE).EQ.0) THEN
        IREA=IREA+1 
        NE=NLM(IRE)
CDIR$ IVDEP
        DO 130 KK=1,NE
         ISP=CN(KK,IRE)
         DSPD(ISP)=DSPD(ISP)+MW(ISP)*FBNAN(ISP,IRE)*DPROG(IREA)
  130   CONTINUE
       END IF
  140 CONTINUE
      TEST=ZERO
      DO 150 ISP=1,NSP
      TEST=MIN(TEST,DSPD(ISP)/SPD(ISP))
  150 CONTINUE
      FACTOR=MIN(ONE,-OMGCM2/MIN(TEST,-SMALL))
      DO 170 ISP=1,NSP
  170 SPD(ISP)=SPD(ISP)+FACTOR*DSPD(ISP)
C--------------------------------------------------------------CHEMEV
  180 IF(ICONV.EQ.0) GO TO 190
      IF(N.LT.NCHEM) GO TO 20
      WRITE(*,900) N
  190 CONTINUE
  200 CONTINUE
C==============================================================CHEMEV
      RETURN
C
  900 FORMAT(' WARNING - NOT CONVERGED IN CHEMEV, ITS=',I3)
  910 FORMAT(' SINGULAR MATRIX AT ITER',I3)
      END
