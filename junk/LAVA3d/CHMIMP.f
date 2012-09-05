*DECK CHMIMP
      SUBROUTINE CHMIMP
C==============================================================CHMIMP
c     version 14 -- Please, God! -- 8/19/93
c     version 4
c     Two-Temperature version 4/27/95
C==============================================================CHMIMP
C     This routine cannot handle reactions with
C     Ke (T,Te) = kf(Te)/kb(T)  (e.g. ArH+ + e- <--> H* + Ar )
C     Further generalization is needed
C==============================================================CHMIMP
C     THIS SUBROUTINE CALCULATES THE CHANGE IN SPECIES DENSITIES AND
C     INTERNAL ENERGY DUE TO IMPLICIT FAST KINETIC AND EQUILIBRIUM
C     CHEMICAL REACTIONS, FOR A GENERAL REACTION SET.
C----------
C     CALLED BY LAVA
C----------
C     CHMIMP CALLS THE FOLLOWING LAPACK SUBROUTINES: DGESVX AND
C     DGELSS WHICH ARE DESCRIBED BELOW WHERE THE CALL IS MADE.
C==============================================================CHMIMP
      INCLUDE 'COML.h'
C----------
      DOUBLE PRECISION TIC,RTIC,RDTIC,ROCV,TA,TALOG,DELTAT,EQC(LNRE),
     1     EQCC,REQCC,DLNKDT(LNRE),KFE(LNRE),EKFE,TEFE,DLKFDT(LNRE),
     2     KFEC,RKFEC,BAM(LNRE,LNRE),RHS(LNRE),FS,RDIAG,
     3     PROG(LNRE),DPROG(LNRE),PRMAX,PRMIN,REF,RP,RPP,RPC,PP,PPC,
     4     TEST,ROM,AMB,AMBROM,PWR,DIR,DIP,PPP,VS,VSPM,DVDOSS,DVDOST,
     5     ASS,BSS,AST,BST,ALFSS,BETASS,ALFST,BETAST,PROGT,US,USPM,
     6     DUDOSS,DUDOST,FACT,DSPDR(NSP),DROER,FACTOR,AVG,SCALE,
     7     PROGTS,WORK(5*LNRE),SWORK(LNRE),RCOND,FERR(1),BERR(1),
     8     RWORK(LNRE),CWORK(LNRE),BAMF(LNRE,LNRE),MCHACC,CONC,ERAT,
     9     ROCVE,ROCVH,DTDWS,DTEDWS,QEFFS,DTDWT,DTEDWT,QEFFT,DROEER
      INTEGER IRE,NE,KK,KREF,AREF,BREF,M,MA,IREA,NREA,ISKIP(LNRE),
     1     IPCON(LNRE),ICONV,N,NDIAG,NZONES,NSUM,NCELLS,NCHEM,ICMAX,
     2     IWORK(LNRE),IPIV(LNRE),INFO,IRANK,ICONR,ISUMA,ISUMB
      CHARACTER*1 EQUED
      EXTERNAL DGESVX,DGELSS
C--------------------------------------------------------------CHMIMP
      DATA NCHEM,NDIAG /50,10/
C---------- MCHACC=100*MACHINE ACCURACY
      DATA MCHACC/1.0D-13/
C==============================================================CHMIMP
      NZONES=0
      NSUM=0
      NCELLS=0
C--------------------------------------------------------------CHMIMP
C     TRIPLE NESTED DO LOOP IS NECESSARY TO AVOID BOUNDARY CELLS
C--------------------------------------------------------------CHMIMP
      DO 200 K=K1,K2
      ICK=(K-1)*NXYT
      DO 200 J=J1,J2
      ICJ=ICK+(J-1)*NXT
      DO 200 I=I1,I2
      IC=ICJ+I
C----------------------------------------- BLOCKAGE LOGIC -----CHMIMP
      IF(RC(IC).LE.SMALL) GO TO 200
C--------------------------------------------------------------CHMIMP
      RDTIC=ONE/DT(IC)
      ROCV=ZERO
      CONC=ZERO
      IT=INT(HUNDTH*TEMP(IC))+1
      DO 5 ISP=1,NSP
       CONC=CONC+SPDN(IC,ISP)*RMW(ISP)
       ROCV=ROCV+SPDN(IC,ISP)*HUNDTH*(EK(IT+1,ISP)-EK(IT,ISP))
    5 CONTINUE
      IF(IELC.NE.0) THEN
       ROCVE=SPDN(IC,IELC)*RMW(IELC)*THRHAF*RGAS
       ROCVH=ROCV-ROCVE
      END IF
      NZONES=NZONES+1
C--------------------------------------------------------------CHMIMP
       NREA=0
      DO 10 IRE=1,NRE
        TIC=(ONE-TFLAGE(IRE))*TEMP(IC)+TFLAGE(IRE)*TE(IC)
        RTIC=ONE/TIC
        TA=THOUTH*TIC
        TALOG=LOG(TA)
C--------------------------------------------------------------CHMIMP
        IF(TIC.GT.TCUTEL(IRE) .AND. TIC.LT.TCUTEH(IRE)) THEN
          ISKIP(IRE)=0
          NREA=NREA+1
          EQC(IRE)=EXP(AS(IRE)*TALOG+BS(IRE)/TA
     &                 +CS(IRE)+TA*(DS(IRE)+ES(IRE)*TA))
          DLNKDT(IRE)=THOUTH
     1             *(AS(IRE)/TA-BS(IRE)/TA**2+DS(IRE)+TWO*ES(IRE)*TA)
C--------------------------------------------------------------CHMIMP
C     IF ANY EQUILIBRIUM CONSTANT CANNOT BE PUT IN STANDARD
C     FORM, CODE EQC AND DLNKDT BY HAND AND PUT THEM HERE
C--------------------------------------------------------------CHMIMP
c------ charge exchange and dissociative recombination of nitrogen
      if(ire.eq.5) then
        eqc(ire)=exp(1.0976d0*talog-6.9610d0/ta
     1               -2.6627d0+0.03124d0*ta-2.8837d0/ta**2)
        dlnkdt(ire)=thouth*(1.0976d0/ta+6.9610d0/ta**2+0.03124d0
     1                      +two*2.8837d0/ta**3)
      endif
      if(ire.eq.6) then
        eqc(ire)=exp(1.552891d0*talog-62.453d0/ta-12.96607d0
     1               +0.475054d0*ta-0.0269699d0*ta**2-2.8837d0/ta**2)
        dlnkdt(ire)=thouth*(1.552891d0/ta+62.453d0/ta**2+0.475054d0
     1                      -two*0.0269699d0*ta+two*2.8837d0/ta**3)
      endif
c------
C------- FORWARD REACTION COEFFICIENTS FOR FAST REACTIONS -----CHMIMP
          EKFE=ONE
          TEFE=ONE
          IF(EFE(IRE).NE.ZERO) EKFE=EXP(-EFE(IRE)*RTIC)
          IF(ZETAFE(IRE).NE.ZERO) TEFE=TIC**ZETAFE(IRE)
          KFE(IRE)=CFE(IRE)*TEFE*EKFE
          DLKFDT(IRE)=(ZETAFE(IRE)+EFE(IRE)*RTIC)*RTIC
C--------------------------------------------------------------CHMIMP
C     IF ANY REACTION RATE CANNOT BE PUT IN STANDARD
C     FORM, CODE KFE AND DLKFDT BY HAND AND PUT THEM HERE
C--------------------------------------------------------------CHMIMP
c----- use Hoffert-Lien model for reaction rates and DLKFDT
      if(ire.eq.1) then
          kfe(ire)=2.2584d8*tic*sqrt(tic)*(1.353d5*rtic+two)
     1             *exp(-1.353d5*rtic)
        dlkfdt(ire)=rtic*(1.5d0+1.353d5*(rtic-one/(1.353d5+two*tic)))
          if(kfe(ire).le.small) iskip(ire)=1
      endif
C==============================================================CHMIMP
C     IF ROCVE IS SMALL, THEN SKIP REACTIONS INVOLVING ELECTRONS
C     THIS PART IS STILL UNDER CONSTRUCTION, AND DO NOT TO GIVE AN
C     IMPRESSION OF GENERALITY
C--------------------------------------------------------------CHMIMP
          IF(ROCVE.LE.SMALL) THEN
            IF(TFLAGE(IRE).EQ.1) iskip(ire)=1
            IF(TFLAGE(IRE).EQ.0) rocve=small
          END IF
C==============================================================CHMIMP
          PROG(NREA)=ZERO
          DPROG(NREA)=ZERO
          IPCON(IRE)=1
        ELSE
          ISKIP(IRE)=1
        END IF
   10 CONTINUE
C--------------------------------------------------------------CHMIMP
      N=0
      FACTOR=ONE
C-------------------------------------- THE ITERATION LOOP ----CHMIMP
   20 N=N+1
      ICONV=0
      IREA=0
      DO 100 IRE=1,NRE
       IF(ISKIP(IRE).EQ.0) THEN
        IREA=IREA+1
        PROG(IREA)=PROG(IREA)+FACTOR*DPROG(IREA)
        IF(NLTE.EQ.0) THEN
         DELTAT=(ROER(IC)*RRC(IC)-ROEN(IC))/ROCV
        ELSE
         DELTAT=(ONE-TFLAGE(IRE))*((ROER(IC)-ROEER(IC))*RRC(IC)
     1                              -ROEN(IC)+ROEEN(IC))/ROCVH
     2          +TFLAGE(IRE)*(ROEER(IC)*RRC(IC)-ROEEN(IC))/ROCVE
        END IF
        EQCC=EXP(DLNKDT(IRE)*DELTAT)*EQC(IRE)
        REQCC=ONE/EQCC
        KFEC=EXP(DLKFDT(IRE)*DELTAT)*KFE(IRE)
        RKFEC=ONE/KFEC
C--------------------------------------------------------------CHMIMP
        PRMAX=ZERO
        PRMIN=ZERO
        REF=ZERO
        RP=ONE
        PP=ONE
        RPC=ONE
        PPC=ONE
        ASS=ZERO
        BSS=ZERO
        ISUMA=0
        ISUMB=0
        ICONR=1
        NE=NLM(IRE)
C--------- FIND REFERENCE SPECIES AND CALCULATE PP AND RP -----CHMIMP
        DO 40 KK=1,NE
         ISP=CN(KK,IRE)
         TEST=SMALL*MW(ISP)*RC(IC)
         IF(N.EQ.1 .AND. SPDR(IC,ISP).LE.TEST) SPDR(IC,ISP)=TEST
         ROM=SPDR(IC,ISP)*RRC(IC)*RMW(ISP)
         AMBROM=-FBNAN(ISP,IRE)/ROM
         PRMAX=MAX(PRMAX,AMBROM)
         PRMIN=MIN(PRMIN,AMBROM)
         IF(IGAS(ISP).EQ.1) THEN
           TEST=ABS(AMBROM)
           IF(TEST.GT.REF) KREF=ISP
           REF=MAX(REF,TEST)
           RP=RP*ROM**AN(ISP,IRE)
           PP=PP*ROM**BN(ISP,IRE)
           FACT=FBNAN(ISP,IRE)*MW(ISP)*RC(IC)/SPDR(IC,ISP)
           ASS=ASS+AN(ISP,IRE)*FACT
           BSS=BSS+BN(ISP,IRE)*FACT
         ELSE
           RPC=RPC*ROM**AN(ISP,IRE)
           PPC=PPC*ROM**BN(ISP,IRE)
           ISUMA=ISUMA+AN(ISP,IRE)
           ISUMB=ISUMB+BN(ISP,IRE)
         END IF
   40   CONTINUE
C--------------------------------------------------------------CHMIMP
        ERAT=-QEQ(IRE)*(ONE-FRE(IRE))*RC(IC)/ROER(IC)
        PRMAX=MAX(PRMAX,ERAT)
        PRMIN=MIN(PRMIN,ERAT)
        IF(NLTE.NE.0) THEN
         ERAT=-(QEQ(IRE)*FEE(IRE)+FBNAN(IELC,IRE)*THRHAF*RGAS*TE(IC))
     1         *RC(IC)/ROEER(IC)
         PRMAX=MAX(PRMAX,ERAT)
         PRMIN=MIN(PRMIN,ERAT)
        END IF
C--------------------------------------------------------------CHMIMP
        PRMAX=ONE/PRMAX
        PRMIN=ONE/PRMIN
        AREF=AN(KREF,IRE)
        BREF=BN(KREF,IRE)
C--------------------------------------------------------------CHMIMP
        IF(NLTE.EQ.0) THEN
         ALFSS=DLKFDT(IRE)*QEQ(IRE)*(ONE-FRE(IRE))/ROCV
         BETASS=DLNKDT(IRE)*QEQ(IRE)*(ONE-FRE(IRE))/ROCV
        ELSE
         DTDWS=(QEQ(IRE)*(ONE-FRE(IRE)-FEE(IRE))
     1          -FBNAN(IELC,IRE)*THRHAF*RGAS*TE(IC))/ROCVH
         DTEDWS=(QEQ(IRE)*FEE(IRE)
     1          +FBNAN(IELC,IRE)*THRHAF*RGAS*TE(IC))/ROCVE
         QEFFS=(ONE-TFLAGE(IRE))*DTDWS+TFLAGE(IRE)*DTEDWS
         ALFSS=DLKFDT(IRE)*QEFFS
         BETASS=DLNKDT(IRE)*QEFFS
        END IF
C--------------------------------------------------------------CHMIMP
        PWR=ONE
        PROGT=RDTIC*RKFEC*PROG(IREA)
        IF(ABS(PROGT).GE.TENTH*REQCC*PP) IPCON(IRE)=0
        PROGTS=RDTIC*RKFEC*(PRMAX-PRMIN)
        SCALE=MAX(RP,REQCC*PP,ABS(PROGT),MCHACC*PROGTS)
        TEST=REQCC*PP+PROGT-RP
        IF(TEST.GT.EPSCHM*SCALE .AND. PPC.LT.MCHACC*CONC**ISUMB)
     1     ICONR=0
        IF(TEST.LT.-EPSCHM*SCALE .AND. RPC.LT.MCHACC*CONC**ISUMA)
     1     ICONR=0
        IF(ICONR.EQ.1 .AND. ABS(TEST).GT.EPSCHM*SCALE) ICONV=1
C--------------------------------------------------------------CHMIMP
        DIR=DBLE(AREF*(AREF-1))*RP
        DIP=DBLE(BREF*(BREF-1))*REQCC*PP
        IF(DIP.GT.DIR) THEN
         US=EQCC*(RP-PROGT)
         PPP=PP
         USPM=ONE
         IF(IPCON(IRE).EQ.1 .AND. BREF.GT.1) THEN
          PWR=ONE/DBLE(BREF)
          PPP=PP**PWR
          USPM=US**(PWR-ONE)
         END IF
         FS=PPP-USPM*US
         DUDOSS=EQCC*(RP*(ASS+BETASS)
     1                +RDTIC*RKFEC*(PROG(IREA)*(ALFSS-BETASS)-ONE))
         RDIAG=ONE/(PWR*(PPP*BSS-USPM*DUDOSS))
        ELSE
         VS=REQCC*PP+PROGT
         RPP=RP
         VSPM=ONE
         IF(IPCON(IRE).EQ.1 .AND. AREF.GT.1) THEN
          PWR=ONE/DBLE(AREF)
          RPP=RP**PWR
          VSPM=VS**(PWR-ONE)
         END IF
         FS=RPP-VSPM*VS
      DVDOSS=REQCC*PP*(BSS-BETASS)+RDTIC*RKFEC*(ONE-PROG(IREA)*ALFSS)
         RDIAG=ONE/(PWR*(RPP*ASS-VSPM*DVDOSS))
        END IF
        FS=FS*DBLE(ICONR)
C--------------------------------------------------------------CHMIMP
        IF(N.LE.NDIAG) THEN
C----------------- ONE-STEP GAUSS-SEIDEL-NEWTON ITERATION -----CHMIMP
         DPROG(IREA)=-OMGCHM*RDIAG*FS
         DPROG(IREA)=MIN(DPROG(IREA),OMGCM2*PRMAX)
         DPROG(IREA)=MAX(DPROG(IREA),OMGCM2*PRMIN)
C--------------------------------------------------------------CHMIMP
CDIR$ IVDEP
         DO 60 KK=1,NE
          ISP=CN(KK,IRE)
   60     SPDR(IC,ISP)=SPDR(IC,ISP)+MW(ISP)*FBNAN(ISP,IRE)
     1                 *DPROG(IREA)*RC(IC)
         ROER(IC)=ROER(IC)+QEQ(IRE)*(ONE-FRE(IRE))*DPROG(IREA)*RC(IC)
         IF(NLTE.NE.0) ROEER(IC)=ROEER(IC)+(QEQ(IRE)*FEE(IRE)
     1        +FBNAN(IELC,IRE)*THRHAF*RGAS*TE(IC))*DPROG(IREA)*RC(IC)
        ELSE
C----------------------------------- FOR NEWTON ITERATION -----CHMIMP
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
            IF(IGAS(ISP).NE.0) THEN
             FACT=FBNAN(ISP,M)*MW(ISP)*RC(IC)/SPDR(IC,ISP)
             AST=AST+AN(ISP,IRE)*FACT
             BST=BST+BN(ISP,IRE)*FACT
            END IF
   80      CONTINUE
C--------------------------------------------------------------CHMIMP
           IF(NLTE.EQ.0) THEN
            ALFST=DLKFDT(IRE)*QEQ(M)*(ONE-FRE(M))/ROCV
            BETAST=DLNKDT(IRE)*QEQ(M)*(ONE-FRE(M))/ROCV
           ELSE
            DTDWT=(QEQ(M)*(ONE-FRE(M)-FEE(M))
     1             -FBNAN(IELC,M)*THRHAF*RGAS*TE(IC))/ROCVH
            DTEDWT=(QEQ(M)*FEE(M)
     1             +FBNAN(IELC,M)*THRHAF*RGAS*TE(IC))/ROCVE
            QEFFT=(ONE-TFLAGE(IRE))*DTDWT+TFLAGE(IRE)*DTEDWT
            ALFST=DLKFDT(IRE)*QEFFT
            BETAST=DLNKDT(IRE)*QEFFT
           END IF
C--------------------------------------------------------------CHMIMP
           IF(DIP.GT.DIR) THEN
            DUDOST=EQCC*(RP*(AST+BETAST)
     1                   +RDTIC*RKFEC*PROG(IREA)*(ALFST-BETAST))
            BAM(IREA,MA)=PWR*(PPP*BST-USPM*DUDOST)*RDIAG*DBLE(ICONR)
           ELSE
            DVDOST=REQCC*PP*(BST-BETAST)-RDTIC*RKFEC*PROG(IREA)*ALFST
            BAM(IREA,MA)=PWR*(RPP*AST-VSPM*DVDOST)*RDIAG*DBLE(ICONR)
           END IF
C--------------------------------------------------------------CHMIMP
          END IF
         END IF
   90    CONTINUE
        END IF
      END IF
  100 CONTINUE
C==============================================================CHMIMP
      IF(N.LE.NDIAG) GO TO 180
C-------------------------------------------- LAPACK CALL -----CHMIMP
C     SUBROUTINE DGESVX SOLVES THE LINEAR SYSTEM A*X=B, USING LU
C     DECOMPOSITION.  HERE, A IS THE MATRIX BAM(L,M), X IS THE
C     VECTOR DPROG(L), AND B IS THE VECTOR RHS(L).
C-------------------------------------------- LAPACK CALL -----CHMIMP
C     SUBROUTINE DGELSS (LAPACK) SOLVES THE SAME SYSTEM A*X=B BY 
C     SINGULAR VALUE DECOMPOSITION, WHEN A IS SINGULAR 
C     (UNDERDETERMINED SYSTEM).  THIS ROUTINE IS NECESSARY FOR A
C     LINEARLY DEPENDENT SYSTEM OF KINETIC CHEMCIAL REACTIONS IN
C     EQUILIBRIUM LIMIT.  THIS ROUTINE RETURNS X IN THE SAME STORAGE
C     USED FOR B.
C--------------------------------------------------------------CHMIMP
      CALL DGESVX('N','N',NREA,1,BAM,LNRE,BAMF,LNRE,IPIV,EQUED,
     1            RWORK,CWORK,RHS,LNRE,DPROG,LNRE,RCOND,FERR,BERR,
     2            WORK,IWORK,INFO)
      IF(RCOND.LE.MCHACC) THEN
C--------------------------------------------------------------CHMIMP
C     THIS WRITE STATEMENT SHOULD BE DISABLED, IF SINGULAR MATRIX IS
C     EXPECTED -- LINEARLY DEPENDENT EQUILIBRIUM REACTIONS
C--------------------------------------------------------------CHMIMP
c       WRITE(19,910) NCYC,IC,N
        CAll DGELSS(NREA,NREA,1,BAM,LNRE,RHS,LNRE,SWORK,MCHACC,IRANK,
     1              WORK,5*LNRE,INFO)
        DO 110 IREA=1,NREA
         DPROG(IREA)=RHS(IREA)
 110    CONTINUE
        IF(INFO.NE.0) THEN
c         WRITE(19,*) NCYC,IC,'-- INFO=',INFO
          GO TO 200
        END IF
      END IF
C--------------------------------------------------------------CHMIMP
      DO 120 ISP=1,NSP
  120 DSPDR(ISP)=ZERO
      DROER=ZERO
      DROEER=ZERO
      IREA=0
      DO 140 IRE=1,NRE
       IF(ISKIP(IRE).EQ.0) THEN
        IREA=IREA+1 
        NE=NLM(IRE)
CDIR$ IVDEP
        DO 130 KK=1,NE
         ISP=CN(KK,IRE)
      DSPDR(ISP)=DSPDR(ISP)+MW(ISP)*FBNAN(ISP,IRE)*DPROG(IREA)*RC(IC)
  130   CONTINUE
       DROER=DROER+QEQ(IRE)*(ONE-FRE(IRE))*DPROG(IREA)*RC(IC)
       IF(NLTE.NE.0) DROEER=DROEER+(QEQ(IRE)*FEE(IRE)
     1        +FBNAN(IELC,IRE)*THRHAF*RGAS*TE(IC))*DPROG(IREA)*RC(IC)
       END IF
  140 CONTINUE
C--------------------------------------------------------------CHMIMP
      TEST=ZERO
      DO 150 ISP=1,NSP
      TEST=MIN(TEST,DSPDR(ISP)/SPDR(IC,ISP))
  150 CONTINUE
      TEST=MIN(TEST,DROER/ROER(IC))
      IF(NLTE.NE.0) TEST=MIN(TEST,DROEER/ROEER(IC))
      FACTOR=MIN(ONE,-OMGCM2/MIN(TEST,-SMALL))
      DO 170 ISP=1,NSP
  170 SPDR(IC,ISP)=SPDR(IC,ISP)+FACTOR*DSPDR(ISP)
      ROER(IC)=ROER(IC)+FACTOR*DROER
      IF(NLTE.NE.0) ROEER(IC)=ROEER(IC)+FACTOR*DROEER
C--------------------------------------------------------------CHMIMP
  180 IF(ICONV.EQ.0) GO TO 190
      IF(N.LT.NCHEM) GO TO 20
      ICMAX=IC
      NCELLS=NCELLS+1
  190 NSUM=NSUM+N
      DO 210 IRE=1,NRE
      NE=NLM(IRE)
      DO 210 KK=1,NE
        ISP=CN(KK,IRE)
        IF(DBLE(1-IGAS(ISP))*SPDR(IC,ISP)*RMW(ISP)*RRC(IC).GT.HUNDRD*
     1    MCHACC*CONC .AND. TIC.GT.TCUTEH(IRE)) WRITE(19,920) IC,NCYC
  210 CONTINUE
  200 CONTINUE
      IF(NZONES.EQ.0 .OR. NCELLS.EQ.0) RETURN
      AVG=DBLE(NSUM)/DBLE(NZONES)
c     WRITE(19,900) NCELLS,ICMAX,NCYC,NZONES,AVG
      WRITE(*,905) NCYC,NCELLS,AVG
C==============================================================CHMIMP
      RETURN
C
  900 FORMAT(' WARNING - ',I3,' CELL(S) NOT CONVERGED IN CHMIMP; LAST
     & IC =',I5,', CYCLE',I9/I5,' ACTIVE CELLS, AVG ITS=',F6.2)
  905 FORMAT(' CYCLE',I9,',',I5,' CELLS NOT CONVERGED, AVG ITS=',F6.2)
  910 FORMAT(' SINGULAR MATRIX IN CYCLE',I8,', CELL',I5,', ITER',I3)
  920 FORMAT(' CONDENSED PHASE PRESENT ABOVE CRITICAL TEMP AT CELL ',
     &       I5,' IN CYCLE ',I8)
      END
