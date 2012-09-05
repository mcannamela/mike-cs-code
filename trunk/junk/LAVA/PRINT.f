*DECK PRINT
      SUBROUTINE PRINT(IPRINT)
C===============================================================PRINT
C     THIS ROUTINE PRINTS NECESSARY INFORMATIONS ALONG THE TIME
C     MARCHING.
C----------
C     IPRINT = 0 :  DEFINE FILES
C            = 1 :  PRINT RUN TIME MESSAGE ON SCREEN (EVERY JUMP CYCLES)
C            = 2 :  PRINT FOR MOVIE FILES ON lava.out (EVERY DTJUMP TIME)
C            = 3 :  PRINT MOVIE FILES AT THE END OF THE JOB
C            = 4 :  PRINT PARTICLE STATISTICS
C            = 11 : PRINT AT END OF JOB
C            =-1 :  NO PRINT
C----------
C     CALLED BY LAVA
C---------------------------------------------------------------PRINT
C     RESERVED UNITS FOR OUTPUT FILES ARE 10 -- 49
C===============================================================PRINT
      INCLUDE 'COML.h'
C----------
      INTEGER IPRINT
      INTEGER NMA
C----------
      DOUBLE PRECISION CTC,CTV
      EXTERNAL CTC,CTV
      CHARACTER*1 TAB
C----------
      double precision ptype,dgom,epdgom,ypn
      double precision rctot,yy,uu
      double precision xAr(n3),xN(n3),xH(n3),xO(n3),dN(n3),dH(n3),
     1     dO(n3),dAri(n3),dNi(n3),dHi(n3),dOi(n3),ta,talog
C
      TAB=CHAR(9)
C============================== NO PRINT OUT               =====PRINT
      IF(IPRINT.EQ.-1) RETURN
C============================== INITIAL STAGE - OPEN FILES =====PRINT
      IF(IPRINT.EQ.0) THEN
        OPEN(UNIT=61,FILE='l.imp50',STATUS='UNKNOWN',
     1       ACCESS='SEQUENTIAL',FORM='FORMATTED')
        OPEN(UNIT=62,FILE='l.imp75',STATUS='UNKNOWN',
     1       ACCESS='SEQUENTIAL',FORM='FORMATTED')
        OPEN(UNIT=63,FILE='l.imp100',STATUS='UNKNOWN',
     1       ACCESS='SEQUENTIAL',FORM='FORMATTED')
        OPEN(UNIT=64,FILE='l.imp125',STATUS='UNKNOWN',
     1       ACCESS='SEQUENTIAL',FORM='FORMATTED')
        OPEN(UNIT=65,FILE='l.imp150',STATUS='UNKNOWN',
     1       ACCESS='SEQUENTIAL',FORM='FORMATTED')
c       OPEN(UNIT=64,FILE='l.imp20',STATUS='UNKNOWN',
c    1       ACCESS='SEQUENTIAL',FORM='FORMATTED')
c       OPEN(UNIT=65,FILE='l.imp25',STATUS='UNKNOWN',
c    1       ACCESS='SEQUENTIAL',FORM='FORMATTED')
c       OPEN(UNIT=66,FILE='l.imp30',STATUS='UNKNOWN',
c    1       ACCESS='SEQUENTIAL',FORM='FORMATTED')
C=================== PRINT RUN-TIME MESSAGES ON THE SCREEN =====PRINT
      ELSEIF(IPRINT.EQ.1) THEN
        WRITE(*,9010) NCYC,TIME,DT(IC1),V(NXYT-NXT-NX-1)
C==================================== GENERATE MOVIE FILES =====PRINT
      ELSEIF(IPRINT.EQ.2) THEN
        WRITE(NFOUT,9020) NCYC,TIME,DT(IC1)
C=================================== AT THE END OF THE JOB =====PRINT
      ELSEIF(IPRINT.EQ.3) THEN
C---------------------------------------------------------------PRINT
      WRITE(NFOUT,9000) 
      WRITE(NFOUT,9030) TIME
      WRITE(NFOUT,*) 'J, I, IC, u, v, T, p, xAr, xH2, xN2, xO2'
      DO 330 K=1,NZT
      ICK=(K-1)*NXYT
      DO 320 J=1,NYT
      ICJ=ICK+(J-1)*NXT
      DO 310 I=1,NXT
      IC=I+ICJ
C---------------------------------------------------------------PRINT
      xAr(ic)=(SPD(ic,1)+SPD(ic,2))*RMW(1)
      xH(ic)=(SPD(ic,3)+SPD(ic,4)+SPD(ic,5))*RMW(4)
     1       +SPD(ic,14)*RMW(14)+two*SPD(ic,15)*RMW(15)
      xN(ic)=(SPD(ic,6)+SPD(ic,7)+SPD(ic,8)+SPD(ic,9))*RMW(8)
      xO(ic)=(SPD(ic,10)+SPD(ic,11)+SPD(ic,12))*RMW(11)
     1       +SPD(ic,14)*RMW(14)+SPD(ic,15)*RMW(15)
       rctot=one/(xAr(ic)+xH(ic)+xN(ic)+xO(ic))
       xAr(ic)=xAr(ic)*rctot
       xH(ic)=xH(ic)*rctot
       xN(ic)=xN(ic)*rctot
       xO(ic)=xO(ic)*rctot
      write(nfout,9032)
     1 j,i,ic,u(ic),v(ic),temp(ic),p(ic),xAr(ic),xH(ic),xN(ic),xO(ic)
c----- degree of non-equilibrium
      ta=THOUTH*temp(ic)
      talog=LOG(TA)
      dH(ic)=(spd(ic,4)*rmw(4))**2/(spd(ic,3)*rmw(3)
     1     *EXP(AKS(1)*TALOG+BKS(1)/TA+CKS(1)+TA*(DKS(1)+EKS(1)*TA)))
      dN(ic)=(spd(ic,6)*rmw(6))**2/(spd(ic,4)*rmw(4)
     1     *EXP(AKS(2)*TALOG+BKS(2)/TA+CKS(2)+TA*(DKS(2)+EKS(2)*TA)))
      dO(ic)=(spd(ic,9)*rmw(9))**2/(spd(ic,8)*rmw(8)
     1     *EXP(AKS(3)*TALOG+BKS(3)/TA+CKS(3)+TA*(DKS(3)+EKS(3)*TA)))
      dAri(ic)=spd(ic,2)*rmw(2)*spd(ic,13)*rmw(13)/(spd(ic,1)*rmw(1)
     1          *EXP(AS(1)*TALOG+BS(1)/TA+CS(1)+TA*(DS(1)+ES(1)*TA)))
      dHi(ic)=spd(ic,5)*rmw(5)*spd(ic,13)*rmw(13)/(spd(ic,4)*rmw(4)
     1          *EXP(AS(2)*TALOG+BS(2)/TA+CS(2)+TA*(DS(2)+ES(2)*TA)))
      dNi(ic)=spd(ic,9)*rmw(9)*spd(ic,13)*rmw(13)/(spd(ic,8)*rmw(8)
     1          *EXP(AS(3)*TALOG+BS(3)/TA+CS(3)+TA*(DS(3)+ES(3)*TA)))
      dOi(ic)=spd(ic,12)*rmw(12)*spd(ic,13)*rmw(13)/(spd(ic,11)*rmw(11)
     1          *EXP(AS(4)*TALOG+BS(4)/TA+CS(4)+TA*(DS(4)+ES(4)*TA)))
c-----
  310 CONTINUE
  320 CONTINUE
  330 CONTINUE
C-----------print info of particles in lava.out    -------------PRINT
      IF(NP.GE.1) THEN
        WRITE(NFOUT,9000)
        WRITE(NFOUT,*) 'IP,ICP,TYPE,XP,YP,ZP,UP,VP,WP,TP,DGOM,RADP'
        DO 360 IP=1,NP
         PTYPE=ONE
         IF(PMSP(IP,1)/PMASS(IP).LE.HALF) PTYPE=ZERO
        EPDGOM=MAX(EP(IP),EPM(IP))
        EPDGOM=MIN(EPDGOM,EPL(IP))
        DGOM=(EPDGOM-EPM(IP))/(EPL(IP)-EPM(IP))
          WRITE(NFOUT,9031) IP,ICP(IP),PTYPE,XP(IP),YP(IP),ZP(IP),
     1                      UP(IP),VP(IP),WP(IP),TP(IP),DGOM,RADP(IP)
  360   CONTINUE
      END IF
C---------------------------------------------------------------PRINT
      OPEN(UNIT=50,FILE='xy.dat',STATUS='UNKNOWN',
     1     ACCESS='SEQUENTIAL',FORM='FORMATTED')
      WRITE(50,*) 'Variables = x y u v w T p x_A_r x_H x_N x_O L'
      WRITE(50,9300) NX,NY
 9300 FORMAT('Zone t="0" I=',I2,' J=',I2)
C---------------------------------------------------------------PRINT
      do 350 j=2,nyt-1
      icj=(j-1)*nxt
      do 340 i=2,nxt-1
      ic=i+icj
      yy=half*(y(j)+y(j-1))
c--- velocity in m/s
      uu=half*(u(ic)+u(ic-1))*hundth
      flh(ic)=half*(v(ic)+v(ic-nxt))*hundth
c--- turbulent length scale
      te(ic)=CMU*tke(ic)*sqrt(tke(ic))/eps(ic)
      write(50,9021) xc(ic),yy,uu,flh(ic),w(ic)*hundth,
     1        temp(ic),p(ic)*1.0d-4,
     2        xAr(ic),xH(ic),xN(ic),xO(ic),te(ic)
  340 continue
  350 continue
C------------------------------------------ AXIAL PROFILES -----PRINT
      OPEN(UNIT=51,FILE='AxialVariable.dat',STATUS='UNKNOWN',
     1     ACCESS='SEQUENTIAL',FORM='FORMATTED')
C------    Field value axial distribution              ---------PRINT
      OPEN(UNIT=97,FILE='AxialDistrib.dat',STATUS='UNKNOWN',
     1       ACCESS='SEQUENTIAL',FORM='FORMATTED')
      WRITE(97,103)TAB,TAB,TAB
      WRITE(51,*) 'Variables = y v T p x_A_r x_H x_N x_O'
      WRITE(51,*) 'ZONE T="primary"'
      do 510 j=2,nyt-1
      ic=(j-1)*nxt+2
      write(51,9021) half*(y(j)+y(j-1)),flh(ic),temp(ic),p(ic)*1.0d-4,
     1  xAr(ic),xH(ic),xN(ic),xO(ic)
      write(97,104) half*(y(j)+y(j-1)),tab,flh(ic),tab,
     &  temp(ic),tab,p(ic)*1.0d-4
  510 continue
103     FORMAT('*'/'Y(cm)',A1,'V(cm/s)',
     &  A1,'Tp(K)',A1,'Press(bar)')
104     FORMAT(1PE10.3,3(A1,1PE10.3))
C---------- degree of non-quilibrium                  ----------PRINT
      OPEN(UNIT=52,FILE='AxialNonEquil.dat',STATUS='UNKNOWN',
     1     ACCESS='SEQUENTIAL',FORM='FORMATTED')
      WRITE(52,*)
     1   'Variables = y `z_H `z_N `z_O `z_A_r_i `z_H_i `z_N_i `z_O_i'
      WRITE(52,*) 'ZONE T="chemical"'
      do 520 j=2,nyt-1
      ic=(j-1)*nxt+2
      write(52,9021) half*(y(j)+y(j-1)),dH(ic),dN(ic),dO(ic),
     1               dAri(ic),dHi(ic),dNi(ic),dOi(ic)
  520 continue
C---------- distrib. of particle partial den.         ----------PRINT
      OPEN(UNIT=53,FILE='AxialPartDen.dat',STATUS='UNKNOWN',
     1     ACCESS='SEQUENTIAL',FORM='FORMATTED')
      WRITE(53,*)
     1  'Variables = y Ar Ar^+ H_2 H H^+ N_2 N_2^+ N N^+ e^- OH H_2O'
      WRITE(53,*) 'ZONE T="number"'
      do 530 j=2,nyt-1
      ic=(j-1)*nxt+2
      write(53,9021) half*(y(j)+y(j-1)),avogad*spd(ic,1)*rmw(1),
     1          avogad*spd(ic,2)*rmw(2),avogad*spd(ic,3)*rmw(3),
     2          avogad*spd(ic,4)*rmw(4),avogad*spd(ic,5)*rmw(5),
     3          avogad*spd(ic,6)*rmw(6),avogad*spd(ic,7)*rmw(7),
     4          avogad*spd(ic,11)*rmw(11),avogad*spd(ic,12)*rmw(12),
     5          avogad*spd(ic,13)*rmw(13),avogad*spd(ic,14)*rmw(14),
     6          avogad*spd(ic,15)*rmw(15)
  530 continue
C----------------------------------------- RADIAL PROFILES -----PRINT
c     OPEN(UNIT=52,FILE='RadialVariable.dat',STATUS='UNKNOWN',
c    1     ACCESS='SEQUENTIAL',FORM='FORMATTED')
C---------------------------------------------------------------PRINT
c     WRITE(52,*) 
c    1         'Variables = x T T_e p `z n_e'
c     xc(1)=zero
c     WRITE(52,*) 'ZONE T="0.5 cm"'
c     do 524 i=1,nxt-1
c     ic=(26-1)*nxt+i
c     write(52,9021) xc(ic),temp(ic),te(ic),p(ic)*1.0d-4,
c    1               eqsc(ic),xelec(ic)
c 524 continue
C================================= FOR PARTICLE STATISTICS =====PRINT
      ELSEIF(IPRINT.EQ.11) THEN
c ----  added by WAN, only to get the plot of particles
c         at end of job
c
        OPEN(UNIT=98,FILE='NiCrALY.dat',STATUS='UNKNOWN',
     1       ACCESS='SEQUENTIAL',FORM='FORMATTED')
        OPEN(UNIT=99,FILE='ZrO2.dat',STATUS='UNKNOWN',
     1       ACCESS='SEQUENTIAL',FORM='FORMATTED')
        WRITE(98,100)TAB,TAB,TAB,TAB,TAB,TAB,TAB,TAB
        WRITE(99,100)TAB,TAB,TAB,TAB,TAB,TAB,TAB,TAB
        DO 1120 IP=1,NP
        PTYPE=ONE
        IF(PMSP(IP,1)/PMASS(IP).LE.HALF) PTYPE=ZERO
        epdgom=max(ep(ip),epm(ip))
        epdgom=min(epdgom,epl(ip))
        dgom=(epdgom-epm(ip))/(epl(IP)-epm(ip))
        IF(PTYPE.EQ.ONE) THEN
          WRITE(98,102)YP(IP),TAB,XP(IP),TAB,ZP(IP),TAB,VP(IP),
     &    TAB,UP(IP),TAB,TP(IP),TAB,RADP(IP),TAB,PARTN(IP),
     &    TAB,DGOM,TAB,PMASS(IP)
        ELSE
          WRITE(99,102)YP(IP),TAB,XP(IP),TAB,ZP(IP),TAB,VP(IP),
     &    TAB,UP(IP),TAB,TP(IP),TAB,RADP(IP),TAB,PARTN(IP),
     &    TAB,DGOM,TAB,PMASS(IP)
        END IF     
 1120   CONTINUE
100     FORMAT('*'/'Y(cm)',A1,'X(cm)',A1,'Z(cm)',A1,'V(cm/s)',
     &  A1,'U(m/s)',A1,'Tp(K)',A1,'R(cm)',A1,'Num_Parti',
     &  A1,'Energ_Statu',A1,'Mass(g)')
102     FORMAT(1PE10.3,9(A1,1PE10.3))
C ----  END WAN
C
c --- Collect particle statistics during simulation
c
      ELSEIF(IPRINT.EQ.4) THEN
      DO 1110 IP=1,NP
        ypn=yp(ip)-dt(ic1)*vp(ip)
        PTYPE=ONE
        IF(PMSP(IP,1)/PMASS(IP).LE.HALF) PTYPE=ZERO
        epdgom=max(ep(ip),epm(ip))
        epdgom=min(epdgom,epl(ip))
        dgom=(epdgom-epm(ip))/(epl(IP)-epm(ip))
c-----
      IF(YPN.LT.FIVE .AND. YP(IP).GE.FIVE) THEN
      write(61,9110) ncyc,ip,ptype,xp(ip),zp(ip),up(ip),wp(ip),vp(ip),
     1               tp(ip),dgom,radp(ip),partn(ip)
      ELSEIF(YPN.LT.7.5d0 .AND. YP(IP).GE.7.5d0) THEN
      write(62,9110) ncyc,ip,ptype,xp(ip),zp(ip),up(ip),wp(ip),vp(ip),
     1               tp(ip),dgom,radp(ip),partn(ip)
      ELSEIF(YPN.LT.TEN .AND. YP(IP).GE.TEN) THEN
      write(63,9110) ncyc,ip,ptype,xp(ip),zp(ip),up(ip),wp(ip),vp(ip),
     1               tp(ip),dgom,radp(ip),partn(ip)
      ELSEIF(YPN.LT.12.5d0 .AND. YP(IP).GE.12.5d0) THEN
      write(64,9110) ncyc,ip,ptype,xp(ip),zp(ip),up(ip),wp(ip),vp(ip),
     1               tp(ip),dgom,radp(ip),partn(ip)
      ELSEIF(YPN.LT.15.0d0 .AND. YP(IP).GE.15.0d0) THEN
      write(65,9110) ncyc,ip,ptype,xp(ip),zp(ip),up(ip),wp(ip),vp(ip),
     1               tp(ip),dgom,radp(ip),partn(ip)
c     ELSEIF(YPN.LT.20.0d0 .AND. YP(IP).GE.20.0d0) THEN
c     write(64,9110) ncyc,ip,ptype,xp(ip),zp(ip),up(ip),wp(ip),vp(ip),
c    1               tp(ip),dgom,radp(ip),partn(ip)
c     ELSEIF(YPN.LT.25.0d0 .AND. YP(IP).GE.25.0d0) THEN
c     write(65,9110) ncyc,ip,ptype,xp(ip),zp(ip),up(ip),wp(ip),vp(ip),
c    1               tp(ip),dgom,radp(ip),partn(ip)
c     ELSEIF(YPN.LT.30.0d0 .AND. YP(IP).GE.30.0d0) THEN
c     write(66,9110) ncyc,ip,ptype,xp(ip),zp(ip),up(ip),wp(ip),vp(ip),
c    1               tp(ip),dgom,radp(ip),partn(ip)
       END IF
c-----
 1110 CONTINUE
C===============================================================PRINT
      END IF
C================================================= FORMATS =====PRINT
 9000 FORMAT(5X,'-----------------------------------------------')
 9010 FORMAT(1X,I8,1X,E12.5,1X,E20.14,1X,E20.14)
 9020 FORMAT(5X,'NCYC = ',I8,', TIME =',E12.5,', DT =',E12.5)
 9021 FORMAT(12(1X,E12.5))
 9030 FORMAT(5X,'TIME = ',E15.8,' SEC')
 9031 FORMAT(2I4,12D12.5)
 9032 FORMAT(3I4,12D12.5)
 9110 FORMAT(2I6,12D12.5)
      RETURN
      END
C===============================================================PRINT
c ----  the next two functions are not in use
      double precision function ctc(y0,i1,i,y,nt,t,n3)
      integer i1,i,nt,n3
      double precision y0,y(nt),t(n3),yh
      yh=0.5d0*(y(i1)+y(i1+1))
      if(y0.le.yh) then
        ctc=t(i-nt)+(t(i)-t(i-nt))
     1              *(2.0d0*y0-y(i1)-y(i1-1))/(y(i1+1)-y(i1))
      else
        ctc=t(i)+(t(i+nt)-t(i))
     1              *(2.0d0*y0-y(i1+1)-y(i1))/(y(i1+2)-y(i1+1))
      endif
      return
      end
C===============================================================PRINT
      double precision function ctv(y0,i1,i,y,nt,t,n3)
      integer i1,i,nt,n3
      double precision y0,y(nt),t(n3)
      ctv=t(i-nt)+(t(i)-t(i-nt))*(y0-y(i1))/(y(i1+1)-y(i1))
      return
      end
C===============================================================PRINT

