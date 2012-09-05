*DECK BLOCK DATA
      BLOCK DATA HOT
C=================================================================HOT
C     INITIALIZE DATA AND CONSTANTS
C----------
C     CALLED BY LAVA
C=================================================================HOT
      INCLUDE 'COML.h'
C----------
      INTEGER N
C
C=============================================== I=O CONTROL =====HOT
      DATA IDUMP /1/
      DATA NFDMP,NFRST,NFINP,NFOUT /11,12,13,14/
C========================================= CONTROL CONSTANTS =====HOT
      DATA TEMPMX /30000.0D0/
C-----------------------------------------------------------------HOT
C---------- CHEMISTRY ----------
      DATA TCHEMI /1.0D-10/
C---------- CHEMEQ ----------
      DATA EPSCHM /0.02D0/
      DATA OMGCHM,OMGCM2 /1.0D0,0.98D0/
C---------- WE MAY NOT USE CH4I'S AND NCMAX,FRAC,MINVRT ----------HOT
C---------- CHEME4 ----------
      DATA CH4I1,CH4I2 /25.0D0,2.0D0/
C---------- CHEMSS ----------
      DATA NCMAX,FRAC,MINVRT /25,0.20D0,0/
C-----------------------------------------------------------------HOT
C     NCNDNS = FLAG FOR CONDENSED PHASE
C-----------------------------------------------------------------HOT
      DATA NGAS,NCNDNS /1,0/
C================== SAFETY FACTORS USED IN TIME STEP CONTROL =====HOT
c ----  at the early time period 5*10&-6, NYCY<500, the controling dt 
c       is time growth. enlarge it to make quicker convergence
c ----  change the value of safedt, safdtd to 10 times of its old values
c       to enforce a quick convergence, by Y.P. WAN 7/9/97 
c       this data is now input in the input file
c      DATA SAFEDT,SAFDTD,DTGROW /0.25D0,0.80D0,1.02D0/
c      DATA SAFEDT,SAFDTD,DTGROW /0.6D0,0.8D0,1.02D0/
C================ SEED CONSTANTS FOR RANDOM NUMBER GENERATOR =====HOT
      DATA RANB,RANS /2396745.0D0,2396745.0D0/
C======================= CONSTANTS USED IN TURBULENCE MODELS =====HOT
      DATA CE1,CE2,CE3,CES,CMU
     1     /1.44D0,1.92D0,-1.0D0,1.50D0,0.09D0/
C---- CPS FOR TURBULENT DISPERSION
      DATA CPS /0.16432D0/
C---- FOR SGS MODEL
      DATA SGSA,SGSD /0.05D0,1.0D0/
c------- numbers suggested by Cloutman - does not seem to work
c     data sgsa,sgsd /0.117D0,1.4D0/
C---------- 1/PR FOR EPS DIFFUSION
      DATA RPRE /0.769231D0/
C---------- K-EPSILON CORRECTION
      DATA CE2CR,CMUCOR /0.0667D0,0.04D0/
      DATA AKEML1,AKEML2 /0.495D0,0.005D0/
C---------- CRITICAL REYNOLDS NUMBER FOR WALL FUNCTION
      DATA REYC /130.3D0/
C---- LAMINAR PRANDTL NUMBER, 1/PR ,1/SCHMIDT NUMBER FOR TURBULENCE
      DATA PRL,RPRT,RSCT /0.74D0,1.42857D0,1.42857D0/
C================================== PHYSICAL CONSTANTS (CGS) =====HOT
      DATA PI     /3.14159 26535 89793 23846D0/
      DATA PIO2   /1.57079 63267 94896 61923D0/
      DATA PI2    /6.28318 53071 79586 47692D0/
      DATA PI4   /12.56637 06143 59172 95384D0/
      DATA PI4O3  /4.18879 02047 86390 98461D0/
      DATA PIO180 /0.01745 32925 19932 18466D0/
C-----------------------------------------------------------------HOT
C     UNIVERSAL GAS CONSTANT = 8.3143D+7 ERG/(MOL*K)
      DATA RGAS /8.3143D+7/
C     AVOGADRO'S NUMBER = MOLECULES/G-MOL
      DATA AVOGAD /6.02252D+23/
C     ELEMENTARY CHARGE (ELECTRON) = 4.8029 ESU
      DATA ELECHG /4.8029D-10/
C     KCAL TO ERG CONVERSION = 4.184D+10 ERG/KCAL
      DATA ERGCAL /4.184D+10/
C     STEFAN-BOLTZMANN CONSTANT = 5.67D-5 ERG/(CM**2*K**4)
      DATA STEBOL /5.67D-5/
C     BOLTZMANN'S CONSTANT = 1.38D-16 ERG/K
      DATA BOLTZ /1.38D-16/
C     PRESSURE UNIT TRANSFER 1ATM=1013250.0d0 g/cm.s2
      DATA P1ATM /1013250.0d0/
C-----------------------------------------------------------------HOT
C     VISCOSITY RATIO
      DATA VISRAT /-0.66666667D0/
C     CONSTANT USED FOR DIFFUSION COEFFICIENT
      DATA ETA0 /200.0D0/
C================================================== NUMBERS ======HOT
      DATA ZERO,ONE,TWO,THREE,FOUR /0.0D0,1.0D0,2.0D0,3.0D0,4.0D0/
      DATA FIVE,SIX,SEVEN,EIGHT,NINE /5.0D0,6.0D0,7.0D0,8.0D0,9.0D0/
      DATA HALF,PNTHRE,PNTSIX,PNTNIN /0.5D0,0.3D0,0.6D0,0.9D0/
      DATA ONEO3,TWOO3 /0.33333333D0,0.66666667D0/
      DATA FOURO3,FIVO3 /1.33333333D0,1.66666667D0/
      DATA THRHAF,FIVHAF /1.5D0,2.5D0/
      DATA TEN,TENTH,HUNDRD,HUNDTH /10.0D0,0.1D0,100.0D0,0.01D0/
      DATA THOUSD,THOUTH /1.0D3,1.0D-3/
      DATA TWENTY /20.0D0/
      DATA SMALL,LARGE /1.0D-80,1.0D+100/
C=================================== MATHEMATICAL FUNCTIONS ======HOT
C     INVERSE OF ERROR FUNCTION, ERF(Y)**-1. INTERVAL OF Y IS 0.05
C-----------------------------------------------------------------HOT
      DATA (RERF(N),N=1,21)
     & /0.0D0, 4.4340387910005496D-02, 8.8855990494257692D-02,
     & 0.1337269216648197D0,0.1791434546212917D0,0.2253120550121781D0,
     & 0.2724627147267544D0,0.3208583217151815D0,0.3708071585935579D0,
     & 0.4226802386475592D0,0.4769362762044699D0,0.5341590877497024D0,
     & 0.5951160814499951D0,0.6608544253423844D0,0.7328690779592169D0,
     & 0.8134198475976184D0,0.9061938024368219D0,1.017902464832028D0,
     & 1.163087153676674D0,1.385903824349678D0,2.0D0/
C=================================================================HOT
C     VALUES DEFINED ABOVE ARE NOT PROBLEM DEPENDENT; USERS SHOULD
C     NOT CHANGE THOSE WITHOUT GOOD REASON.
C     VALUES DEFINED BELOW ARE PROBLEM DEPENDENT; USERS SHOULD DEFINE
C     THESE.
c     all this data are moved to the input file read by RINPUT by
c       Y.P. WAN, 10-21-97
C================================================== SET GRID =====HOT
C      DATA NX,X(1),XL,IGRIDX / 56,0.0D0, 6.00000D0,1/
C      DATA NY,Y(1),YL,IGRIDY / 65,0.0D0,15.00000D0,1/
C      DATA NZ,Z(1),ZL,IGRIDZ /  1,0.0D0, 1.00000D0,0/
C============================================ DEFAULT VALUES =====HOT
c  --- p = 1 atm for Brossa and Pfender, 0.85 atm for Fincke
c     DATA PAMB,TEMAMB /1013250.0D0,300.0D0/
C      data pamb,temamb /852000.0D0,300.0D0/
C---------- GRAVITIES
C      DATA GX,GY,GZ /0.0D0,0.0D0,0.0D0/
C=================================================================HOT
C     FILES CONTAINING TRANSPORT PROPERTY TABLE ARE IN THE
C     SUB-DIRECTORY NAMED XPTY. COMPONENT INDEX MUST BE CONSISTENT
C     WITH THE SUBROUTINE XPORT.
C-----------------------------------------------------------------HOT
C     INTERVALS ARE T=100(N-1). UNITS ARE IN CGS.
C=================================================================HOT
C     H2 AND N2 PLASMA TRANSPORT PROPERTY DATA,
C     OBTAINED FROM THE UNIVERSITY OF MINNESOTA
C=================================================================HOT
      INCLUDE 'XPTY/Ar.V'
c      INCLUDE 'XPTY/Ar.K'
      INCLUDE 'XPTY/mikesArgon.K'
      INCLUDE 'XPTY/Ar.R'
      INCLUDE 'XPTY/H2.V'
c      INCLUDE 'XPTY/H2.K'
      INCLUDE 'XPTY/mikesHydrogen.K'
      INCLUDE 'XPTY/H2.R'
      INCLUDE 'XPTY/Air.V'
c      INCLUDE 'XPTY/Air.K'
      INCLUDE 'XPTY/mikesAir.K'
      INCLUDE 'XPTY/Air.R'
C=================================================================HOT
C     HK ARRAYS ARE THE ENTHALPIES OF THE SPECIES, TAKEN FROM THE
C     JANAF THERMOCHEMICAL TABLES.  INTERVALS ARE T=100(N-1),
C     AND UNITS ARE KCAL/MOLE (LATER CONVERTED TO ERGS/GM IN RINPUT).
C     REFERENCE TEMPERATURE IS 0 K INSTEAD OF 298K OF JANAF TABLE.
C=================================================================HOT
      INCLUDE 'ENTHALPY/Ar.E'
      INCLUDE 'ENTHALPY/Ar+.E'
      INCLUDE 'ENTHALPY/H2.E'
c     INCLUDE 'ENTHALPY/H2+.E'
      INCLUDE 'ENTHALPY/H.E'
      INCLUDE 'ENTHALPY/H+.E'
      INCLUDE 'ENTHALPY/N2.E'
      INCLUDE 'ENTHALPY/N2+.E'
      INCLUDE 'ENTHALPY/N.E'
      INCLUDE 'ENTHALPY/N+.E'
      INCLUDE 'ENTHALPY/O2.E'
c     INCLUDE 'ENTHALPY/O2+.E'
      INCLUDE 'ENTHALPY/O.E'
      INCLUDE 'ENTHALPY/O+.E'
      INCLUDE 'ENTHALPY/e-.E'
      INCLUDE 'ENTHALPY/HO.E'
      INCLUDE 'ENTHALPY/H2O.E'
C=================================================================HOT
      END
