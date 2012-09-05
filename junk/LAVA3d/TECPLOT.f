      SUBROUTINE TECPLOT
C=============================================================TECPLOT
C
C   TECPLOT IS USED TO PREPARE INPUT DATA FOR THE GRAPHICS
C   PROGRAM "TECPLOT" ON UNIX SYSTEMS.
C   Version 1.0 :  ONLY FOR THE 2-D CYLINDRICAL COORDINATE PLASMA FIELD
C     BY. Y.P. WAN, 10-22-97
C     Modified RLW, 7-8-98
C     Modified RLW, 9-28-98
C
C     CALLED BY LAVA
C=============================================================TECPLOT
      INCLUDE 'COML.h'
      INTEGER NX1,NY1
      DOUBLE PRECISION EPDGOM
      DOUBLE PRECISION DGOM(NPAR)
C
      CHARACTER*13 FILENAM1, FILENAM2
      CHARACTER*8 PLAS, PART
      CHARACTER*4 PLT
      DATA PLAS,PART,PLT /'Plas_','Part_','.plt'/
C=============================================================TECPLOT
C
C     First the plasma flow data
C
      NX1=NX+1
      NY1=NY+1
      WRITE (FILENAM1,2000) PLAS,NMJ,PLT
c 
      OPEN(96,FILE=FILENAM1)
      WRITE(96,883) TIME
      WRITE(96,884)
      WRITE(96,886)NX1,NY1
      WRITE(96,888)((X(I),I=1,NX1),J=1,NY1)
      WRITE(96,888)((Y(J),I=1,NX1),J=1,NY1)
      WRITE(96,888)((U((J-1)*NXT+I),I=1,NX1),J=1,NY1)
      WRITE(96,888)((V((J-1)*NXT+I),I=1,NX1),J=1,NY1)
      WRITE(96,888)((TEMP((J-1)*NXT+I),I=1,NX1),J=1,NY1)
      WRITE(96,888)((RO((J-1)*NXT+I),I=1,NX1),J=1,NY1)
C      DO 10 ISP=1,NSP
C      WRITE(96,888)((SPD((J-1)*NXT+I,ISP),I=1,NX1),J=1,NY1)
C  10  CONTINUE
       ISP = 6
       WRITE(96,888)((SPD((J-1)*NXT+I,ISP),I=1,NX1),J=1,NY1)
c
      CLOSE(96)
C
C     Then the particle data (all particles in domain)
C
      WRITE (FILENAM2,2000) PART,NMJ,PLT
c 
      OPEN(97,FILE=FILENAM2)
      WRITE(97,983) TIME
      WRITE(97,984)
      WRITE(97,986) np
C
        DO 20 IP=1,NP
        EPDGOM=MAX(EP(IP),EPM(IP))
        EPDGOM=MIN(EPDGOM,EPL(IP))
        DGOM(IP)=(EPDGOM-EPM(IP))/(EPL(IP)-EPM(IP))
        WRITE(97,988) XP(IP),YP(IP),ZP(IP),
     1         UP(IP)*1.0D-2,VP(IP)*1.0D-2,WP(IP)*1.0D-2,
     2         TP(IP),DGOM(IP),RADP(IP)*2.0D4,PARTN(IP)
   20   CONTINUE
C
        CLOSE(UNIT=97)
C
C
      TIMOVI = TIMOVI + DTJUMP
      NMJ = NMJ + 1
C
C
  883 FORMAT(1x,'title="Plasma Flow Data   Time = ',e12.5,' s"')
  884 FORMAT(1x,'variables = "x","y","u","v","temp","dens","N2"')
  886 FORMAT(1x,'zone t= "zone 1" , i=',i5,' , j=',i3, ', f= block')
  888 FORMAT(6(1pe12.5,1x))
c
  983 FORMAT(1x,'title="Particle Data   Time = ',e12.5,' s"')
  984 FORMAT(1x,'variables = "x (cm)" "y (cm)" "z (cm)" ',
     1 '"u velocity (m/s)" "v velocity (m/s)" "w velocity (m/s)" ',
     2 '"Particle temperature (K)" "Particle Melt Fraction" ',
     3 '"Particle diameter (micron)" "PARTN"')
  986 FORMAT(1x,'zone',2x,'I = ',i5,2x,'f= point')
  988 FORMAT(10(1pe12.5,1x))
C
 2000 FORMAT(A5,I4.4,A4)
C
      RETURN
      END
