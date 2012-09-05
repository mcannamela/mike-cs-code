      SUBROUTINE HDFOUT
C=============================================================TECPLOT
C 
C   HDFOUT IS USED TO PREPARE DATA IN HDF FORMAT. 
C   Version 1.0 :  ONLY FOR THE 2-D CYLINDRICAL COORDINATE PLASMA FIELD
C     BY. Y.P. WAN, 10-24-97
C
C     CALLED BY LAVA
C=============================================================TECPLOT
      INCLUDE 'COML.h'
      INTEGER NX1,NY1,NPLANES,NROWS,NCOLS, narr
      DOUBLE PRECISION VMAX,VMIN
      CHARACTER(7) fstr
C=============================================================TECPLOT
c   
      NX1=NX+1
      NY1=NY+1
C ----- PRINT OUT DATA ACCORDING TO THE FORMAT OF fp2hdf
c  ---  only the temperature field is printed out as a test
      write(fstr, *) NCYC
      OPEN(95,FILE='Temp_'//fstr//'.dat')
     
      NPLANES=1
      NROWS=NX1
      NCOLS=NY1
c     number of arrays written to this file. 
      nArr = 7
	  
      WRITE(95,1005)TIME 
      WRITE(95,1002)NPLANES
      WRITE(95,1002)NROWS
      WRITE(95,1002)NCOLS
      WRITE(95,1002)nArr
      
      WRITE(95,1004)(X(I),I=1,NX1)
      WRITE(95,1004)(Y(J),J=1,NY1)
      WRITE(95,1004)((TEMP((J-1)*NXT+I),J=1,NY1),I=1,NX1)
      WRITE(95,1004)((V((J-1)*NXT+I),J=1,NY1),I=1,NX1)
      WRITE(95,1004)((U((J-1)*NXT+I),J=1,NY1),I=1,NX1)
c      WRITE(95,1004)((RO((J-1)*NXT+I),J=1,NY1),I=1,NX1)
c      WRITE(95,1004)((COND((J-1)*NXT+I),J=1,NY1),I=1,NX1)
c      WRITE(95,1004)((VISC((J-1)*NXT+I),J=1,NY1),I=1,NX1)

c     Argon mole fraction    
      WRITE(95,1004)(( spd((J-1)*NXT+I,1)*rmw(1)+spd((J-1)*NXT+I,7)
     1  *rmw(2) ,J=1,NY1),I=1,NX1)

c     Hydrogen mole fraction
      WRITE(95,1004)((spd((J-1)*NXT+I,3)*rmw(3)+
     1  half*(spd((J-1)*NXT+I,4)*rmw(4)
     2   + spd((J-1)*NXT+I,5)*rmw(5) ),J=1,NY1),I=1,NX1)

c     Nitrogen mole fraction
      WRITE(95,1004)((spd((J-1)*NXT+I,6)*rmw(6)+
     1 spd((J-1)*NXT+I,7)*rmw(7)
     2    +half*(spd((J-1)*NXT+I,8)*rmw(8)+
     3   spd((J-1)*NXT+I,9)*rmw(9)),
     4   J=1,NY1),I=1,NX1)

c     Oxygen mole fraction
      WRITE(95,1004)((spd((J-1)*NXT+I,10)*rmw(10)+half*
     1  (spd((J-1)*NXT+I,11)*rmw(11)+spd((J-1)*NXT+I,12)*rmw(12))
     2  ,J=1,NY1),I=1,NX1)
	 

C
 1001 FORMAT('TEXT')
 1002 FORMAT(I3)
 1003 FORMAT(E12.5)
 1004 FORMAT(6E12.5)
 1005 FORMAT(10E15.8)
C
      RETURN
      END
