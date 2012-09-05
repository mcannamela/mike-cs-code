*DECK FRAN
      DOUBLE PRECISION FUNCTION FRAN(DUMMY)
C----------------------------------------------------------------FRAN
C     PORTABLE PSEUDO-RANDOM NUMBER GENERATOR, FROM TONY WARNOCK,
C     CRAY RESEARCH INSTITUTE.
C     USE SINGLE PRECISION ON LONG WORD LENGTH MACHINES. QQ MUST
C     BE GIVEN TO MACHINE SIGNIFICANCE (5.9604644775390625E-8)
C     IF 2.**(-24) WON'T COMPILE ON CDC, USE 16704000000000000000B
C----------------------------------------------------------------FRAN
C     FRAN IS CALLED BY:  BREAK, COLIDE, INJECT, AND PMOVTV
C----------------------------------------------------------------FRAN
C               ANYONE WHO CONSIDERS ARITHMETICAL
C               METHODS OF PRODUCING RANDOM
C               DIGITS IS, OF COURSE, IN A STATE OF SIN.
C                         --- JOHN VON NEUMANN  (1951) ---
C----------------------------------------------------------------FRAN
      INCLUDE 'COML.h'
      DOUBLE PRECISION PP,QQ,QQQ,GB,GS,A,B,DUMMY
      PARAMETER (PP=2.0D0**24,QQ=2.0D0**(-24),QQQ=2.0D0**(-48),
     1           GB=1.0D0,GS=14390069.0D0)
C----------------------------------------------------------------FRAN
      A=GS*RANS+ONE
      B=GB*RANS+GS*RANB+INT(A*QQ)
      RANB=B-INT(B*QQ)*PP
      RANS=A-INT(A*QQ)*PP
      FRAN=(RANB*PP+RANS)*QQQ
C----------------------------------------------------------------FRAN
      RETURN
      END
