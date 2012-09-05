*DECK REPACK
      SUBROUTINE REPACK
C==============================================================REPACK
C     THIS ROUTINE DESTROYS PARTICLES THAT HAVE LEFT THE SYSTEM
C     THROUGH AN OUTFLOW BOUNDARY.
C--------------------------------------------------------------REPACK
C     REAPCK AND DESTROTY ANY PARTICLES OUT OF THE SYSTEM.
C     HOWEVER, IF THERE IS ANY SYMMETRY PLANE, WE NEED TO PUT
C     PARTICLES BACK INTO THE DOMAIN.
C----------
C     CALLED BY LAVA
C==============================================================REPACK
      INCLUDE 'COML.h'
      INCLUDE 'COMMPTC.h'
C----------
      DOUBLE PRECISION RRP
      double precision ptype,dsplat,thicksp,re,we,ssol
C==============================================================REPACK
      NPN=0
C----------------------------- FOR AXISYMMETRIC WITH AXIS -----REPACK
       IF(ICYL.EQ.1 .AND. I3DP2D.EQ.0 .AND. X(1).LE.SMALL) THEN
        DO 100 IP=1,NP
        IF(YP(IP).GE.Y(NYS) .AND. YP(IP).LE.Y(NYL) .AND.
     1     ZP(IP).GE.Z(NZS) .AND. ZP(IP).LE.Z(NZL)) THEN
          IF(XP(IP).GE.X(NXS) .AND. XP(IP).LE.X(NXL)) THEN
            NPN=NPN+1
            XP(NPN)=XP(IP)
            YP(NPN)=YP(IP)
            ZP(NPN)=ZP(IP)
            UP(NPN)=UP(IP)
            VP(NPN)=VP(IP)
            WP(NPN)=WP(IP)
            ICP(NPN)=ICP(IP)
            RADP(NPN)=RADP(IP)
            PARTN(NPN)=PARTN(IP)
            TP(NPN)=TP(IP)
            EP(NPN)=EP(IP)
            DTPDEP(NPN)=DTPDEP(IP)
            AREAP(NPN)=AREAP(IP)
            PMASS(NPN)=PMASS(IP)
            TMM(NPN)=TMM(IP)
            TML(NPN)=TML(IP)
            EPM(NPN)=EPM(IP)
            EPL(NPN)=EPL(IP)
            RPSPHS(NPN)=RPSPHS(IP)
            RPSPHL(NPN)=RPSPHL(IP)
            EMSSP(NPN)=EMSSP(IP)
            MTYPE(NPN)=MTYPE(IP)
            DO 110 IPSP=1,NPSP
              PMSP(NPN,IPSP)=PMSP(IP,IPSP)
  110       CONTINUE
C --- ADD THE ARRAYS USED IN PARTICLE HEATING MODEL BY WAN 4/30/99
            RD(NPN)=RD(IP)
            RM(NPN)=RM(IP)
            RRS(NPN)=RRS(IP)
            VI(NPN)=VI(IP)
            VIRS(NPN)=VIRS(IP)
            SOLVS(NPN)=SOLVS(IP)
            SOLVL(NPN)=SOLVL(IP)
            SOLVRS(NPN)=SOLVRS(IP)
            DRD(NPN)=DRD(IP)
            DO 111 IPSP=1,NGRID
            TPCD(NPN,IPSP)=TPCD(IP,IPSP)
  111       CONTINUE
C---------- FOR PARTICLES PASSING THROUGH THE AXIS
          ELSEIF(XP(IP).LT.X(NXS) .AND. XP(IP).GE.-X(NXL)) THEN
            NPN=NPN+1
            XP(NPN)=-XP(IP)
            YP(NPN)=YP(IP)
            ZP(NPN)=ZP(IP)
            UP(NPN)=-UP(IP)
            VP(NPN)=VP(IP)
            WP(NPN)=WP(IP)
            ICP(NPN)=ICP(IP)
            RADP(NPN)=RADP(IP)
            PARTN(NPN)=PARTN(IP)
            TP(NPN)=TP(IP)
            EP(NPN)=EP(IP)
            DTPDEP(NPN)=DTPDEP(IP)
            AREAP(NPN)=AREAP(IP)
            PMASS(NPN)=PMASS(IP)
            TMM(NPN)=TMM(IP)
            TML(NPN)=TML(IP)
            EPM(NPN)=EPM(IP)
            EPL(NPN)=EPL(IP)
            RPSPHS(NPN)=RPSPHS(IP)
            RPSPHL(NPN)=RPSPHL(IP)
            EMSSP(NPN)=EMSSP(IP)
            MTYPE(NPN)=MTYPE(IP)
            DO 120 IPSP=1,NPSP
              PMSP(NPN,IPSP)=PMSP(IP,IPSP)
  120       CONTINUE
C --- ADD THE ARRAYS USED IN PARTICLE HEATING MODEL BY WAN 4/30/99
            RD(NPN)=RD(IP)
            RM(NPN)=RM(IP)
            RRS(NPN)=RRS(IP)
            VI(NPN)=VI(IP)
            VIRS(NPN)=VIRS(IP)
            SOLVS(NPN)=SOLVS(IP)
            SOLVL(NPN)=SOLVL(IP)
            SOLVRS(NPN)=SOLVRS(IP)
            DRD(NPN)=DRD(IP)
            DO 121 IPSP=1,NGRID
            TPCD(NPN,IPSP)=TPCD(IP,IPSP)
  121       CONTINUE
          END IF
        END IF
  100   CONTINUE
C------------------ FOR PSEUDO 3-D, WITH AND WITHOUT AXIS -----REPACK
      ELSEIF(I3DP2D.EQ.1) THEN
        DO 200 IP=1,NP
        RRP=SQRT(XP(IP)**2+ZP(IP)**2)
        IF(RRP.GE.X(NXS) .AND. RRP.LE.X(NXL) .AND.
     1     YP(IP).GE.Y(NYS) .AND. YP(IP).LE.Y(NYL)) THEN
          NPN=NPN+1
          XP(NPN)=XP(IP)
          YP(NPN)=YP(IP)
          ZP(NPN)=ZP(IP)
          UP(NPN)=UP(IP)
          VP(NPN)=VP(IP)
          WP(NPN)=WP(IP)
          ICP(NPN)=ICP(IP)
          RADP(NPN)=RADP(IP)
          PARTN(NPN)=PARTN(IP)
          TP(NPN)=TP(IP)
          EP(NPN)=EP(IP)
          DTPDEP(NPN)=DTPDEP(IP)
          AREAP(NPN)=AREAP(IP)
          PMASS(NPN)=PMASS(IP)
          TMM(NPN)=TMM(IP)
          TML(NPN)=TML(IP)
          EPM(NPN)=EPM(IP)
          EPL(NPN)=EPL(IP)
          RPSPHS(NPN)=RPSPHS(IP)
          RPSPHL(NPN)=RPSPHL(IP)
          EMSSP(NPN)=EMSSP(IP)
          MTYPE(NPN)=MTYPE(IP)
          DO 210 IPSP=1,NPSP
            PMSP(NPN,IPSP)=PMSP(IP,IPSP)
  210     CONTINUE
C --- ADD THE ARRAYS USED IN PARTICLE HEATING MODEL BY WAN 4/30/99
          RD(NPN)=RD(IP)
          RM(NPN)=RM(IP)
          RRS(NPN)=RRS(IP)
          VI(NPN)=VI(IP)
          VIRS(NPN)=VIRS(IP)
          SOLVS(NPN)=SOLVS(IP)
          SOLVL(NPN)=SOLVL(IP)
          SOLVRS(NPN)=SOLVRS(IP)
          DRD(NPN)=DRD(IP)
          DO 211 IPSP=1,NGRID
          TPCD(NPN,IPSP)=TPCD(IP,IPSP)
  211     CONTINUE
        ELSE
          WRITE(*,*)'NOTE, PARTICLE OUT OF DOMAIN,IP=',IP
        END IF
  200   CONTINUE
C----- FOR PLANAR 2-D, 3-D, AND AXISYMMETRIC WITHOUT AXIS -----REPACK
      ELSE
        DO 300 IP=1,NP
        IF(XP(IP).GE.X(NXS) .AND. XP(IP).LE.X(NXL) .AND.
     1     YP(IP).GE.Y(NYS) .AND. YP(IP).LE.Y(NYL) .AND.
     2     ZP(IP).GE.Z(NZS) .AND. ZP(IP).LE.Z(NZL)) THEN
          NPN=NPN+1
          XP(NPN)=XP(IP)
          YP(NPN)=YP(IP)
          ZP(NPN)=ZP(IP)
          UP(NPN)=UP(IP)
          VP(NPN)=VP(IP)
          WP(NPN)=WP(IP)
          ICP(NPN)=ICP(IP)
          RADP(NPN)=RADP(IP)
          PARTN(NPN)=PARTN(IP)
          TP(NPN)=TP(IP)
          EP(NPN)=EP(IP)
          DTPDEP(NPN)=DTPDEP(IP)
          AREAP(NPN)=AREAP(IP)
          PMASS(NPN)=PMASS(IP)
          TMM(NPN)=TMM(IP)
          TML(NPN)=TML(IP)
          EPM(NPN)=EPM(IP)
          EPL(NPN)=EPL(IP)
          RPSPHS(NPN)=RPSPHS(IP)
          RPSPHL(NPN)=RPSPHL(IP)
          EMSSP(NPN)=EMSSP(IP)
          MTYPE(NPN)=MTYPE(IP)
          DO 310 IPSP=1,NPSP
            PMSP(NPN,IPSP)=PMSP(IP,IPSP)
  310     CONTINUE
C --- ADD THE ARRAYS USED IN PARTICLE HEATING MODEL BY WAN 4/30/99
          RD(NPN)=RD(IP)
          RM(NPN)=RM(IP)
          RRS(NPN)=RRS(IP)
          VI(NPN)=VI(IP)
          VIRS(NPN)=VIRS(IP)
          SOLVS(NPN)=SOLVS(IP)
          SOLVL(NPN)=SOLVL(IP)
          SOLVRS(NPN)=SOLVRS(IP)
          DRD(NPN)=DRD(IP)
          DO 311 IPSP=1,NGRID
          TPCD(NPN,IPSP)=TPCD(IP,IPSP)
  311     CONTINUE
        END IF
  300   CONTINUE
      END IF
C--------------------------------------------------------------REPACK
      NP=NPN
C==============================================================REPACK
      RETURN
      END
*DECK REPACK
      SUBROUTINE REPACKEVAP
C==============================================================REPACKEVAP
C     THIS ROUTINE DESTROYS PARTICLES THAT HAVE BEEN EVAPORATED
C     BY Y.P. WAN 1-6-98
C--------------------------------------------------------------REPACKEVAP
C     CALLED BY LAVA
C==============================================================REPACKEVAP
      INCLUDE 'COML.h'
      INCLUDE 'COMMPTC.h'
C----------
      DOUBLE PRECISION RRP
C==============================================================REPACKEVAP
      NPN=0
C----------------------------- FOR AXISYMMETRIC WITH AXIS -----REPACKEVAP
      DO 100 IP=1,NP
        IF(RADP(IP).GT.1.D-7) THEN
            NPN=NPN+1
            XP(NPN)=XP(IP)
            YP(NPN)=YP(IP)
            ZP(NPN)=ZP(IP)
            UP(NPN)=UP(IP)
            VP(NPN)=VP(IP)
            WP(NPN)=WP(IP)
            ICP(NPN)=ICP(IP)
            RADP(NPN)=RADP(IP)
            PARTN(NPN)=PARTN(IP)
            TP(NPN)=TP(IP)
            EP(NPN)=EP(IP)
            DTPDEP(NPN)=DTPDEP(IP)
            AREAP(NPN)=AREAP(IP)
            PMASS(NPN)=PMASS(IP)
            TMM(NPN)=TMM(IP)
            TML(NPN)=TML(IP)
            EPM(NPN)=EPM(IP)
            EPL(NPN)=EPL(IP)
            RPSPHS(NPN)=RPSPHS(IP)
            RPSPHL(NPN)=RPSPHL(IP)
            EMSSP(NPN)=EMSSP(IP)
            MTYPE(NPN)=MTYPE(IP)
            DO 110 IPSP=1,NPSP
              PMSP(NPN,IPSP)=PMSP(IP,IPSP)
  110       CONTINUE
C --- ADD THE ARRAYS USED IN PARTICLE HEATING MODEL BY WAN 4/30/99
            RD(NPN)=RD(IP)
            RM(NPN)=RM(IP)
            RRS(NPN)=RRS(IP)
            VI(NPN)=VI(IP)
            VIRS(NPN)=VIRS(IP)
            SOLVS(NPN)=SOLVS(IP)
            SOLVL(NPN)=SOLVL(IP)
            SOLVRS(NPN)=SOLVRS(IP)
            DRD(NPN)=DRD(IP)
            DO 111 IPSP=1,NGRID
            TPCD(NPN,IPSP)=TPCD(IP,IPSP)
  111       CONTINUE
        ELSE
          WRITE(*,*)'NOTE, PARTICLE VAPORIZED,IP=',IP
        END IF
  100 CONTINUE
C--------------------------------------------------------------REPACKEVAP
      NP=NPN
C==============================================================REPACKEVAP
      RETURN
      END

