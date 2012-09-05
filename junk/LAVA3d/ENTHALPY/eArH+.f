c234567
      program eArH
C +++ INTERVALS ARE T=100(N-1), AND UNITS ARE KCal/mol.
C +++ -----------------------------------------------------              
      implicit none
      character*7 afile
      double precision h1(201),Rgas,ERGCAL
      integer i,j,n
      external oupu
c-----
      data Rgas/8.3143D+7/
C     KCAL TO ERG CONVERSION = 4.184D+10 ERG/KCAL
      DATA ERGCAL /4.184D+10/
c----- in Kcal/mol
c----- for ArH+
      do 400 i=1,201
      h1(i)=3.5d0*Rgas*dble(i-1)*100.0d0/ERGCAL
  400 continue
c-----
      afile='ArH+.E'
      call oupu (afile,h1)
c-----
      stop
      end
c-----
      subroutine oupu(afile,array)
      double precision array(201)
      character*7 afile
c-----
      open(unit=21,file=afile,status='unknown')
      write(21,2500)
      write(21,2100)
      write(21,2000) (array(n),n=1,50)
      write(21,2200)
      write(21,2000) (array(n),n=51,100)
      write(21,2300)
      write(21,2000) (array(n),n=101,150)
      write(21,2400)
      write(21,2000) (array(n),n=151,201)
      close(unit=21)
c-----
      return
 2500 format('C------------------------------------ ENTHALPTY T',
     1       'ABLE OF H2+ -----HOT',/,
     2       'C     EXTRAPOLATED ABOVE 6000 K',/,
     3       'C------------------------------------------------',
     4       '-----------------HOT')
 2000 format(41(5x,'&',1x,5(1x,p1d10.4,','),/))
 2100 format('      DATA (HK(N,10),N=1,50)')
 2200 format('      DATA (HK(N,10),N=51,100)')
 2300 format('      DATA (HK(N,10),N=101,150)')
 2400 format('      DATA (HK(N,10),N=151,201)')
      end

