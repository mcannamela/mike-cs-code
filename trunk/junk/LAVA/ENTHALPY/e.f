c234567
      program e
C +++                                                                    
C +++ HK ARRAYS ARE THE ENTHALPIES OF THE SPECIES, TAKEN FROM THE        
C +++ JANAF THERMOCHEMICAL TABLES.  INTERVALS ARE T=100(N-1), AND UNITS  
C +++ ARE KJ/MOLE (GET CONVERTED TO KCal/mol).            
C +++                                                                    
C +++ -----------------------------------------------------              
C +++                                                                    
      implicit none
      character*7 afile
      real hk(61,3)
      double precision hkd(201,3),h1(201),h2(201),h3(201),
     1     dhk1,dhk2,dhk3
      integer i,j,n
      external oupu
c----- in KJ/mol
c----- for H2+
      DATA (HK(N,1),N=1,61)
     1 /-8.583, -5.739, -2.870, 0.054, 2.993, 5.968, 9.003, 12.115,
     2  15.310, 18.586, 21.938, 25.361, 28.848, 32.393, 35.990, 39.636,
     3  43.327, 47.060, 50.833, 54.643, 58.491, 62.374, 66.293, 70.247,
     4  74.236, 78.260, 82.318, 86.411, 90.537, 94.696, 98.886,103.105,
     5  107.352, 111.624, 115.919, 120.232, 124.562, 128.904, 133.256,
     6  137.612, 141.971, 146.328, 150.679, 155.021, 159.351, 163.666,
     7  167.962, 172.237, 176.489, 180.715, 184.913, 189.081, 193.218,
     8  197.322, 201.393, 205.428, 209.428, 213.392, 217.319, 221.209,
     9  225.063/
c----- for N2+
      DATA (HK(N,2),N=1,61)
     1 /-8.669, -5.768, -2.858, 0.054, 2.975, 5.927, 8.933, 12.006,
     2  15.153, 18.369, 21.650, 24.988, 28.379, 31.816, 35.297, 38.819,
     3  42.381, 45.984, 49.627, 53.311, 57.039, 60.811, 64.628, 68.492,
     4  72.402, 76.360, 80.364, 84.414, 88.510, 92.650, 96.883,101.056,
     5  105.318, 109.617, 113.951, 118.316, 122.710, 127.132, 131.578,
     6  136.048, 140.537, 145.045, 149.568, 154.106, 158.656, 163.216,
     7  167.786, 172.362, 176.944, 181.531, 186.122, 190.714, 195.308,
     8  199.902, 204.495, 209.088, 213.679, 218.268, 222.854, 227.436,
     9  232.015/
c----- for O2+
      DATA (HK(N,3),N=1,61)
     1 /-8.674, -5.771, -2.860, 0.054, 2.990, 5.979,
     2  9.045, 12.195, 15.426, 18.729, 22.095, 25.515, 28.981, 32.485,
     3  36.022, 39.588, 43.178, 46.789, 50.419, 54.065, 57.726, 61.399,
     4  65.085, 68.781, 72.486, 76.201, 79.923, 83.653, 87.390, 91.134,
     5  94.883, 98.639, 102.400, 106.166, 109.938, 113.714, 117.496,
     6  121.282, 125.073, 128.869, 132.670, 136.476, 140.288, 144.105,
     7  147.927, 151.756, 155.592, 159.434, 163.284, 167.142, 171.009,
     8  174.885, 178.771, 182.667, 186.575, 190.495, 194.428, 198.375,
     9  202.338, 206.316, 210.310/
c----- convert KJ/mol to Kcal/mol
      do 200 j=1,3
      do 100 i=1,61
      hkd(i,j)=dble(hk(i,j))/4.184d0
  100 continue
  200 continue
c----- extroplate above 6000K
      dhk1=hkd(61,1)-hkd(60,1)
      dhk2=hkd(61,2)-hkd(60,2)
      dhk3=hkd(61,3)-hkd(60,3)
      do 300 i=62,201
      hkd(i,1)=hkd(61,1)+dhk1*dble(i-61)
      hkd(i,2)=hkd(61,2)+dhk2*dble(i-61)
      hkd(i,3)=hkd(61,3)+dhk3*dble(i-61)
  300 continue
c----- output
      do 400 i=1,201
      h1(i)=hkd(i,1)-hkd(1,1)
      h2(i)=hkd(i,2)-hkd(1,2)
      h3(i)=hkd(i,3)-hkd(1,3)
c     h1(i)=h1(i)*4.184
c     h2(i)=h2(i)*4.184
c     h3(i)=h3(i)*4.184
c     write(21,*) float(i-1)*100.0,h1(i),h2(i),h3(i)
  400 continue
c     print*,dhk1*4.184d0,dhk2*4.184d0,dhk3*4.184d0
      print*,dhk1,dhk2,dhk3
      print*,h1(201)-h1(200),h2(201)-h2(200),h3(201)-h3(200)
c-----
      afile='H2+.E'
      call oupu (afile,h1)
      afile='N2+.E'
      call oupu (afile,h2)
      afile='O2+.E'
      call oupu (afile,h3)
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
 2100 format('      DATA (HK(N,11),N=1,50)')
 2200 format('      DATA (HK(N,11),N=51,100)')
 2300 format('      DATA (HK(N,11),N=101,150)')
 2400 format('      DATA (HK(N,11),N=151,201)')
      end
