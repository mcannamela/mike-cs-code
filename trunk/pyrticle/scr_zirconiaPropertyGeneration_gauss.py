import scipy
scipy.pkgload()
from numpy import *
from pylab import *

t = float64([0, 300, 900 ,1500 ,2100, 2700, 2951, 3026 ,3326, 20000])
conductivity = ones_like(t)*2.2
density = ones_like(t)*5819.
enthalpy = zeros_like(t)
cp = array([444, 444, 607, 676,944, 1227, 1346, 1381, 1400, 1400 ])
viscosity = ones_like(t)*100
heatOfFusion = 706300.
solidus_liquidus = array([2951., 3026.])
vapTemp = array([4500., 4700.])
heatOfVaporization = 5227000.

normdist = scipy.stats.norm(loc = mean(solidus_liquidus), scale = diff(solidus_liquidus)/3.)
vapDist = scipy.stats.norm(loc = mean(vapTemp), scale = diff(vapTemp)/4.)

props = array([t, conductivity, density, enthalpy, cp, viscosity])


n = 8000
P = zeros((6,n))
P[0] = linspace(0,20000, n)
N = 10.
smooth = lambda x:scipy.ndimage.filters.uniform_filter(x, N, mode = 'nearest')
for i in range(1,6):
	P[i] = scipy.interpolate.interp1d(props[0], props[i])(P[0])
	if i==4:
		#account for heat of fusion with specific heat
		#P[i][logical_and(P[0]>=solidus_liquidus[0],P[0]<=solidus_liquidus[1])] = 0 
		P[i] += heatOfFusion*normdist.pdf(P[0])
		
#		#account for vaporization with specific heat
#		P[i][logical_and(P[0]>=vapTemp[0],P[0]<=vapTemp[1])] = 0 
#		P[i] += heatOfVaporization*vapDist.pdf(P[0])
#		P[3] = cumsum(P[i])*(P[0][1]-P[0][0])
	
	P[i] = smooth(P[i])
	
# P[4]*=0
# P[4]+=1400
#P[3] = cumsum(P[i])*(P[0][1]-P[0][0])
P[3] = r_[atleast_1d([0]),scipy.integrate.cumtrapz(P[4], P[0])]
close('all')

figure()
plot(P[0], P[4], 'k-', label = 'with phase change')	
plot(t, cp, 'g--', label = 'without phase change')
plot(2951*ones(2),[0,12000],'b:', label = 'solidus')
plot(3026*ones(2),[0,12000],'r-.', label = 'liquidus')
xlim((0,5000))
legend(loc = 'upper left')
xlabel('temperature, K')
ylabel(r'specific heat, $\frac{J}{kg\cdot K}$')

figure()
plot(P[0], P[3])
M = 8
f = open('zirconiaGauss.txt','w')
for i in range(6):
	f.write(str(P[i][0:-1:M]).replace(',', ' ').replace('[', '').replace(']', '').replace('\n', '')+'\n')
	
f.close()