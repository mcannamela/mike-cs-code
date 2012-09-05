import cPickle as pkl
import helpers as hp
from interpolators import cubicInterpolator as interpolator
import gasMixture as gm
import numpy as np
import scipy
scipy.pkgload()
import scipy.stats as stat
import mlabMacros as mlm
from plotMacros import *


with open('twoPeaksMappingSharpT.pkl', 'rb') as f:
    DD = pkl.load(f)

D = DD['D']
V = DD['V']
d = DD['d']
v = DD['v']

z = DD['z']
r = DD['r']
S = DD['S']

axesMin = array([d[0],v[0] ])
axesMax = array([d[-1],v[-1] ])
axesRange = axesMax-axesMin
f.close()

zirconia = gm.defaultMaterial


#diameter varies with row, injection velocity with column
Tmap = scipy.ndimage.gaussian_filter(S['surfaceTemperature'],1.,mode = 'nearest')
Hmap = scipy.ndimage.gaussian_filter(S['enthalpy']/(4*((D/2)**3)/(3*pi)),1.,mode = 'nearest')

temperatureMapper = interpolator(axesMin, axesMax, Tmap)
enthalpyMapper = interpolator(axesMin, axesMax, Hmap)

priorMeans = array([50e-6,20. ])
priorStdDev = array([22e-6, 3 ])

N = 50000
x = zeros((2,N))
x[0] = stat.norm(loc = priorMeans[0], scale = priorStdDev[0]).rvs(N)
#x[0] = hp.randomLognormal(N, priorMeans[0], priorMeans[0])
x[1] = stat.norm(loc = priorMeans[1], scale = priorStdDev[1]).rvs(N)

m = any(logical_and(x<axesMax[...,newaxis],
                    x>axesMin[...,newaxis]), axis = 0)
x = x[:,m]
print x.shape

T = temperatureMapper(x)
H = enthalpyMapper(x)

T+= 25*randn(T.shape[0])
close('all')

#fig()
#subplot(1,2,1)
#hist(Tmap.ravel(), 50)
#subplot(1,2,2)
#hist(Hmap.ravel(), 50)

#temperature histogram
figure(figsize = (20, 10))
#subplot(1,2,1)
hist(T,100, facecolor = 'k')
xlab('temperature, K')
ylab('count')
axisFontSize()
sf('twoPeaksTemperatureHistogram_computed')

#fig()
#hist(H/1000.,100, facecolor = 'k')
#xlab('enthalpy, J/g')
#ylab('count')
#axisFontSize()
#sf('twoPeaksEnthalpyHistogram_computed')

#prior distribution and temperature map
figure()
#subplot(1,2,2)
h,b0,b1 = histogram2d(x[0],x[1], bins = [100,101]) #prior distributions
h =scipy.ndimage.gaussian_filter(h, 1,mode = 'nearest')

contourf(d*1e6, v,  Tmap.T, 30, label = 'Temperature, K')
#colorbar()
contour(b0[:-1]*1e6, b1[:-1],  h.T, colors = 'k', linewidths = 2, label = 'p(v0,d)')
axis([3,d[-1]*1e6,1,v[-1]])
axisFontSize()

ylabel('injection velocity, m/s')
xlabel(r'diameter, $\mu$m')
title('Mapping of Particle States to Temperature:\nContours of Particle State Distribution Superposed on Temperature Mapping ')

sf('distributionSuperposedOnTemperature')

figure()
contourf(d*1e6, v,  Hmap.T, 50)
ylabel('injection velocity, m/s')
xlabel(r'diameter, $\mu$m')
title('Specific Enthalpy as a Function \n of Particle Injection Conditions' )
sf('SpecificEnthalpyField')

#fig()
#h,b0,b1 = histogram2d(x[0],x[1], bins = [100,101]) #prior distributions
#h =scipy.ndimage.gaussian_filter(h, 1,mode = 'nearest')
#contourf(v, d*1e6, Hmap, 30, label = 'Temperature, K')
#colorbar()
#contour(b1[:-1], b0[:-1]*1e6, h, colors = 'k', linewidths = 2, label = 'p(v0,d)')
#axis([1,v[-1],3,d[-1]*1e6])
#axisFontSize()
#
#xlab('injection velocity, m/s')
#ylab('diameter, microns')
#tit('Mapping of Particle States to Enthalpy:\nContours of Particle State Distribution Superposed on Temperature Mapping ')
#
#sf('distributionSuperposedOnEnthalpy')



mf = (Tmap-2950.)/75.
mf[mf>1.1]=1.1
mf[mf<-.1]=-.1


#fig()
#contourf(v, d*1e6, mf, 30)
#colorbar()
#axis([1,25,3,125])
#axisFontSize()
#
#xlab('injection velocity, m/s')
#ylab('diameter, microns')
#tit('Mapping of Particle States to Molten Fraction')
#
#sf('moltenFraction')
