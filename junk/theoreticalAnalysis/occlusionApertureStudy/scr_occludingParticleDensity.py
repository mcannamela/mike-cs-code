#estimate the number of particles in the measurement volume at any given instant using mass flux, mean particle volume, density, mean particle velocity
from numpy import *
#needed for plume linear density
massFlux = .5e-3 #kg/s, ranges from .33 - .7
meanVelocity = 150 #m/s
meanVolumeDiameter = linspace(5e-6, 100e-6)
meanVolume = 4*pi*(meanVolumeDiameter)**3/3  #m^3

density = 6100 #kg/m^3
particleFlux = massFlux/(density*meanVolume)
linearParticleDensity = particleFlux/meanVelocity

#needed to get the sampling volume
plumeWidth = 3e-3 #m
aperture = linspace(

#occluding particles is number of particles in the sampling depth minus one for the occluded particle
Nstar = (linearParticleDensity*occludingDepth - 1)

close('all')
plot(meanVolumeDiameter, Nstar)
