from numpy import * 
from plotMacros import *

from scipy.integrate import cumtrapz

def planck(wavelength, T=None):
    """
    spherical intensity of radiated light at each wavelength and temperature, in W/steradian
    if T = None, the particle's current temperature profile is used
    
    wavelength - array of wavelengths
    T - array of temperatures
    """

    
    #physical constants, planck, speed o' light, boltzmann
    h=6.626e-34; #J*s
    c=3e8; #m/s
    kB= 1.381e-23;#J/K
    const=8*pi*h*c;
    
    #T varies along axis 0, wavelength along axis 1
    radiance=const*pow(expand_dims(wavelength,0),-5)*pow(exp(h*c/(outer(T,wavelength)*kB )) -1,-1) #W/m^3/steradians
    return radiance

#wavelengths are defined here, be sure they are wide enough to capture most of the energy
wave = 1e-9*linspace(1.,10000, 5000 )

#change the temperature here, in K
T = array([ 2500, 3000, 3500])


sr = squeeze(planck(wave, T = T))
csr = cumtrapz(sr,wave)
csr/=csr[:,-1][...,newaxis]

#change these to print the fraction of energy in a particular band
bandWavelengths = 1e-9*array([300, 1000])

fig()
fig()
for i in range(T.shape[0]):
    figure(1)
    plot(1e9*wave,sr[i], linewidth = 2, label = 'T = %d'%T[i])
    xlab('wavelength')
    ylab('spectral radiance')

    figure(2)
    plot(1e9*wave[:-1],csr[i], linewidth = 2,label = 'T = %d'%T[i])
    xlab('wavelength')
    ylab('spectral radiance')

figure(1)
legend()
figure(2)
legend()

print 'fractional energy in your band is: '
print str(csr[:,argmin(abs(wave-bandWavelengths[1]))]-csr[:,argmin(abs(wave-bandWavelengths[0]))])
