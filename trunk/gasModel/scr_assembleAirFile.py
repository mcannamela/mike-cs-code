import cPickle as pkl
from scipy.interpolate import interp1d
from scipy.integrate import cumtrapz
from numpy import *
import os 
from pylab import *
outName = 'air.txt'

gasDir = 'E:\\labCODE\\gasModel'
pklNames = ['airConductivity_p=1_nT=201.pkl',
                        'airDensity_p=1_nT=201.pkl',
                        'dummyEnthalpy.pkl',
                        'airSpecificHeat_p=1_nT=201.pkl',
                        'airViscosity.pkl']

fnames = [os.path.join(gasDir, name) for name in pklNames]

propNames = ['conductivity', 'density', 'enthalpy', 'specificHeat', 'viscosity']

with open('E:\\labCODE\\gasModel\\helium.txt') as f:
    T = float64(f.readline().split())
X = [T]
for i in range(len(fnames)):
    with open(fnames[i], 'rb') as f:
        D = pkl.load(f)
    X+= [interp1d(D['temperature'], D[propNames[i]]) (T) ]

    
X = array(X)
X[3 ]= r_[array([0]), cumtrapz(X[4], T)]

for i in range(1, X.shape[0]):
    figure()
    plot(T, X[i])
    xlabel('T')
    ylabel(propNames[i-1])

with  open(outName, 'wb') as f:
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            f.write(str(X[i,j])+' ')
        f.write('\n')
        
