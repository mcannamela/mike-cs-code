from numpy import *
from pylab import *
import os
import pdb
from gasses import *

def readLavaProperty(fname):
    with open(fname, 'r') as f:
        L = f.readlines()
    X = []
    for x in L:
        if x.find('&')!=-1:
            X+=x.replace('&','').replace(r'/','').replace('\n','').replace(',','').replace('D','e').split(' ')
    
    Y = []
    for x in X:
        if x=='':
            continue
        Y+=[x]
    T = arange(len(Y))*100.0
    
    return T, float64(Y)
    
if __name__=="__main__":
    oxyFile = os.path.join(os.curdir, 'ENTHALPY', 'O2.E')
    N2File = os.path.join(os.curdir, 'ENTHALPY', 'N2.E')
    airKFile = os.path.join(os.curdir, 'XPTY', 'Air.K')
    
    hFactor = 4.184e6 #from converting kcal/mol to J/kg, still needs molar mass in kg/kmol on bottom
    Moxy = 32.0
    MN2 = 28.0
    
    T, hoxy = readLavaProperty(oxyFile)
    hoxy*= hFactor/Moxy
    T, hN2 = readLavaProperty(N2File)
    hN2*= hFactor/MN2
    
    figure()
    plot(T, hoxy, 'r',label = 'LAVA, $O_2$')
    plot(T, hN2, 'b', label = 'LAVA, $N_2$')
    plot(T, .22*hoxy+.78*hN2, 'g', label ='LAVA, $O_2+N_2$')
    plot(T, Air.TFuns['enthalpy'](T), 'y', label = "D'Angola, air")
    xlabel('Temperature, K')
    ylabel('Enthalpy, J/kg')
    title('enthalpy of air in LAVA and in our code')
    legend(loc = 'upper left')
    
    T, kAir = readLavaProperty(airKFile)
    kAir*=1e-5
    figure()
    plot(T, kAir, 'r', label = 'LAVA, air')
    plot(T, Air.TFuns['thermalConductivity'](T), 'y', label = "D'Angola, air")
    xlabel('Temperature, K')
    ylabel('thermal conductivity, W/m/K')
    legend(loc = 'upper left')
    title('thermal conductivity of air in LAVA and in our code')
    
    
    figure()
    plot(T, Air.TFuns['specificHeat'](T))
    xlabel('Temperature, K')
    ylabel('specific heat, J/kg/K')
    title('specific heat of air')
    show()
    