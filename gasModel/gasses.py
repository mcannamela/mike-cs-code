from gasModel import *

with open('/media/raidArray/CODE/gasModel/argonPropertiesDict.pkl', 'rb') as f:
        DAr = pkl.load(f)
        DAr['enthalpy']-=DAr['enthalpy'][0]
Ar = gas(40.0, DAr, 10000)

    
with open('/media/raidArray/CODE/gasModel/airPropertiesDict.pkl', 'rb') as f:
    DA = pkl.load(f)
    DA['enthalpy']-=DA['enthalpy'][0]
    
Air = gas(29.0, DA, 10000)

with open('/media/raidArray/CODE/gasModel/ArH05PropertiesDict.pkl', 'rb') as f:
    DArH05 = pkl.load(f)
    DArH05['enthalpy']-=DArH05['enthalpy'][0]
ArH05 = gas(38.0, DArH05, 10000)

if __name__=="__main__":
    Ar.plotForPaper()
    Air.plotForPaper()
    
    figure()
    T = Ar.DT['temperature']
    idx = flatnonzero(T<=25000)[-1]
    plot(T[:idx]/1000., Ar.DT['specificHeat'][:idx], 'k')
    xlabel('Temperature, kK')
    ylabel(r'specific heat, $\frac{J}{kg \cdot K}$')
    xlim((0,25))    
    
    
