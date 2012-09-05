from gasModel import *


if __name__=="__main__":
    with open('ArH05PropertiesDict.pkl', 'rb') as f:
        DAr = pkl.load(f)
        DAr['enthalpy']-=DAr['enthalpy'][0]
    Ar = gas(38.0, DAr, 10000)
  
    with open('airPropertiesDict.pkl', 'rb') as f:
        DA = pkl.load(f)
        DA['enthalpy']-=DA['enthalpy'][0]
    A = gas(29.0, DA, 10000)
    
    f = linspace(0,1, 100)
    hMin = -Inf
    hMax = Inf
    for frac in f:
        G = Ar.blend(A, frac)
#        plot(G.T, G.DT['enthalpy'])
        hMin = {True:G.h[0], False:hMin}[G.h[0]>hMin]
        hMax = {True:G.h[-1], False:hMax}[G.h[-1]<hMax]
    
    h = linspace(hMin, hMax, Ar.nLUT)
    T = array([Ar.blend(A, frac).hFuns['temperature'](h) for frac in f])
    
    writeLUT2()('ArH05AirTemperature.LUT', f, h, T)
    
#    figure()
#    [plot(h, t) for t in T]
#    
    show()

