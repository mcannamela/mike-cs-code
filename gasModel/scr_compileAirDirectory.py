
from numpy import *
from scr_dangola_air_parameters_pickling import *
from scipy.interpolate import interp1d
import os
import cPickle as pkl

def writeLUT(fname, X, Y):
    with open(fname, 'w') as f:
        f.write(str(X[0])+' '+str(X[-1])+'\n')
        for i,y in enumerate(Y):
            f.write(str(y))
            if i<(Y.shape[0]-1):
                f.write(' ')
            else:
                f.write('\n')
        

DTC = dangolaThermalCond("dangola_fitParameters_air_thermalConductivity.csv")
DSH = dangolaSpecificHeat("dangola_fitParameters_air_specificHeat.csv")
DSE = dangolaSpecificEnthalpy("dangola_fitParameters_air_specificEnthalpy.csv")
DEC = dangolaElectricalConductivity("dangola_fitParameters_air_electricalConductivity.csv")
DV = dangolaViscosity("dangola_fitParameters_air_viscosity.csv")
#DD = dangolaDensity("dangola_fitParameters_air_meanMolarMass.csv")

DD = dict(zip(['thermalConductivity', 'specificHeat', 
          'specificEnthalpy', 'electricalConductivity', 
          'viscosity'],#, 'density'],
          [DTC, DSH, DSE, DEC, DV]))#, DD]))
          
n = 5000
T = linspace(0, 40000, n)
p = 1.

d = dict(zip(DD.keys(), [DD[k](p,T) for k in DD.keys()]))

N = 1600
h = linspace(0, d['specificEnthalpy'][-1], N)


Th = interp1d(d['specificEnthalpy'], 
                        T, bounds_error = False, 
                        fill_value = 0)(h)
        

dh = dict(zip(d.keys()+['temperature'], 
              [interp1d(d['specificEnthalpy'], 
                        d[k], bounds_error = False, 
                        fill_value = 0)(h)  for k in d.keys()]+[Th]))
h-= h[0]#h[1]+963
dh['specificHeat'][0] = dh['specificHeat'][1]
                        
del dh['specificEnthalpy']
#
for k in dh.keys():
    figure()
    if k =='electricalConductivity':
        semilogy(h, dh[k])
    else:
        plot(h, dh[k])
    title(k)
    xlabel('enthalpy')
    
figure()
plot(h, dh['thermalConductivity']/dh['specificHeat'])
xlabel('enthalpy')
title('thermal diffusivity, k/c')
    
airDir = os.path.join(os.path.curdir, 'air')
if not os.path.isdir(airDir):
    os.mkdir(airDir)
    
for k in dh.keys():
    writeLUT(os.path.join(airDir, k+'.LUT'), h, dh[k])
    print sum(isnan(dh[k]))
#    with open(os.path.join(airDir, k+'.LUT'), 'w') as f:
#        f.write(str(h[0])+' '+str(h[-1])+'\n')
#        for i,x in enumerate(dh[k]):
#            f.write(str(x))
#            if i<(dh[k].shape[0]-1):
#                f.write(' ')
#            else:
#                f.write('\n')

T_ = linspace(T[0], T[-1], N)
H_ = DSE(p, T_)
H_-=H_[0]#-h[0]
writeLUT(os.path.join(airDir, 'enthalpy'+'.LUT'), T_, H_)
#with open(os.path.join(airDir, 'enthalpy.LUT'), 'w') as f:
#    f.write(str(T_[0])+' '+str(T_[-1])+'\n')
#    for i,x in enumerate(H_):
#        f.write(str(x))
#        if i<(H_.shape[0]-1):
#            f.write(' ')
#        else:
#            f.write('\n')

            
pDict = dict(zip('enthalpy temperature viscosity specificHeat electricalConductivity thermalConductivity'.split(),
                 [h, Th, dh['viscosity'],dh['specificHeat'],dh['electricalConductivity'],dh['thermalConductivity']]))
with open('airPropertiesDict.pkl', 'wb') as f:
    pkl.dump(pDict, f, -1)
    

    
