######################################
##see shaw, "regular perturbation solution of the elenbaas-heller equation"
## in Journal of Applied Physics
######################################
from numpy import *
from scr_dangola_air_parameters_pickling import *
from scipy.interpolate import interp1d
import os
from scipy.integrate import cumtrapz
from scipy.stats import linregress

def regress(x,y, xLim):
    valid = (x>=min(xLim))*(x<=max(xLim))
    slope, intercept, r_value, p_value, std_err = linregress(x[valid],y[valid])
    
    def y_of_x(x_):
        return slope*x_+intercept
    
    return (slope, intercept, y_of_x, r_value)
    

DTC = dangolaThermalCond("dangola_fitParameters_air_thermalConductivity.csv")
DSH = dangolaSpecificHeat("dangola_fitParameters_air_specificHeat.csv")
DSE = dangolaSpecificEnthalpy("dangola_fitParameters_air_specificEnthalpy.csv")
DEC = dangolaElectricalConductivity("dangola_fitParameters_air_electricalConductivity.csv", minVal = 0.0)
DV = dangolaViscosity("dangola_fitParameters_air_viscosity.csv")
DD = dangolaDensity("dangola_fitParameters_air_meanMolarMass.csv")

DD = dict(zip(['thermalConductivity', 'specificHeat', 
          'specificEnthalpy', 'electricalConductivity', 
          'viscosity', 'density'],
          [DTC, DSH, DSE, DEC, DV, DD]))
          
n = 1600
T = linspace(1, 30000, n)
p = 1.

d = dict(zip(DD.keys(), [DD[k](p,T) for k in DD.keys()]))

K = d['thermalConductivity']
S = d['electricalConductivity']



TLin = float64([2400, 6600])

mT, bT, yTFun, RT = regress(1.0/T, log(S), sort(1.0/TLin))

Ti = -mT
S_Ahrr_T = exp(bT)*exp(-Ti/T)

figure()
plot(1.0/T, log(S),'k.')
plot(1.0/T, yTFun(1.0/T), 'y')
xlabel('inverse temperature')
ylabel('log sigma')
xlim(.8*min(1/TLin), 1.2*max(1/TLin))
ylim(.9*amin(log(S)), 1.1*amax(log(S)))
        
XI = r_[atleast_1d([.01]), cumtrapz(K, T)]

xiLin = float64([1/1.65e-5, 1/1.3e-4])
#xiLin = float64([1/1.3e-4, 1/1.2e-3])
#xiLin = float64([1/1.65e-5, 1/1.2e-3])
#xiLin = float64([1/1.2e-3, 1/.005])

mXI, bXI, yXIFun, RXI = regress(1.0/XI, log(S), sort(1.0/xiLin))

XIi = -mXI
S_Ahrr_XI = exp(bXI)*exp(-XIi/XI)


figure()
plot(1.0/XI[1:], log(S)[1:], 'k.')
plot(1.0/XI[1:], yXIFun(1.0/XI[1:]), 'y')
xlabel(r'inverse $\xi$')
ylabel('log sigma')
#xlim(.8*min(1/xiLin), 1.2*max(1/xiLin))
ylim(.9*amin(log(S)), 1.1*amax(log(S)))


figure()
plot(T, S_Ahrr_T,'y', label = 'T')
plot(T, S_Ahrr_XI,'r', label = r'$\xi$')
plot(T, S, 'k.')
xlabel('Temperature, K')
ylabel('electrical conductivity')
title('Ahrrenius fits for electrical conductivity')
#      
#
print "characteristic ionization xi is %.2e"%XIi
print "characteristic ionization temperature is %.2e"%Ti
print "R-squareds are %.2f, %.2f"%(RT**2, RXI**2)

#
#figure()
#plot(T, 1.0/XI)
#xlabel('T')
#ylabel(r'1/$\xi$')
#
#figure()
#plot(1.0/XI, log(S))
#xlabel(r'1/$\xi$')
#ylabel(r'ln($\sigma$)')
#
#
#figure()
#plot(x,y, 'y.')
#plot(x, intercept+slope*x,'k')
#xlabel(r'1/$\xi$')
#ylabel(r'ln($\sigma$)')
#
#
#
figure()
plot(T, XI)
xlabel('T')
ylabel(r'$\xi$')
title('thermal potential')
#
#figure()
#subplot(1,2,1)
#plot(XI, log(S))
#xlabel(r'$\xi$')
#ylabel(r'ln($\sigma$)')
#title('log electrical conductivity')
#subplot(1,2,2)
#plot(T, log(S))
#xlabel('T')
#ylabel(r'ln($\sigma$)')
#title('log electrical conductivity')
#
show()
