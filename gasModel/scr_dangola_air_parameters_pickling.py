#read in the fit parameters from a text file, rearrange, and pickle
import cPickle as pkl
from numpy import *
import pdb
from pylab import *
#from D'angola and capitelli, "thermodynamic and transport properties in 
#equilibrium air plasmas in a wide pressure and temperature range"
#eur. phys. j. d 46, 129-150
class dangolaThermalCond(object):
    def __init__(self, fname):
        """
        read in the polynomial coeffs from the csv 
        """
        with open(fname, 'r') as f:
            self.a0 = float64(f.readline().split(',')[0])
            self.alpha = array([float64(L.split(',')) for L in f.readlines()])
        
        self.sigmaIdx = 5

        
    def __call__(self,pressure, temperature = linspace(0,20000, 201)):
        """
        compute the parameters at the given temperature and pressure
        pressure in atmospheres
        """
        self.pressure = pressure
        self.temperature = temperature
        self.conductivity = self.a0+zeros_like(temperature)
        n = self.alpha.shape[1]
        
        self.p = sum(self.alpha*(log(self.pressure)**(arange(n)))[newaxis,...], axis =1)
               
        self.a = self.p[slice(0, None, 3)]
        self.c = self.p[slice(1, None, 3)]
        self.delta = self.p[slice(2, None, 3)]
        
        self.a[1:] = exp(self.a[1:])
        self.c[1:] = exp(self.c[1:])
        self.delta[1:] = exp(self.delta[1:])
        
        self.a[2]*=-1
        self.a[7]*=-1
        
        for i in range(self.sigmaIdx):
            self.conductivity+= self.a[i]*self.sigma(log(temperature), self.c[i], self.delta[i])
        for i in range(self.sigmaIdx, self.a.shape[0]):
            self.conductivity+= self.a[i]*self.gamma(log(temperature), self.c[i], self.delta[i])
        
        self.conductivity = exp(self.conductivity)
        
        with open("airConductivity_p=%1.f_nT=%d.pkl"%(pressure, temperature.shape[0]), 'wb') as f:
            pkl.dump(dict(zip(['pressure', 'temperature', 'conductivity'],              
            [pressure, temperature, self.conductivity])), f, -1)
        
        return self.conductivity
        
    def gamma(self, T, c, delta):
        q = (T-c)/delta
        return exp(-q**2)
        
    def sigma(self, T, c, delta):
        q = (T-c)/delta
        return exp(q)/(exp(q)+exp(-q))

class dangolaSpecificHeat(object):
    def __init__(self, fname):
        """
        read in the polynomial coeffs from the csv 
        """
        with open(fname, 'r') as f:
            self.alpha = array([float64(L.split(',')) for L in f.readlines()])
        
        self.sigmaIdx = 5

        
    def __call__(self,pressure, temperature = linspace(0,20000, 201)):
        """
        compute the parameters at the given temperature and pressure
        pressure in atmospheres
        """
        self.pressure = pressure
        self.temperature = temperature
        
        n = self.alpha.shape[1]
        
        self.p = sum(self.alpha*(log(self.pressure)**(arange(n)))[newaxis,...], axis =1)
        
        self.const  = self.p[:2]
        self.a = self.p[slice(2+0, None, 3)]
        self.c = self.p[slice(2+1, None, 3)]
        self.delta = self.p[slice(2+2, None, 3)]
        
        self.a = exp(self.a)
        self.c = exp(self.c)
        self.delta = exp(self.delta)
        
        self.specificHeat = sum(self.const[:,newaxis]*temperature[newaxis,...]**arange(self.const.shape[0])[:,newaxis],axis = 0)
        for i in range(self.sigmaIdx):
            self.specificHeat+= self.a[i]*self.sigma(temperature, self.c[i], self.delta[i])
        for i in range(self.sigmaIdx, self.a.shape[0]):
            self.specificHeat+= self.a[i]*self.gamma(temperature, self.c[i], self.delta[i])
        
        self.specificHeat *= 4.184*1000.
        with open("airSpecificHeat_p=%1.f_nT=%d.pkl"%(pressure, temperature.shape[0]), 'wb') as f:
            pkl.dump(dict(zip(['pressure', 'temperature', 'specificHeat'],              
            [pressure, temperature, self.specificHeat])), f, -1)
        
        self.specificHeat[0] = self.specificHeat[1]
        return self.specificHeat
        
    def gamma(self, T, c, delta):
        q = (T-c)/delta
        return exp(-q**2)
        
    def sigma(self, T, c, delta):
        q = (T-c)/delta
        return exp(q)/(exp(q)+exp(-q))
        
class dangolaSpecificEnthalpy(dangolaSpecificHeat):
    def __init__(self, fname):
        """
        read in the polynomial coeffs from the csv 
        """
        with open(fname, 'r') as f:
            self.alpha = array([float64(L.split(',')) for L in f.readlines()])
        
        self.sigmaIdx = 7
        
    def __call__(self,pressure, temperature = linspace(0,20000, 201)):
        """
        compute the parameters at the given temperature and pressure
        pressure in atmospheres
        """
        self.pressure = pressure
        self.temperature = temperature
        
        n = self.alpha.shape[1]
        
        self.p = sum(self.alpha*(log(self.pressure)**(arange(n)))[newaxis,...], axis =1)
        
        self.const  = self.p[:2]
        self.a = self.p[slice(2+0, None, 3)]
        self.c = self.p[slice(2+1, None, 3)]
        self.delta = self.p[slice(2+2, None, 3)]
        
        self.a = exp(self.a)
        self.c = exp(self.c)
        self.delta = exp(self.delta)
        
        self.specificEnthalpy = sum(self.const[:,newaxis]*temperature[newaxis,...]**arange(1,self.const.shape[0]+1)[:,newaxis],axis = 0)
        for i in range(self.sigmaIdx):
            self.specificEnthalpy+= self.a[i]*self.sigma(temperature, self.c[i], self.delta[i])
        
        self.specificEnthalpy *= 4.184*1000.
        with open("airSpecificEnthalpy_p=%1.f_nT=%d.pkl"%(pressure, temperature.shape[0]), 'wb') as f:
            pkl.dump(dict(zip(['pressure', 'temperature', 'specificEnthalpy'],              
            [pressure, temperature, self.specificEnthalpy])), f, -1)
        
        return self.specificEnthalpy
        
class dangolaElectricalConductivity(dangolaSpecificHeat):
    def __init__(self, fname, minVal = 1e-4):
        """
        read in the polynomial coeffs from the csv 
        """
        with open(fname, 'r') as f:
            self.alpha = array([float64(L.split(',')) for L in f.readlines()])
        
        self.sigmaIdx = 7
        self.minVal = minVal        
        
    def __call__(self,pressure, temperature = linspace(0,20000, 201)):
        """
        compute the parameters at the given temperature and pressure
        pressure in atmospheres
        """
        self.pressure = pressure
        self.temperature = temperature
        
        n = self.alpha.shape[1]
        
        self.p = sum(self.alpha*(log(self.pressure)**(arange(n)))[newaxis,...], axis =1)
        
        self.xiPars  = self.p[:4]
        
        self.a = self.p[slice(4+0, None, 3)]
        self.c = self.p[slice(4+1, None, 3)]
        self.delta = self.p[slice(4+2, None, 3)]
        
        self.xiPars = exp(self.xiPars)        
        
        self.a[0] = exp(self.a[0])
        self.a[2:5]*=-1
        self.a[-1] = sum(self.a[2:5])- sum(self.a[:2])-self.a[5]
        self.c = exp(self.c)
        self.delta = exp(self.delta)
        
        self.electricalConductivity = self.xi(log(temperature), self.xiPars[0], self.xiPars[1],
                                                  self.xiPars[2],self.xiPars[3])
        for i in range(self.sigmaIdx):
            self.electricalConductivity+= self.a[i]*self.sigma(temperature, self.c[i], self.delta[i])
        
        #add minVal to the ec to keep 1/ec reasonable!
        self.electricalConductivity = self.minVal+exp(self.electricalConductivity)
        
        with open("airElectricalConductivity_p=%1.f_nT=%d.pkl"%(pressure, temperature.shape[0]), 'wb') as f:
            pkl.dump(dict(zip(['pressure', 'temperature', 'electricalConductivity'],              
            [pressure, temperature, self.electricalConductivity])), f, -1)
        
        return self.electricalConductivity
        
    def xi(self,T, a,c,delta, w ):
        return a-c*exp(-(T/delta)**w)
        
class dangolaViscosity(dangolaElectricalConductivity):
    def __init__(self, fname):
        """
        read in the polynomial coeffs from the csv 
        """
        with open(fname, 'r') as f:
            self.alpha = array([float64(L.split(',')) for L in f.readlines()])
        
        self.sigmaIdx = 5
        self.logSigmaIdx = 5
        
        
    def __call__(self,pressure, temperature = linspace(0,20000, 201)):
        """
        compute the parameters at the given temperature and pressure
        pressure in atmospheres
        """
        self.pressure = pressure
        self.temperature = temperature
        
        n = self.alpha.shape[1]
        
        self.p = sum(self.alpha*(log(self.pressure)**(arange(n)))[newaxis,...], axis =1)
        
        self.xiPars  = self.p[:4]
        
        self.a = self.p[slice(4+0, None, 3)]
        self.c = self.p[slice(4+1, None, 3)]
        self.delta = self.p[slice(4+2, None, 3)]
        
        self.a[1:5] = exp(self.a[1:5])
        self.a[6:9] = exp(self.a[6:9])
        self.a[5:9]*=-1        
        
        self.c[1:] = exp(self.c[1:])
        self.delta[1:-1] = exp(self.delta[1:-1])
        
        XI = self.xi(temperature, self.xiPars[0], self.xiPars[1],
                                                  self.xiPars[2],self.xiPars[3])
                                                  
        SIG = 0
        for i in range(self.logSigmaIdx):
            SIG += self.a[i]*self.sigma(temperature, self.c[i], self.delta[i])
            
        self.viscosity = log(XI+SIG)
            
        for i in range(self.logSigmaIdx, self.logSigmaIdx+self.sigmaIdx):
            self.viscosity+= self.a[i]*self.sigma(temperature, self.c[i], self.delta[i])
            
        self.viscosity = exp(self.viscosity)
        
        with open("airViscosity_p=%1.f_nT=%d.pkl"%(pressure, temperature.shape[0]), 'wb') as f:
            pkl.dump(dict(zip(['pressure', 'temperature', 'viscosity'],              
            [pressure, temperature, self.viscosity])), f, -1)
        
        return self.viscosity

class dangolaDensity(object):
    def __init__(self, fname):
        """
        read in the polynomial coeffs from the csv 
        """
        with open(fname, 'r') as f:
            self.alpha = array([float64(L.split(',')) for L in f.readlines()])
                
    def __call__(self,pressure, temperature = linspace(0,20000, 201)):
        """
        compute the parameters at the given temperature and pressure
        pressure in atmospheres
        """
        self.pressure = pressure
        self.temperature = temperature
        
        n = self.alpha.shape[1]
        
        self.p = sum(self.alpha*(log(self.pressure)**(arange(n)))[newaxis,...], axis =1)
        
        self.const  = atleast_1d(self.p[0])
        self.a = self.p[slice(1+0, None, 3)]
        self.c = self.p[slice(1+1, None, 3)]
        self.delta = self.p[slice(1+2, None, 3)]
        
        self.a = exp(self.a)
        self.c = exp(self.c)
        self.delta = exp(self.delta)
        
        self.density = sum(self.const[:,newaxis]*temperature[newaxis,...]**arange(self.const.shape[0])[:,newaxis],axis = 0)
        for i in range(self.a.shape[0]):
            self.density-= self.a[i]*self.sigma(temperature, self.c[i], self.delta[i])
        
        #convert from mean molar mass to density via ideal gas law
        self.density*=101325*pressure/(temperature*8.314172)
        
        #absolute zero correction
        self.density[isinf(self.density)] = np.max(self.density[logical_not(isinf(self.density))])
        
        with open("airDensity_p=%1.f_nT=%d.pkl"%(pressure, temperature.shape[0]), 'wb') as f:
            pkl.dump(dict(zip(['pressure', 'temperature', 'density'],              
            [pressure, temperature, self.density])), f, -1)
        
        return self.density
        
    def gamma(self, T, c, delta):
        q = (T-c)/delta
        return exp(-q**2)
        
    def sigma(self, T, c, delta):
        q = (T-c)/delta
        return exp(q)/(exp(q)+exp(-q))
        
if __name__=='__main__':
    DTC = dangolaThermalCond("dangola_fitParameters_air_thermalConductivity.csv")
    DSH = dangolaSpecificHeat("dangola_fitParameters_air_specificHeat.csv")
    DSE = dangolaSpecificEnthalpy("dangola_fitParameters_air_specificEnthalpy.csv")
    DEC = dangolaElectricalConductivity("dangola_fitParameters_air_electricalConductivity.csv")
    DV = dangolaViscosity("dangola_fitParameters_air_viscosity.csv")
    DD = dangolaDensity("dangola_fitParameters_air_meanMolarMass.csv")
   
    p = 1.
    T = linspace(0,30000,201)
   # T = linspace(0,50000,501)
    k = DTC(p,T)
    c = DSH(p,T)
    r = DD(p,T)
    mu = DV(p,T)
    sig = DEC(p,T)
    h = DSE(p,T)
    
    close('all')
    plot(T, k)
    title("thermal conductivity of air plasma")
    xlabel("T, K")
    ylabel(r"k, $\frac{W} {m^{2}K}$")
    
    figure()
    plot(T, c)
    title("specific heat of air plasma")
    xlabel("T, K")
    ylabel(r"c$_p$, $\frac{J} {kg-K}$")
    
    figure()
    plot(T, r)
    title("density of air plasma")
    xlabel("T, K")
    ylabel(r"$\rho$, $\frac {kg} {m^{3}}$")
    
    figure()
    plot(T, h)
    title("specific enthalpy of air plasma")
    xlabel("T, K")
    ylabel(r"h, $\frac {J} {kg}$")
    
    figure()
    plot(T, sig)
    title("electrical conductivity of air plasma")
    xlabel("T, K")
    ylabel(r"$\sigma$, $\frac {S} {m}$")
    
    figure()
    plot(T, mu)
    title("viscosity of air plasma")
    xlabel("T, K")
    ylabel(r"$\mu$, $\frac {kg} {m*s}$")