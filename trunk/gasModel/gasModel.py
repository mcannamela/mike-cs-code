from numpy import *
from pylab import *
from interpolators import linear1DGridInterpolator 
from scipy.interpolate import interp1d
import pdb
import cPickle as pkl

def polator(x,y):
    return interp1d(x,y)

def massFraction2molFraction(f, Mf, M):
    """
    convert mass fraction to a mole fraction 
    Mf - molar mass of the gas with mass fraction f
    M - molar mass of the other gas
    """
    return f/( (Mf/M) + f*(1-Mf/M) )
    
def writeLUT(fname, X, Y):
    with open(fname, 'w') as f:
        f.write(str(X[0])+' '+str(X[-1])+'\n')
        for i,y in enumerate(Y):
            f.write(str(y))
            if i<(Y.shape[0]-1):
                f.write(' ')
            else:
                f.write('\n')
                
class writeLUT2:
    def __call__(self,fname, x0, x1, Y):
        with open(fname, 'w') as f:
            f.write('%.2e   %.2e   %d\n'%(x0[0], x0[-1], len(x0)))
            f.write('%.2e   %.2e   %d\n'%(x1[0], x1[-1], len(x1)))
            for i,y in enumerate(Y.T):
                f.write(self.line(y))
    def line(self, X):
        """
        format a 1-d array into one line of text
        """
        s = ''
        for x in X:
            s+='%e '%x
        s+='\n'
        return s

class gas(object):
    def __init__(self, molarMass, propertiesDict, nLUT):
        """
        initialize with 
        molarMass - molar mass of the gas in kg/kmol
        propertiesDict - a dictionary of T, h, c_p, sigma, k, mu
        nLUT - number of entries in the lookup-tables
        """
        self.M = molarMass
        self.R = 8314.0        
        self.rhoConstant = 100000.0*self.M/self.R
        
        self.nLUT = nLUT
        self.rawDict = propertiesDict
        self.Dh, self.DT, self.TFuns, self.hFuns = self.cookProperties(self.rawDict)
    
    def __getstate__(self):
        """
        need to exclude all lambda objects
        """
        black =['TFuns', 'hFuns']
        S = []
        K = []
        for k in self.__dict__.keys():
            if k in black:
                continue
            K+=[k]
            S+=[self.__dict__[k]]
            
        return dict(zip(K,S))
        
    def __setstate__(self, D):
        """
        gonna have to call the methods that reconstruct the lambdas, which cannot be pickled
        """
        self.__dict__=D.copy()
        self.Dh, self.DT, self.TFuns, self.hFuns = self.cookProperties(self.rawDict)
        
    def cookProperties(self, D):
        Traw = D['temperature']
        hraw = D['enthalpy']
        
        self.T = linspace(Traw[0], Traw[-1], self.nLUT)
        self.h = linspace(hraw[0], hraw[-1], self.nLUT)
        
        TFuns = []
        K = []
        for k in D.keys():
            if k=='temperature':
                continue
            K+=[str(k)]
            p = polator(Traw, D[k])
            TFuns+=[polator(self.T, p(self.T))]
        
        DT = dict(zip(K+['temperature'], [f(self.T) for f in TFuns]+[self.T]))
        DTFun = dict(zip(K, TFuns))
        
        
        hFuns = []
        K = []
        for k in D.keys():
            if k=='enthalpy':
                continue
            K+=[str(k)]
            p = polator(hraw, D[k])
            hFuns+=[polator(self.h, p(self.h))]
        Dh = dict(zip(K+['enthalpy'], [f(self.h) for f in hFuns]+[self.h]))
        DhFun = dict(zip(K, hFuns))
        return Dh, DT, DTFun, DhFun
        
            
    
    def blend(self, other, f):
        """
        blend this gas with another according to mass fraction f
        other - another gas object the properties of which we will blend
        f - mass fraction of this gas. the fraction of other is therefore 1-f
        """
        TB = temperatureBlender(self.TFuns['enthalpy'], other.TFuns['enthalpy'])   
        
        if self.T[-1] > other.T[-1]:
            T   = other.T
            DT  = other.DT
            fun = self.TFuns
            otherFlag = True
        else:
            T   = self.T
            DT  = self.DT
            fun = other.TFuns
            otherFlag = False

        h = TB(T, f)
        Mbar = (f/self.M+(1-f)/other.M)**-1
        
        K = []
        P = []
        for k in fun.keys():
            K+=[k]
            if otherFlag:
                P+=[(1-f)*DT[k]+f*fun[k](T)]                
            else:
                P+=[f*DT[k]+(1-f)*fun[k](T)]                
                
        D = dict(zip(K+['temperature'],P+[T.copy()]))
        D['enthalpy'] = h
        return gas(Mbar, D, self.nLUT)
    
    def molarMass(self):
        return double(self.M)
    def Density(self, T):
        return self.rhoConstant/T
    def Conductivity(self, T):
        return self.TFuns['thermalConductivity'](T)
    def SpecificHeat(self, T):
        return self.TFuns['specificHeat'](T)
    def Viscosity(self, T):
        return self.TFuns['viscosity'](T)

    def plotT(self):
        figure()
        T = self.DT['temperature']
        for cnt, k in enumerate(self.DT.keys()):
            if k!='temperature':
                subplot(2,3,cnt+1)
                plot(T/1000., self.DT[k], 'g')
                xlabel('T, kK')
                ylabel(k)
                
    def plotForPaper(self):
        figure()
        T = self.DT['temperature']
        K = ['enthalpy',  'specificHeat', 'viscosity', 'electricalConductivity', 'thermalConductivity']
        yLabels = [r'specific enthalpy, $\frac{J}{K}$',  r'specific heat, $\frac{J}{kg \cdot K}$', 
                   r'viscosity, $\frac{kg}{m\cdot s}$', r'electrical conductivity, $\frac{S}{m}$', 
                   r'thermal conductivity, $\frac{W}{m\cdot K}$']
        idx = flatnonzero(T<=25000)[-1]
        for cnt, k in enumerate(K):
            if k!='temperature':
                subplot(3,2,cnt+1)
                plot(T[:idx]/1000., self.DT[k][:idx], 'k')
                xlabel('T, kK')
                ylabel(yLabels[cnt])
                xlim((0,25))
                
    def ploth(self):
        figure()
        h = self.Dh['enthalpy']
        for cnt, k in enumerate(self.Dh.keys()):
            if k!='enthalpy':
                subplot(2,3,cnt+1)
                plot(h, self.Dh[k], 'g')
                xlabel('enthalpy, J/kg')
                ylabel(k)
        
        
        
        
class temperatureBlender(object):
    def __init__(self, h1Fun, h2Fun):
        """
        h1Fun, h2Fun - callables that take as an argument a vector of temperatures
                        and return a vector of specific enthalpies for the component gasses 1 and 2
        M1, M2 - molar masses of each component gas
                        
        """
        self.h1 = h1Fun
        self.h2 = h2Fun
        
    def __call__(self, T, f):
        """
        T - vector of temperatures at which to compute the specific enthalpy
        f - mass fraction of component 1
        """
        
        
        h = f*self.h1(T)+(1-f)*self.h2(T)
        
        return h

class binaryGasMixture(object):
    def __init__(self, g1, g2):
        
        self.f = linspace(0,1, 100)
        self.g1 = g1
        self.g2 = g2
        self.G = array([g1.blend(g2, frac) for frac in self.f])

    def __call__(self, T,f):
        g = self.g(f)
        return (g.Density(T),
                g.Viscosity(T),
                g.Conductivity(T),
                g.SpecificHeat(T),
                g.molarMass())
    
    def f2idx(self, f):
        return amin([amax([int32(f*100)]),99])
    def g(self, f):
        return self.G[self.f2idx(f)]
        
    def meanMolarMass(self, f = 0):
        return self.G[self.f2idx(f)].molarMass()
    def Density(self, T, f):
        return self.g(f).Density(T)
    def Conductivity(self, T, f):
        return self.g(f).Conductivity(T)
    def SpecificHeat(self, T, f):
        return self.g(f).SpecificHeat(T)
    def Viscosity(self, T, f):
        return self.g(f).Viscosity(T)
        
       
    
if __name__=="__main__":
    with open('argonPropertiesDict.pkl', 'rb') as f:
        DAr = pkl.load(f)
        DAr['enthalpy']-=DAr['enthalpy'][0]
    Ar = gas(40.0, DAr, 10000)
    
#    Ar.plotT()
#    Ar.ploth()
    
    with open('airPropertiesDict.pkl', 'rb') as f:
        DA = pkl.load(f)
        DA['enthalpy']-=DA['enthalpy'][0]
        
    A = gas(29.0, DA, 10000)
    
#    A.ploth()
#    A.plotT()
    
    f = .5
    G = Ar.blend(A, f)
#    G.plotT()

    figure()
    plot(A.T, Ar.TFuns['enthalpy'](A.T), 'b')
    plot(A.T, A.DT['enthalpy'], 'r')
    plot(G.T, G.DT['enthalpy'], 'g')
    plot(f*Ar.hFuns['temperature'](Ar.h)+(1-f)*A.hFuns['temperature'](Ar.h), Ar.h, 'y')

    show()
    