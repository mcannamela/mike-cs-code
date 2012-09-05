from numpy import *
#from scr_compileAirDirectory import writeLUT
from scipy.interpolate import interp1d
import os
from pylab import *
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
                
def correctDiscontinuity(x, a,b):
    x[(b+1):] -= (x[b+1]-x[b])-(x[b+2]-x[b+1])
    x[(a+1):] -= (x[a+1]-x[a])-(x[a+2]-x[a+1])
    return x
    

    
def fixBump(x,c):
    x[c] = .5*(x[c-1]+x[c+1])
    return x
    
def extend(X, T, T1):
    return r_[X[0], X, X[-1]+(T1-T[-1])*(X[-1]-X[-2])/(T[-1]-T[-2])]
                
gasFolder = "/media/raidArray/CODE/gasModel"
murphyFolder = 'murphy_Ar_properties'
LUTFolder = 'argon'
if not os.path.isdir(LUTFolder):
    os.mkdir(os.path.join(os.path.curdir, LUTFolder))

thermoFile = 'ar_1atm.thd'
transportFile = 'ar_1atm.tra'

dataStartLine = 6
N = 1000

thermoIdxDict = {'T':1, 'enthalpy':5, 'specificHeat':6}
transportIdxDict = {'T':1, 'viscosity':2, 'thermalConductivity':3, 'electricalConductivity':4 }

Textend = 50000

with open(os.path.join(os.path.curdir,murphyFolder, thermoFile),'r') as f:
    D = f.readlines()[5:]
    X = float64([L.split() for L in D]).T
    T = X[thermoIdxDict['T']]
    H = extend(X[thermoIdxDict['enthalpy']], T, Textend)
   
    Cp = extend(X[thermoIdxDict['specificHeat']], T, Textend)
    Cp[-1] = Cp[-2]
    
    H[0] = H[1]-300*Cp[1]

    
    
#correct what appears to be a mistake in the enthalpy...
idx15300 = flatnonzero(T>=15300)[0]
idx10100 = flatnonzero(T>=10100)[0]

a = 99
b = 151
c = 160


H = fixBump(H,c)
H = correctDiscontinuity(H, a, b)
Cp = fixBump(Cp, c)
Cp = correctDiscontinuity(Cp, a, b)

H-=H[0]

#H[1] = H[3]-Cp[1]*200



#h300 = 277235 #from younglove 1982
#h300 = 573000 #sum of heat of fusion, vaporization, e(84) from Asger, specific heat term
#H[1:]+= h300 - H[1]


Cp[0] = (H[1]-H[0])/300.



with open(os.path.join(os.path.curdir,murphyFolder, transportFile),'r') as f:
    D = f.readlines()[5:]
    X = float64([L.split() for L in D]).T
    MU = extend(X[transportIdxDict['viscosity']], T, Textend)
    
    K = extend(X[transportIdxDict['thermalConductivity']], T, Textend)
    S = extend(X[transportIdxDict['electricalConductivity']], T, Textend)

T = r_[0, T, Textend]
hFun = interp1d(T,H)

K = fixBump(K, c)
K = correctDiscontinuity(K, a, b)



S = fixBump(S, c)
S = correctDiscontinuity(S, a, b)

S+= linspace(1e-6,1, len(S))


if __name__=="__main__":
        
    writeLUT(os.path.join(LUTFolder, 'enthalpy.LUT'), linspace(T[0], T[-1],N), hFun(linspace(T[0], T[-1],N)))
    
    h = linspace(H[0], H[-1], N)    
    Th = interp1d(H,T, bounds_error = False, fill_value = 0)(h)
    writeLUT(os.path.join(LUTFolder, 'temperature.LUT'), h, Th)
    
    P = [Cp, MU, K, S]
    Ph = [interp1d(H,p, bounds_error = False, fill_value = 0)(h) for p in P]
    Pnames = 'specificHeat viscosity thermalConductivity electricalConductivity'.split()
    [writeLUT(os.path.join(LUTFolder, Pnames[i]+'.LUT'), h, ph) for i,ph in enumerate(Ph)]
    
    pDict = dict(zip('enthalpy temperature'.split()+Pnames, [h, Th]+Ph))
    with open('argonPropertiesDict.pkl', 'wb') as f:
        pkl.dump(pDict, f, -1)

    close('all')
    figure()
    subplot(2,3,1)
    plot(T,H)
    xlabel('temperature')
    subplot(2,3,2)
    plot(h, Th)
    ylabel('temperature')
    for i,ph in enumerate(Ph):
        subplot(2,3,3+i)    
        plot(h, ph)
        ylabel(Pnames[i])
    


