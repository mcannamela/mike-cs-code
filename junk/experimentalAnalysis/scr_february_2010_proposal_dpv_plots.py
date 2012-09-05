from numpy import *
from plotMacros import *

def viscosity(T):
    #return 3.24e-6*exp(2.83e4/T)
    return .0037*exp(6100/T)
    #return 5e-6

def surfaceTension(T):
    #return .65*(1-3.095e-4*(T-4000)/.65)
    return .43 - 3.095e-4*(T-3452)


density = 5819.

a = 1.29
alpha = 1


LS = array([2950, 3025])

fnames = ['Condition1.prt', 'Condition2.prt']
close('all')
for fname in fnames:
    with open(fname, 'r') as f:
        f.readline()
        X= float64([float64(L.split()[4:7]) for L in f.readlines()]).T
        
    V = X[0]
    T = X[1]+273
    D = X[2]*1e-6
    
    ht, bt = histogram(T,150)
    tMelt = bt[argmax(ht)]
    print 'T_melt'
    print tMelt
    
    midx = logical_and(T>tMelt, V>50)
    
    fig()
    hist(T, 150)
    xlab('T')
    
    T = T[midx]
    D = D[midx]
    V = V[midx]
    
    mu = viscosity(T)
    gamma = surfaceTension(T)
    
    Re = density*V*D/(mu)
    We = density*V**2*D/gamma
    
    #splatThickness = (2*D/(3*(a*Re**.2)**2))*1e6
    splatThickness = (D/(.89*Re**.2))*1e6
    N = alpha*We**.5*Re**.25
    sBins = linspace(1,25,150)
    nBins = linspace(0,2200,150)
    
    fig()
    hist(splatThickness, sBins, facecolor = 'k')
    xlim(sBins[0], sBins[-1])
    xlab('Splat Thickness, microns')
    ylab('count')
    tit('Splat Thickness, number weighted')
    axisFontSize(fsz = 18)
    sf('Splat Thickness, number weighted'+fname[:-4])
    
#    fig()
#    hist(N, nBins, facecolor = 'k')
#    xlim(nBins[0], nBins[-1])
#    xlab('Splashing Index')
#    ylab('count')
#    tit('Splashing Index, number weighted')
#    axisFontSize(fsz = 18)
#    sf('Splashing Index, number weighted'+fname[:-4])
    
    hs, bs = histogram(splatThickness, sBins, weights = D**3/sum(D**3))
    hn, bn = histogram(N, nBins, weights = D**3/sum(D**3))
    
    fig()
    bar(bs[:-1], hs, width = diff(bs[:2])[0], facecolor = 'k')
    xlim(sBins[0], sBins[-1])
    xlab('Splat Thickness, microns')
    ylab('volume fraction')
    tit('Splat Thickness, volume weighted')
    axisFontSize(fsz = 18)
    sf('Splat Thickness, volume weighted'+fname[:-4])
    
#    fig()
#    bar(bn[:-1], hn, width = diff(bn[:2])[0], facecolor = 'k')
#    xlim(nBins[0], nBins[-1])
#    xlab('Splashing Index')
#    ylab('volume fraction')
#    tit('Splashing Index , volume weighted')
#    axisFontSize(fsz = 18)
#    sf('Splashing Index , volume weighted'+fname[:-4])