from numpy import *
from pylab import *
import matplotlib as mpl
import cPickle as pkl
import pdb

def catRuns(fn, gn):
    with open(fn, 'rb') as f:
        with open(gn, 'rb') as g:
            d1 = pkl.load(f)
            d2 = pkl.load(g)
    d = {}        
    for k in d1.keys():
        if k=='diamteter':
            K='diameter'
        else:
            K = k
        if type(d1[k])==type(zeros(1)):
            d[K] = r_['0,1', d1[k], d2[k]]

        elif type(d1[k])==type(list()):
            d[K] = d1[k]+d2[k]
    return d
    
def meanTemperature(d):
    return sum(d['finalRadius']**3*d['finalTemperature'], axis = -1)/sum(d['finalRadius']**3, axis = -1)
            

if __name__=="__main__":
    
    needsCat = False
    
    if needsCat:
        fn = "timeVaryingDistributionStudy_7.000000_2000.pkl"
        gn = "timeVaryingDistributionStudy_7.000000_1500.pkl"
        D7 = catRuns(fn,gn)
        with open("timeVaryingDistributionStudy_7_3500.pkl", 'wb') as f:
            pkl.dump(D7, f, -1)
    #    for k in D.keys():
    #        if type(D[k])==type(zeros(1)):
    #            print D[k].shape
    #        else:
    #            print len(D[k])
    
        fn = "timeVaryingDistributionStudy_5.000000_2000.pkl"
        gn = "timeVaryingDistributionStudy_5.000000_1500.pkl"
        D5 = catRuns(fn,gn)
        with open("timeVaryingDistributionStudy_5_3500.pkl", 'wb') as f:
            pkl.dump(D5, f, -1)
    
    else:
        with open("timeVaryingDistributionStudy_7_3500.pkl", 'rb') as f:
            D7 = pkl.load(f)
        with open("timeVaryingDistributionStudy_5_3500.pkl", 'rb') as f:
            D5 = pkl.load(f)
    
    
    D = [D5, D7]
    C = 'k r'.split()
    M = 'o x'.split()
    L = '5slm 7slm'.split()
    nHist = 200
    figure()
    for i,d in enumerate(D):
        hist(d['diameter']*1e6, nHist, facecolor = C[i], alpha = .4, label = L[i])
    xlabel('diameter')
    legend()
    savefig('diameterHist.png')
    
    
    figure()
    for i,d in enumerate(D):
        hist(d['initialVelocity'], nHist, facecolor = C[i], alpha = .4, label = L[i])
    xlabel('initial velocity, m/s')
    legend()
    savefig('v0Hist.png')
    
    figure()
    for i,d in enumerate(D):
        hist(d['finalVelocity'], nHist, facecolor = C[i], alpha = .4, label = L[i])
    xlabel('final velocity, m/s')
    legend()
    savefig('vFinalHist.png')     
    
    figure()
    for i,d in enumerate(D):
        hist(d['finalTemperature'][:,-1], nHist, facecolor = C[i], alpha = .4, label = L[i])
    xlabel('final surface temperature, K')
    legend()
    savefig('TsFinalHist.png')
    
    figure()
    for i,d in enumerate(D):
        hist(meanTemperature(d), nHist, facecolor = C[i], alpha = .4, label = L[i])
    xlabel('final mean temperature, K')
    legend()
    savefig('TbarFinalHist.png')
    
    figure()
    for i,d in enumerate(D):
        plot(d['diameter']*1e6,meanTemperature(d), C[i]+M[i], alpha = .3, label = L[i] )
    xlabel('diameter, microns')
    ylabel('temperature, K')
    legend()
    savefig('D-T_scatter.png')
    
    figure()
    for i,d in enumerate(D):
        plot(d['diameter']*1e6,d['finalVelocity'], C[i]+M[i], alpha = .3, label = L[i] )
    xlabel('diameter, microns')
    ylabel('final velocity, m/s')
    legend()
    savefig('D-V_scatter.png')
    
    figure()
    for i,d in enumerate(D):
        plot(d['finalVelocity'],meanTemperature(d), C[i]+M[i], alpha = .3, label = L[i] )
    xlabel('final velocity, m/s')
    ylabel('temperature, K')
    legend()
    savefig('V-T_scatter.png')

    