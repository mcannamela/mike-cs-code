import pysces
import pymc
from numpy import *
from numpy.random import randn 
from matplotlib import *
import pdb
from pylab import *
import time
import cPickle as pkl
from scipy import stats
from inversePysces import loadSpeciesDict

def bounceSpecies(S, fname = 'speciesBounce.spc'):
    with open(fname, 'w') as f:
        for i,s in enumerate(S):
            f.write('s%d = %.3e\n'%(i+1, s[-1]))
def spc(pcMod):
    return pcMod.data_sim.getSpecies()[:,1:].T

if __name__=='__main__':
    #-------------------------------------------------------
    #             initialize the pysces model
    #-------------------------------------------------------
    modName = 'R-ERK_s.psc'
    pcMod = pysces.model(modName)
    pcMod.__settings__["lsoda_mxstep"] = 2000
#    pcMod.sim_end = 5e2
#    pcMod.Simulate()
#    bounceSpecies(spc(pcMod))
#    pcMod.__settings__['mode_sim_init'] = 3
    pcMod.sim_end = 1e5
    pcMod.sim_points = 200
    
    idxDic, nmDic, sDic =  loadSpeciesDict(modName[:-4]+'_speciesDictionary.txt')
    #-------------------------------------------------------
    
    #u,sv,vh = linalg.svd(pcMod.nmatrix)
    #figure()
    #plot(sv)
    #title('singular values of N')
    
    
    
    pcMod.Simulate()
    t = pcMod.data_sim.getSpecies()[:,0]
    S = pcMod.data_sim.getSpecies()[:,1:].T
    
    figure()
    
    for i, s in enumerate(S):
        plot(t, log10(s),color = cm.copper(double(i)/S.shape[0]),   label = 's%d'%(i+1))
    
    xlabel('time, s')
    ylabel('log10 species concentrations')
    
    
    figure()
    plot(t, log10(S[0]),'k',   label = 'input')
    plot(t, log10(S[-1]),'k',   label = 'output')
    
    xlabel('time, s')
    ylabel('log10 species concentrations')
    legend()
    
    show()