from numpy import *
from pylab import *

def ofCsvReader(fname):
    with open(fname, 'r') as f:
        L = f.readlines()
    Z = [x.strip().split(',') for x in L]
    X = float64(Z[1:]).T
    
    D = dict()
    
    for i,h in enumerate(Z[0]):
        D[h.replace('"','')] = X[i]
        
    return D
    

def readNetPowerCSV(fname, withTurb = False):
    D = ofCsvReader(fname)
    
    
    P = D['jouleHeating']-abs(D['turbTDiffusion'])*{True:1,False:0}[withTurb]-abs(D['TDiffusion'])-abs(D['Qrad'])
    
    return D,P
    
if __name__=="__main__":
    D,P = readNetPowerCSV('/media/fatMan/taggedOpenFoamRuns/sg100_ArPilot_40_500_highDensity/jouleHeating0.00.csv')
    plot(D['Time'], P,'k.-')