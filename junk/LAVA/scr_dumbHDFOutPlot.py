from plotMacros import *
import mlabMacros as mlm
import time
import pdb
import cPickle as pkl
from interpolators import regularize2d

def HDFOutRead(fname, nArrMax = None):
        
    with  open(fname) as f:
        
        try:
            t = double(f.readline())
        except:
            t = 0
        
        #get dimensions of arrays
        nplanes = int(f.readline())
        nrow = int(f.readline())
        ncol = int(f.readline())
        nArr = int(f.readline())
        
        #take all we're given if we don't specify a max
        if nArrMax==None or nArrMax>nArr:
            nArrMax = nArr
        
        #preallocate
        X = zeros(nrow+ncol+nArr*nrow*ncol)
        T = zeros((nArrMax,ncol,nrow))
        cnt = 0
        
        #read into a giant array
        for L in f.readlines():
            for i in range(len(L)/12):
                s = L[(1+12*i):(12*(i+1))]
                try:
                    d = float64(s)
                except ValueError:
                    s = s[:7]+'E'+s[7:]
                    d = float64(s)
                X[cnt] = d    
                cnt+=1
            
    
    #now peel off and reshape individual arrays
    r = X[:nrow].copy(); X = X[nrow:]
    z = X[:ncol].copy(); X = X[ncol:]
    for i in range(nArrMax):
        T[i] = X[:nrow*ncol].copy().reshape(nrow, ncol).T
        X = X[nrow*ncol:]
        
    return (z,r,T,t)




def cookTemperature(x):
    return x
def cookVel(x):
    return x/100.
def cookDist(x):
    return x/100.
def cookDensity(x):
    return x*1000
def cookConductivity(x):
    return x*1e-5
def cookViscosity(x):
    return .1*x
def cookSpeciesDensity(x):
    return x*1000
def cookSpeciesMoleFraction(x):
    return x

labs = ['temperature','axialVelocity','radialVelocity','ArMols', 'HMols',
        'NMols','OMols']
cookList = [cookTemperature, cookVel, cookVel, cookSpeciesMoleFraction, cookSpeciesMoleFraction, 
            cookSpeciesMoleFraction, cookSpeciesMoleFraction]


if __name__=='__main__':

    nArr = 7
    fnameT = "c:\\code\\LAVA\\I=600_V=60_Ar=40_H2=10.dat"
    try:
        z,r,T,t = HDFOutRead(fnameT, nArr)
    except:
        print 'exception hit!'
        fnameT = "e:\\labCODE\\LAVA\\I=600_V=60_Ar=40_H2=10.dat"
        z,r,T,t = HDFOutRead(fnameT, nArr)
    
    T[3:]/=sum(T[3:], axis = 0)
    
    
#    nr = 200
#    nz = 300
#    rMax = r[-1]
#    zMax = z[-1]
#    Treg = zeros((T.shape[0], nz,nr))
#    rReg = linspace(0, rMax, nr)
#    zReg = linspace(0,zMax, nz)
#    g = meshgrid(rReg, zReg)
    n = (300,200)
    nCut = (60,40)
    Treg = zeros((T.shape[0],)+n)
    showAll = False
    for i in range(T.shape[0]):
        #Treg[i] = pol(array([0.,0.]),array([z[-1],r[-1]]), cookList[i](T[i]))(r_['0,3', g[1], g[0]])
        zReg, rReg, Treg[i] = regularize2d(z,r,cookList[i](T[i]), n,nCut)
        
        if showAll:
          #  mlm.surf3(x = z, y = r, z = T[i],axeLab = ['radial, cm','axial, cm',labs[i]], f = mlm.fig(labs[i]))
            mlm.surf3(x = zReg, y = rReg, z = Treg[i],axeLab = ['radial, cm','axial, cm',labs[i]], f = mlm.fig(labs[i]), clf = False)
        elif not showAll and i<3:
         #   mlm.surf3(x = z, y = r, z = T[i],axeLab = ['radial, cm','axial, cm',labs[i]], f = mlm.fig(labs[i]))
            mlm.surf3(x = zReg, y = rReg, z = Treg[i],axeLab = ['radial, cm','axial, cm',labs[i]], f = mlm.fig(labs[i]), clf = False)

  


    cut = 6
    dumplist = [Treg[i] for i in range(5)]+[Treg[5]+Treg[6]]
    with open('e:\\labcode\\pyrticle\\lavaJet.pkl','wb') as f:
        pkl.dump(dict(zip(['radialCoord', 'axialCoord']+labs[:5]+['airMols'],
                           [cookDist(rReg),cookDist(zReg)]+dumplist )),f, -1)

    mlm.ml.show()

    
    