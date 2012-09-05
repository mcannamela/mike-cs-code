import os
from numpy import *
from scipy.interpolate import Rbf
import pdb 
#import mlabMacros as mlm
from pylab import *
import time
import shutil

foamHeader = r"""/*--------------------------------*- C++ -*----------------------------------*\ 
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.0.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;"""

ptHeaders = ["Point Coordinates:0","Point Coordinates:1","Point Coordinates:2"]
timeHeader = "Time"
uHeaders = ["U:0","U:1","U:2"]
THeader = "Temperature"
HHeader = "H"
def readParaFile(fname):
    with open(fname, 'r') as f:
        L = f.readlines()
    H = L[0].replace('"', '').replace('\n','').split(',')
    X = float64([float64(x.split(',')) for x in L[1:]]).T
    return dict(zip(H, [x for x in X]))

def point(pFileDict):
    return float64([pFileDict[k][0] for k in ptHeaders])
def pTime(pFileDict):
    return pFileDict[timeHeader]
    
    
def readParaFolder(F):    
    return [readParaFile(os.path.join(F,fname)) for fname in os.walk(F).next()[2]]

class paraviewOutputArray(object):
    def __init__(self, folder, downSample = 1):
        self.D = readParaFolder(folder)
        sly = slice(0,None, downSample)
        self.t = pTime(self.D[0])[sly]
        self.p = float64([point(d) for d in self.D]).T
        
        self.T = float64([d[THeader][sly] for d in self.D])
        self.H = float64([d[HHeader][sly] for d in self.D])
        U = [[] for h in uHeaders]
        for i, h in enumerate(uHeaders):
            U[i] = float64([d[h][sly] for d in self.D])
        self.U = float64(U)
        
        
    def x(self):
        return self.p[0].copy()
    def y(self):
        return self.p[1].copy()
    def z(self):
        return self.p[2].copy()
        
class dummyParaviewOutputArray(paraviewOutputArray):
    def __init__(self, p, t, T, H,U, downSample = 1):
        self.t = t
        self.p = p
        self.T = T 
        self.H = H
        self.U = U
        sly = slice(0,None, downSample)
        
class torchOutlet(paraviewOutputArray):
    def setup(self, TFill, HFill, D, N = None ):
        """
        set up interpolants of the outlet data onto square grids
        
        TFill - temperature to fill in where there is no interp data
        HFill - enthalpy to fill in where there is no interp data
        
        D - diameter of the nozzle in m
        N - number of meshpoints to use
        """
        self.TFill = TFill
        self.HFill = HFill
        self.UFill = 0
        if N ==None:
            self.N = int((D**2*(len(self.x())/(pi*self.R()**2)))**.5)
        else:
            self.N = N
        self.G = .5*D*mgrid[-1:1:self.N*1j,-1:1:self.N*1j]
        self.Tg = zeros((len(self.t),)+ self.G[0].shape) 
        self.Hg = zeros((len(self.t),)+ self.G[0].shape) 
        self.Ug0 = zeros((len(self.t),)+ self.G[0].shape) 
        self.Ug1 = zeros((len(self.t),)+ self.G[0].shape) 
        self.Ug2 = zeros((len(self.t),)+ self.G[0].shape) 
        self.Ug = float64([zeros((len(self.t),)+ self.G[0].shape)  for i in range(3)])
        self.interpolate()
        
    def interpolate(self):
        Y = [self.T, self.H] + [u for u in self.U]
        
        I = 'Tg Hg Ug0 Ug1 Ug2'.split()
        fill = [self.TFill, self.HFill, self.UFill, self.UFill, self.UFill]
        start = time.time()
        print 'beginning interpolation'        
        for i,y in enumerate(Y):      
            print 'completed %d of %d in %d s'%(i,len(Y), time.time()-start)
            for j, f in enumerate(y.T):
                rbf = Rbf(self.y(), self.z(), f, function ='cubic')
                self.__dict__[I[i]][j] = self.mask()*rbf(self.G[0], self.G[1]).copy()+(1-self.mask())*fill[i]
        self.Ug[0] = self.Ug0            
        self.Ug[1] = self.Ug1            
        self.Ug[2] = self.Ug2            
        
    def R(self):
        return amax(self.y()**2+self.z()**2)**.5
    def magUg(self):
        return sum(self.Ug**2, axis = 0)**.5
    def mask(self):
        return int32((self.G[0]**2+self.G[1]**2)<self.R()**2)
    def GMin(self):
        return [amin(g) for g in [self.G[0], self.G[1], self.t-1e-10]]
    def GMax(self):
        return [amax(g) for g in [self.G[0], self.G[1], self.t+1e-10]]
    def n(self):
        return [self.G[0].shape[0], self.G[0].shape[1], self.t.shape[0]]
        
class dummyTorchOutlet(torchOutlet):
    def __init__(self, p, t, T, H,U):
        self.t = t
        self.p = p
        self.T = T 
        self.H = H
        self.U = U
        
    

class torchLUTWriter(object): 
    def __call__(self, TOutlet, oFile = "default", zeroTime = True):
        self.tMin = TOutlet.GMin()
        self.tMax = TOutlet.GMax()
        if zeroTime:
            dt = self.tMax[-1]-self.tMin[-1]
            self.tMin[-1] = -1e-10
            self.tMax[-1] = dt
        self.n = TOutlet.n()
        self.T = TOutlet.Tg
        self.H = TOutlet.Hg
        self.U = TOutlet.magUg()
        self.oFile = oFile
        self.write()
    
    def write(self):
        F = [self.oFile+'.T.lut', self.oFile+'.H.lut', self.oFile+'.U.lut']
        XX = [self.T, self.H, self.U]
        
        [self.writeFile(f, XX[i]) for i, f in enumerate(F)]
    
    def writeFile(self, fname, XX):
        with open(fname, 'w') as f:
            f.write(self.tLines())
            for k,X in enumerate(XX):
#                print "k = %d"%k 
                for x in X.T:
                    f.write(self.line(x))  
                    
    def line(self, X):
        s = ''
        for x in X:
            s+='%e '%x
        s+='\n'
        return s
        
    def tLines(self):
        tn = self.tMin
        tx = self.tMax
        n = self.n        
        T = []
        for i in range(len(tn)):
            T += [tn[i], tx[i], n[i]]
             
        return '%e %e %d\n%e %e %d\n%e %e %d\n\n'%tuple(T)
    
        
        

class torchOutletWriter(object):
    def __call__(self,T, BDFolder, xVal = 0.0, noSwirl = True):
        """
        write the torchOutlet object T to the boundaryData folder BDFolder
        """
        self.F = BDFolder
        self.T = T
        self.xVal = double(xVal)
        self.noSwirl = noSwirl
        self.writePoints()
        [self.writeTimestepFolder(i) for i in range(len(self.T.t))]
    def nPoints(self):
        return self.T.G[0].size
    def vecStr(self, u):
        return '(%e %e %e)\n'%(u[0], u[1], u[2])
    def writePoints(self):
        with open(os.path.join(self.F, 'points'), 'w') as f:
            self.writeFoamHeader(f, 'vectorField', 'points')
            f.write('\n%d\n(\n'%self.nPoints())
            for i,y in enumerate(self.T.G[0].flatten()):
                f.write(self.vecStr(r_[self.xVal, y, self.T.G[1].flatten()[i]]))
            f.write(')')
   
    def writeTimestepFolder(self, i):
        tStr = str(self.T.t[i])
        F = os.path.join(self.F,tStr )
        if os.path.isdir(F):
            shutil.rmtree(F)
        os.mkdir(F)
        fnames = 'T H U'.split()
        L = [self.T.Tg[i], self.T.Hg[i], self.T.Ug[:,i].copy()]
        if self.noSwirl:
            L[-1][1]*=0
            L[-1][2]*=0
        for j, fn in enumerate(fnames):
            with open(os.path.join(F, fn), 'w') as f:
                if L[j].ndim == 2:
                    self.writeScalar(f, L[j])
                else:
                    self.writeVector(f, L[j])
                    
    def writeScalar(self, f,X):
        self.writeFoamHeader(f, 'scalarAverageField', 'values')
        f.write('\n//Average\n 0.0\n')
        f.write('\n%d\n(\n'%self.nPoints())
#        f.write('\n(\n')
        for i,x in enumerate(X.flatten()):
            f.write('%e\n'%x)
        f.write(')')
        
    def writeVector(self, f, U):
        self.writeFoamHeader(f, 'vectorAverageField', 'values')
        f.write('\n//Average\n (0 0 0)\n')
        f.write('\n%d\n(\n'%self.nPoints())
#        f.write('\n(\n')
        for i,ux in enumerate(U[0].flatten()):
            f.write(self.vecStr(r_[ux, U[1].flatten()[i], U[2].flatten()[i]]))
        f.write(')')
        
    def writeFoamHeader(self,f, c, ob):
        """
        write the foam file header with class c and object ob
        to the file object f
        """        
        f.write(foamHeader)
        f.write('\n    class       '+c+';')        
        f.write('\n    object      '+ob+';')
        f.write('\n}\n // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n')
        
        
if __name__=="__main__":
    #pfolder = 'example/outlet_data'
    pfolder = '/media/fatMan/taggedOpenFoamRuns/jet_sg100_ArPilot_40_300_highDensity/outlet_data'
    #oFolder = 'example/boundaryData/outlet'
    TO = torchOutlet(pfolder, downSample = 1)
    TO.setup(300, 305042, .0805, N = 200)
#    TOW = torchOutletWriter()
#    TOW(TO, oFolder)

    outFile = '/media/fatMan/taggedOpenFoamRuns/jet_sg100_ArPilot_40_300_highDensity/inlet'
    #outFile = 'default'
    TLW = torchLUTWriter()
    TLW(TO, outFile)
    
    with open(outfile+'.pkl','wb') as f:
        pkl.dump(TO, f, -1)
    
        
    
