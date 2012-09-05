from numpy import *
from pylab import *
import os 
import pdb
import blockMeshTools as bm
import time
import mlabMacros as mlm

class readBoundary(bm.blockMeshFileReader):
    def __call__(self, meshFolder):
        bm.blockMeshFileReader.__call__(self, meshFolder)
        self.find('// * * * * * *')
        self.idx+=1
        while True:
            try:
                N = int(self.allLines[self.idx])
                break
            except ValueError:
                self.idx+=1
        P = [self.getPatch() for i in range(N)]
        
        return P
            

    def getName(self, meshFolder):
        return os.path.join(meshFolder, 'boundary')
        
    def getPatch(self):
        self.find('{')
        name = self.allLines[self.idx-1].replace(' ','').replace('\t','').replace('\n','')
        self.idx+=2
        nFaces = int(self.allLines[self.idx].split()[-1].replace(';',''))
        startFace = int(self.allLines[self.idx+1].split()[-1].replace(';',''))
        
        return patch(name, nFaces, startFace)
        
class readFaces(bm.readPoints):
    def __call__(self, meshFolder):
        bm.blockMeshFileReader.__call__(self, meshFolder)

        self.find('// * * * * * *')
        self.idx+=1
        while True:
            try:
                N = int(self.allLines[self.idx])
                break
            except ValueError:
                self.idx+=1
        
        self.pL = self.allLines[self.idx+2:self.idx+2+N]
        self.getPoints()        
        return list(self.points)
        
    def getName(self, meshFolder):
        return os.path.join(meshFolder, 'faces')
    def getPoints(self):
        self.points = [int32(L[1:].replace('(','').replace(')','').split()) for L in self.pL]
        
class patch(object):
    def __init__(self, name, nFaces, startFace):
        self.name = str(name)
        self.nFaces = int(nFaces)
        self.startFace = int(startFace)
        
    def faces(self):
        return int32(arange(self.startFace, self.startFace+self.nFaces))
        
def readerArray():
    return [readBoundary(), bm.readPoints(), readFaces(), bm.readOwner(), bm.readNeighbor()]
    
class polyFieldWriter(bm.foamFieldWriter):
    def __init__(self):
        pass
    def __call__(self, X, patchNames, fname):
        """
        X - (3 x) nCells array holding the field data, could be vector or 
            scalar
        """
        
        self.objName = os.path.split(fname)[1]
        self.patchNames = patchNames
        if X.ndim==1:
            self.nPoints = X.size
            self.writeScalarField(fname,X)
        else:
            self.nPoints = X[0].size
            self.writeVectorField(fname,X)
    
    def writeScalarField(self, fname,X):
        with open(fname, 'w') as f:
            self.writeFoamHeader(f, 'volScalarField', self.objName)
            f.write('\ndimensions      [0 0 0 0 0 0 0];\n')
            f.write('\ninternalField   nonuniform List<scalar>')
            f.write('\n%d\n(\n'%self.nPoints)
            for x in X:
                f.write('%e\n'%x)
            f.write(');')
            self.writeBoundary(f)
        
    def writeVectorField(self, fname, U):
        with open(fname, 'w') as f:
            self.writeFoamHeader(f, 'volVectorField', self.objName)
            f.write('\ndimensions      [0 0 0 0 0 0 0];\n')
            f.write('\ninternalField   nonuniform List<vector>')
            f.write('\n%d\n(\n'%self.nPoints)
            for i in range(U.shape[1]):
                try:
                    f.write(self.vecStr(U[:,i]))
                except IndexError:
                    pdb.set_trace()
            f.write(');')
            self.writeBoundary(f)
            
class polyMesh(bm.blockMesh):
    def __init__(self, meshFolder):
        print 'reading raw files'
        patches, P, F, own, neigh = tuple([R(meshFolder) for R in readerArray()])
        
        self.patches = patches
        self.patchNames = [p.name for p in patches]
        self.patchDict = dict(zip(self.patchNames, self.patches))
        
        
        #nPoints x 3 array of vectors giving the location of each point
        self.points = P.copy()
        
        #nFaces x 4 array of indices into points, giving 4 points that make each face
        self.faces = F
        
        #nFaces array of cell labels, telling which cell owns each face 
        #(self.own[i]=n means that cell n owns face i)
        self.own = own.copy()
        self.neigh = neigh.copy()
        
        print 'allocating cells'
        self.cells = [bm.cell(i,self) for i in range(amax(self.own)+1)]
        print 'populating cells'        
        self.populateCells()
                
        self.fieldReader = bm.readField()
        self.fieldWriter = polyFieldWriter()
        self.fieldNames = []
        self.rawFields = []
        self.fields = []
    def pointIndices(self, F):
        """
        take a list of faces F and convert it to a list of (unique) indices into
        points
        """
        pidx = []#int32(zeros(len(F)*4))
        for i,f in enumerate(F):
            pidx+= self.faces[f].tolist()
        return unique(int32(pidx))
        
    def setField(self,fname):
        pass
    def cc2sub(self,cc):
        pass
    
class polyMeshDistanceFunction(object):
    def __init__(self, pm):
        self.pm = pm
    
    def __call__(self, patchName):
        cc = self.pm.cellCenters()
        pp = self.patchPoints(patchName)
        D = zeros(cc.shape[1])
        print '    finding points for cells'
        start = time.time()
        for i,x in enumerate(cc.T):
            if mod(i,1000)==0:
                print 'completed %d of %d in %.2f min'%(i, len(D), (time.time()-start)/60.)
            D[i] = amin( sum((x-pp)**2, axis = 1) )**.5
            
        return D
            
 
    def patchPoints(self, patchName):
        p = self.pm.patchDict[patchName]
        faceIdx = p.faces()
        pointIdx = self.pm.pointIndices(faceIdx)
        
        return self.pm.points[pointIdx]
        
if __name__=="__main__":
#    #meshFolder = os.path.join(os.curdir, 'examplePolyMesh', 'polyMesh')
##    meshFolder = os.path.join('/media/fatMan/taggedOpenFoamRuns/sg100_ArPilot_40_500_sigmaTe','constant', 'polyMesh')
#
#    
#    PM = polyMesh(meshFolder)
#    
#    DF = polyMeshDistanceFunction(PM)
#    
#    D = DF('anode')
#    
#    #PM.bounceField(D, os.path.join(os.curdir, 'examplePolyMesh','cathodeDistance'))
#    PM.bounceField(D, os.path.join( '/media/fatMan/taggedOpenFoamRuns/sg100_ArPilot_40_500_sigmaTe','anodeDistance'))
#    
#    cc = PM.cellCenters()
#    r = (cc[1]**2+cc[2]**2)**.5
#    plot(r, D,'.')
    
    runFolder = '/media/raidArray/HVAC_runs/myBuoyantBoussinesqPimpleFoam/hotRoom'
    meshFolder = os.path.join(runFolder, 'constant','polyMesh')
    PM = polyMesh(meshFolder)    
    DF = polyMeshDistanceFunction(PM)
    D = [DF(pn) for pn in PM.patchNames]
    wallDistance = .3
    isWall = zeros_like(D[0])
    for d in D:
        isWall[d<wallDistance] = 1
        
    PM.bounceField(isWall, os.path.join( runFolder,'isWall'))

    