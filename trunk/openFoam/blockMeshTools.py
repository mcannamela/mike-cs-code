from numpy import *
import os
import pdb
from pylab import *
import time

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

def uniqueRows(A, atMost = inf):
    U = []
    cnt = 0
    
    for a in A:
        if cnt>atMost:
            break
        isUnique = True
        for u in U:
            if all(a==u):
                isUnique = False
                break
        if isUnique:
            U+=[a.copy()]
            cnt+=1
        
    return array(U)
            

class blockMeshFileReader(object):
    def __call__(self, meshFolder):
        self.fname = self.getName(meshFolder)
        self.readFile()
        self.idx = 0
    
    def getName(self, meshFolder):
        pass 
        
    def readFile(self):
        with open(self.fname,'r') as f:
            self.allLines = f.readlines()
    
    def find(self, s):
        idx = int(self.idx)
        while True:
            #print (self.idx, self.allLines[self.idx], s, self.allLines[self.idx].find(s))
            if self.allLines[self.idx].find(s)!=-1:
                return True
            else:
                if self.idx==(len(self.allLines)-1):
                    self.idx = idx
                    return False
                self.idx+=1
                
        
            
        
class readBlockMeshDict(blockMeshFileReader):
    def __call__(self, meshFolder):
        blockMeshFileReader.__call__(self, meshFolder)

        self.find('convertToMeters')
        self.convertToMeters = double(self.allLines[self.idx].split()[-1][:-1])
        
        self.find('vertices')
        vStart = int(self.idx)+2
        self.find(');')
        vStop = int(self.idx)
        pL = self.allLines[vStart:vStop]
        self.vertices = float64([float64(L.replace('(','').replace(')','').split()) for L in pL])*self.convertToMeters
        
        self.find('hex')
        self.n = int32(self.allLines[self.idx].split(')')[1].replace('(','').split())
        
        xMin = [amin(self.vertices[:,i]) for i in range(3)]
        xMax = [amax(self.vertices[:,i]) for i in range(3)]
        
        found = self.find('boundary')
        if not found:
            found = self.find('patches')
        assert found, "found neither 'boundary', nor 'patches', so cannot get patchnames"
        self.find('(')                
        self.idx+=1
        
        self.patchNames = []
        idx = int(self.idx)
        self.find('mergePatchPairs')
        mergeIdx = int(self.idx)
        self.idx = int(idx)
        
        print 'entering patch reading loop'
        while True:
            if self.find('patch ') or self.find('wall ') or self.find('patch;') or self.find('wall;'):
                if self.idx<mergeIdx:
                    
                    try:
                        self.patchNames += [self.allLines[self.idx-2].split()[1].replace(' ','').replace('\t','').replace('\n','')]
                    except IndexError:
                        self.patchNames += [self.allLines[self.idx-2].split()[0].replace(' ','').replace('\t','').replace('\n','')]
                self.idx+=1
            else:
                break
            self.idx+=1
            
        print 'done, patchnames are '
        print self.patchNames

        return (xMin, xMax, self.n.copy(), list(self.patchNames))
        
    def getName(self, meshFolder):
        return os.path.join(meshFolder, 'blockMeshDict')

    def findNextName(self):
        while True:        
            if self.allLines[self.idx].replace(' ','').replace('\t','').replace('\n','')=='':
                self.idx+=1
                continue
            else:
                return self.allLines[self.idx]

class patch(object):
    def __init__(self, name, nFaces, startFace):
        self.name = str(name)
        self.nFaces = int(nFaces)
        self.startFace = int(startFace)
        
    def faces(self):
        return int32(arange(self.startFace, self.startFace+self.nFaces))
                
class readBoundary(blockMeshFileReader):
    def __call__(self, meshFolder):
        blockMeshFileReader.__call__(self, meshFolder)
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
        
class readPoints(blockMeshFileReader):
    def __call__(self, meshFolder):
        blockMeshFileReader.__call__(self, meshFolder)

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
        return self.points.copy()
    
    def getName(self, meshFolder):
        return os.path.join(meshFolder, 'points')
        
    def getPoints(self):
        self.points = float64([float64(L.replace('(','').replace(')','').split()) for L in self.pL])
        
class readFaces(readPoints):
    def getName(self, meshFolder):
        return os.path.join(meshFolder, 'faces')
    def getPoints(self):
        self.points = int32([int32(L[1:].replace('(','').replace(')','').split()) for L in self.pL])
        
class readOwner(readPoints):
    def getName(self, meshFolder):
        return os.path.join(meshFolder, 'owner')
    def getPoints(self):
        self.points = int32(self.pL)
class readNeighbor(readOwner):
    def getName(self, meshFolder):
        return os.path.join(meshFolder, 'neighbour')

class readField(readPoints):
    def getName(self, meshFolder):
        """
        we pass a whole filename, since in general it could be anywhere
        """
        return meshFolder
    def getPoints(self):
        if len(self.pL[0].replace('(','').replace(')','').split())>1:
            self.points = float64([float64( L.replace('(','').replace(')','').split() ) for L in self.pL]).T
        else:
            self.points = float64(self.pL)
        
def readerArray():
    return [readBoundary(), readBlockMeshDict(), readPoints(), readFaces(), readOwner(), readNeighbor()]
    
class cell(object):
    def __init__(self, idx, BM):
        self.BM = BM
        self.idx = idx
        self.owned = []
        self.points = float64(zeros((8,3)))
    def assignOwned(self, f):
        self.owned+=[f]
    def getPoints(self):
        pIdx = self.BM.pointIndices(self.owned)
        self.points = self.BM.points[pIdx]
        self.center = mean(self.points, axis = 0)
        
class foamFieldWriter(object):
    def __init__(self, cellSubs):
        """
        cellSubs - 3 x nCells array of subscripts into the fields we'll be writing
                s.t. cellSubs[:,i] gives us the indices of the i'th cell
        """
        self.S = cellSubs
    def __call__(self, X, patchNames, fname):
        """
        X - (3 x) nx x ny x nz array holding the field data, could be vector or 
            scalar
        """
        
        self.objName = os.path.split(fname)[1]
        self.patchNames = patchNames
        if X.ndim==3:
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
            for s in self.S.T:
                f.write('%e\n'%X[s[0],s[1],s[2]])
            f.write(');')
            self.writeBoundary(f)
        
    def writeVectorField(self, fname, U):
        with open(fname, 'w') as f:
            self.writeFoamHeader(f, 'volVectorField', self.objName)
            f.write('\ndimensions      [0 0 0 0 0 0 0];\n')
            f.write('\ninternalField   nonuniform List<vector>')
            f.write('\n%d\n(\n'%self.nPoints)
            for s in self.S.T:
                try:
                    f.write(self.vecStr(U[:,s[0],s[1],s[2]]))
                except IndexError:
                    pdb.set_trace()
            f.write(');')
            self.writeBoundary(f)
    def vecStr(self, u):
        return '(%e %e %e)\n'%(u[0], u[1], u[2])    
    def writeBoundary(self, f):
        f.write('boundaryField\n{\n')        
        for p in self.patchNames:
            f.write(p+'\n{\n        type       zeroGradient;\n}\n')
        f.write('\n}')
    def writeFoamHeader(self,f, c, ob):
        """
        write the foam file header with class c and object ob
        to the file object f
        """        
        f.write(foamHeader)
        f.write('\n    class       '+c+';')        
        f.write('\n    object      '+ob+';')
        f.write('\n}\n // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n')
        
    
class blockMesh(object):
    def __init__(self, meshFolder):
        print 'reading raw files'
        patches, block, P, F, own, neigh = tuple([R(meshFolder) for R in readerArray()])
        
        #minimum values for all three dimensions        
        self.xMin = block[0]
       
        #maximum values for all three dimensions
        self.xMax = block[1]
       
        #number of cells in all three dimensions
        self.nx = block[2]
        
        self.patchNames = block[3]
        self.patches = patches
        self.patchNames = [p.name for p in patches]
        self.patchDict = dict(zip(self.patchNames, self.patches))
        
         
        
        #nPoints x 3 array of vectors giving the location of each point
        self.points = P.copy()
        
        #nFaces x 4 array of indices into points, giving 4 points that make each face
        self.faces = F.copy()
        
        #nFaces array of cell labels, telling which cell owns each face 
        #(self.own[i]=n means that cell n owns face i)
        self.own = own.copy()
        self.neigh = neigh.copy()
        
        print 'allocating cells'
        self.cells = [cell(i,self) for i in range(amax(self.own)+1)]
        print 'populating cells'        
        self.populateCells()
        
        self.cc = self.getCoords()
        
        print 'computing cell subscripts'
        self.cellSubs = int32(self.cc2sub(self.cellCenters()))
        
        self.fieldReader = readField()
        self.fieldWriter = foamFieldWriter(self.cellSubs.copy())
        self.fieldNames = []
        self.rawFields = []
        self.fields = []
    
    def getCoords(self):
            try:
                #raw center coords of every cell, 3 x nCells
                cc = self.cellCenters()
                x0 = array(self.xMin)
                xN = array(self.xMax)
                
                #round to 1/100 of the equal spaced increment
                roundFac = self.nx*100/( xN-x0 )
                ccInts = int64((cc-x0[...,newaxis])*roundFac[...,newaxis])
                
                coords = [x0[i]+unique(c)/roundFac[i] for i,c in enumerate(ccInts)]
            except:
                pdb.set_trace()
            return coords
    def cc2sub(self,cc):
        """
        convert cell centers to subscripts
        
        cc - 3 x nCells array of cell centers
        """
        S = zeros_like(cc)
#        for i in range(3):
#            x0 = self.cc[i][0]
#            x1 = self.cc[i][-1]
#            S[i] = int32(.5+(self.nx[i]-1)*(cc[i]-x0)/(x1-x0))

        for i,c in enumerate(cc):
            S[i] = argmin(abs(c[newaxis,...]-self.cc[i][...,newaxis]), axis = 0)

        return S

    def pointIndices(self, F):
        """
        take a list of faces F and convert it to a list of (unique) indices into
        points
        """
        pidx = int32(zeros(len(F)*4))
        for i,f in enumerate(F):
            pidx[4*i:(4*i+4)] = self.faces[f]
        return unique(pidx)
        
        
    def populateCells(self):
        """
        assign the cells their owned faces
        """
        print '    assigning faces to cells'
        for f,odx in enumerate(self.own):
            self.cells[odx].assignOwned(f)
        for f,odx in enumerate(self.neigh):
            self.cells[odx].assignOwned(f)
        print '    finding points for cells'
        start = time.time()
        for i,c in enumerate(self.cells):
            if mod(i,1000)==0:
                print 'completed %d of %d in %.2f min'%(i, len(self.cells), (time.time()-start)/60.)
            c.getPoints()
        #[c.getPoints() for c in self.cells]

    def cellCenters(self):
        return array([c.center.copy() for c in self.cells]).T
    
    def readField(self, fname):
        return self.fieldReader(fname)
    def setField(self,fname):
        nm = os.path.split(fname)[1]
        S = self.cellSubs
        raw = self.readField(fname)
        if nm not in self.fieldNames:
            self.fieldNames+=[nm]
            self.rawFields+=[raw]
            idx = -1
        else:
            idx = flatnonzero(array([nm==n for n in self.fieldNames]))[0]
            self.rawFields[idx] = raw
        
        
        if raw.ndim ==1:
            self.fields += [zeros(tuple(self.nx))]
            self.fields[idx][tuple(S[0]), tuple(S[1]), tuple(S[2])] = raw.copy()

            
        else:
            self.fields += [zeros((3,)+tuple(self.nx))]
            for i in range(3):
                self.fields[idx][i,tuple(S[0]), tuple(S[1]), tuple(S[2])] = raw[i].copy()
        self.fieldDict = dict(zip(self.fieldNames, self.fields))
    def bounceField(self, X, fname):
        """
        write the field X to the file fname in OpenFOAM format
        """
        self.fieldWriter(X,self.patchNames, fname)
        
class blockMeshDistanceFunction(object):
    def __init__(self, bm):
        self.bm = bm
    
    def __call__(self, patchName):
        cc = self.bm.cellCenters()
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
        p = self.bm.patchDict[patchName]
        faceIdx = p.faces()
        pointIdx = self.bm.pointIndices(faceIdx)
        
        return self.bm.points[pointIdx]
        

        
        
if __name__=="__main__":
    meshWriteTest = False
    distanceTest = False    
    if meshWriteTest:
        meshFolder = os.path.join(os.curdir,'exampleBlockMesh','polyMesh')
        #block, P, F, own = tuple([R(meshFolder) for R in readerArray()])
        BM = blockMesh(meshFolder)
        print 'blockmesh is read!'
        BM.setField(os.path.join(os.path.split(meshFolder)[0], 'f'))
        BM.setField(os.path.join(os.path.split(meshFolder)[0], 'U'))
        
        f = BM.fieldDict['f'].copy()
        U = BM.fieldDict['U'].copy()
        
        BM.bounceField(1-f, os.path.join(os.path.split(meshFolder)[0],'0','rf'))
        BM.bounceField(-U,os.path.join(os.path.split(meshFolder)[0], '0','mU'))
        
    if distanceTest:
        runFolder = '/media/raidArray/HVAC_runs/myBuoyantBoussinesqPimpleFoam/hotRoom'
        meshFolder = os.path.join(runFolder, 'constant','polyMesh')
        BM = blockMesh(meshFolder)
        
        DF = blockMeshDistanceFunction(BM)
        D = [DF(pn) for pn in BM.patchNames]
        wallDistance = .3
        isWall = zeros_like(D[0])
        for d in D:
            isWall[d<wallDistance] = 1
            
        BM.bounceField(isWall, os.path.join( runFolder,'isWall'))
