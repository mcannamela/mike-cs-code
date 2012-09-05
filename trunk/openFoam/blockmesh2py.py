from numpy import *
import os
import sys
import pdb
from helpers import ngrid, x2gamma
from pylab import *

#case directory passed as arg or using current directory
if len(sys.argv)>1:
    foamDir = sys.argv[1]
else:
    foamDir = os.curdir

#function to read in the mesh cell locations
def readMesh(fdir):
    meshdir = os.path.join(fdir, 'constant', 'polyMesh')
    black = ['blockMeshDict', 'boundary']
    K = []
    X = []
    
    for fname in  os.walk(meshdir).next()[2]:
        if fname not in black:
            K+=[fname]
            with open(os.path.join(meshdir, fname),'r') as f:
                X+=[readFoamFile(f)]
    rawMeshDict = dict(zip(K, X))
    D = rawMeshDict
    
    D['owner'] = int32(D['owner'])
    D['neighbour'] = int32(D['neighbour'])
    D['faces'] = int32(D['faces'])
    
    nCells = int32(max([amax(D['owner']), 
                        amax(D['neighbour'])]))+1
    
    
    cellFaces = [ int32([]) for i in range(nCells) ]
    
    #pdb.set_trace()
    for i in range(D['owner'].shape[0]):
        cellFaces[D['owner'][i]] = r_[cellFaces[D['owner'][i]], i]
        
    for i in range(D['neighbour'].shape[0]):
        cellFaces[D['neighbour'][i]] = r_[cellFaces[D['neighbour'][i]], i]
    
    for c in cellFaces:
        c =  unique(c)

    cellPts = int32(zeros((nCells, 8)))
    for i in range(nCells):
        cellPts[i] = unique(D['faces'][cellFaces[i]])
    
    cells = zeros((nCells, 3))
    #pdb.set_trace()
    for i in range(nCells):
        cells[i] = mean(D['points'][cellPts[i], :], axis = 0)
    
    return cells

    
def readTimestepDir(tdir, black = ['']):
    K = []
    X = []
    for fname in  os.walk(tdir).next()[2]:
        if fname not in black:
            K+=[fname]
            with open(os.path.join(tdir, fname),'r') as f:
                #print 'reading '+fname
                X+=[readFoamFile(f)]
                
    return dict(zip(K, X))
    
def readCase(fdir, black = ['x2sigma']):
    cells = readMesh(fdir)
    tList = []
    D = []
    black = ['0', 'constant', 'system']
    for dn in os.walk(fdir).next()[1]:
        dirName = os.path.join(fdir, dn)
        if dn not in black:
            try:
                tList+=[float64(dn)]
            except:
                continue
            print 'reading t = '+dn
            D += [readTimestepDir(dirName, black = black)]
            
            

    
    t = float64(tList)
    tidx = argsort(t)
    t = sort(t)
    
    
    X = dict(zip(D[0].keys(),
                   [zeros(t.shape+D[0][k].shape) for k in D[0].keys()]))
    
    for i in range(t.shape[0]):
        for k in D[0].keys():
            X[k][tidx[i]] = D[i][k].copy()
    
    for k in D[0].keys():
        if k in ['p', 'U']:        
            print k        
            for i in range(t.shape[0]):      
                if i>0:
                    print sum((X[k][i]-X[k][i-1])!=0)        
    
            
    return dict(zip(['x', 't', 'u'], [cells, t, X]))
    
def gridify(D):
    """
    D - a foam case dictionary produced by readCase
    """
    X = D['x'].T
    x = [unique(X[i]) for i in range(3)]
    
    nt = D['t'].shape[0]
    G = ngrid(x)
    
    gamma = x2gamma(X.T)    
    Ggamma = x2gamma(r_['0,2', G[0].flatten(), 
                        G[1].flatten(), G[2].flatten()].T)
    
    
    idx = argsort(gamma)
    gidx = argsort(Ggamma)
    
    

    
#    pdb.set_trace()
    
    U = []
    u = D['u']
    for k in u.keys():     
        if u[k].ndim>2:
            n = u[k].shape[-1]
        else:
            n=1
        U += [zeros( (nt,n)+shape(G[0]) )]
        for i in range(n):
            for j in range(nt):
                try:
                    U[-1][j][i].ravel()[gidx] = u[k][j][idx,i]#.reshape(G[0].shape)
                except:
                    U[-1][j][i].ravel()[gidx] = u[k][j][idx]#.reshape(G[0].shape)
                    
    
#        print k
#        print sum(U[-1]!=0)
#        print sum(u[k]!=0)

    #pdb.set_trace()
    DG = dict(zip(['X','G','U'], [x,G,dict(zip(u.keys(),U))]))
    return DG
        

def readFoamFile(f): 
    #burn header
    #print '\t burning header'
    while True: 
        
        if f.readline()[:6]=='// * *':
            f.readline()
            f.readline()
            break
    
    #determine number of entries  
    #print '\t finding number of entries'
    while True:
        theLine = f.readline()
        try:
            nPts = int(theLine)
            break
        except:
            continue
    
    #burn a line
    f.readline()
    
    #determine number of dimensions
    #print '\t finding nDim'
    line = f.readline()
    if line.find('(')!=-1:
        sly = slice(line.find('(')+1, line.find(')')-len(line))
        nDim = len(line[sly].split())
    else:
        nDim = 1
    
    #print '\t\t nDim = %d'%nDim
    #read in the entries according to dimension
    #print '\tparsing entries'
    
    if nDim==1:
        L = [line]+f.readlines()[:nPts-1]
        pts = float64(L)
    else:
        pts = zeros((nPts, nDim))
        L = [line]+f.readlines()[:nPts-1]
        for i in range(nPts):
            pts[i] = float64(L[i][sly].split())
                
    return pts
    
if __name__=='__main__':
    #foamDir = '/home/wichtelwesen/openFoam/tutorials/compressible/rhoPisoFoam/ras/cavity'
    #foamDir = '/home/wichtelwesen/openFoam/tutorials/compressible/rhoPisoFoam/ras/gradedDuct_lowRe'    
    foamDir = '/home/wichtelwesen/Dropbox/bigProposal_spring_2011/hartmann3d_loMag'
    case = readCase(foamDir)
    x = case['x']
    u = case['u']['U']
    t = case['t']
    
    D = gridify(case)
    
    U = squeeze(D['U']['U'][:,0])
    p = squeeze(D['U']['p'])
    G = squeeze(D['G'])
    X = D['X']
    
    
    
    
    #Umag = sum(U**2, axis = 0)**.5
#    quiver(squeeze(D['G'][0]), squeeze(D['G'][1]), U[0], U[1])


    