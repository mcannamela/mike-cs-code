from numpy import *
from plotMacros import *
import os
import pdb
import mlabMacros as mlm

delaunayPath = 'E:\\labCODE\\qhull-2010.1'
if not os.path.isdir(delaunayPath):
    delaunayPath = 'D:\\CODE\\qhull-2010.1'
    
tmpPath = 'E:\\labCODE\\helpers'
if not os.path.isdir(tmpPath):
    tmpPath = 'D:\\CODE\\helpers'

def delaunayInput(X, fname = 'delaunayInput.tmp'):
    """
    prepare input for qdelaunay taking X.shape[0] as the dimension of the data
    """
    IDX = []
    with open(os.path.join(tmpPath,'delaunayInput.tmp'), 'w') as f:
        f.write('%d\n%d\n'%(X.shape[0], X[0].size))
        for idx, x in ndenumerate(X[0]):
            f.write(str(X[(slice(None),)+idx]).replace('[','').replace(']',''))
            f.write('\n')
            IDX+= [idx]
    tup = tuple([])
    IDX = array(IDX).T
    for i in range(IDX.shape[0]):
        tup+= (IDX[i],)
    return (fname, tup)

def delaunayOutput(fname):
    """
    open fname and read the output of delaunay triangulation
    """
    with open(fname, 'r') as f:
        n = int(f.readline())
        tri = int32([int32(L.split()) for L in f.readlines()])
    return tri

def delaunay(X, fout = 'delaunayOut.tmp'):
    """
    compute the delaunay triangulation of X taking as the dimension X.shape[0]
    """
    fname, idx = delaunayInput(X)
    fin = os.path.join(tmpPath, fname)
    comm = os.path.join(delaunayPath, 'qdelaunay')+' Qt i '
    
    ret = os.system(comm+ 'TI ' +fin +'> '+ os.path.join(tmpPath,fout))
    assert ret==0, 'delaunay did not return 0! check delaunayPath...currently is '+delaunayPath
    tri = delaunayOutput(os.path.join(tmpPath, fout))
    Y = X[(slice(None),)+idx]
    return (tri, idx, Y)

if __name__=='__main__':
    X = rand(3,10,10)
    tri, idx, Y = delaunay(X)
    f = mlm.fig('the points')
    mlm.ml.points3d(X[0], X[1], X[2], figure = f)
    
    mlm.ml.triangular_mesh(Y[0], Y[1], Y[2], tri[:,:3], figure = f)
    mlm.ml.show()
    
    
    
    
    
    
    