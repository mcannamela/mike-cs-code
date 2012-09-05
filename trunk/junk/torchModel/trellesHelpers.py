import os 
import scipy
scipy.pkgload('io')
from numpy import *
import numpy as np
import cPickle as pkl
from pylab import *
#import mlabMacros as mlm
import pdb


faces = array([[0,1,2,3], 
                         [0,1,5,4],
                            [0,3,7,4],
                       [2,3,7,6],
                     [1,2,6,5],
                     [4,5,6,7],
                     ])
varDict = dict(zip(['magneticPotential', 
                    'electricPotential', 
                    'pressure', 
                    'electronTemperature', 
                    'heavyTemperature', 
                    'velocity'],
                    
                    [('A1, A2, A3'), 
                     'E',
                     'P', 
                     'Te',
                     'Th',
                     ('U1', 'U2', 'U3')]))
                     
def fetchMesh(outputDir = "D:\\00PLASMA_DATA\\juansMatlabData", elements=False, boundaries = False, precision = 'single'):
    D = scipy.io.loadmat(os.path.join(outputDir,'MESHFE_torch5.mat'))
    N = float32(D['NODES'].copy().T)
    E = None
    bNodes = None
    bFaces = None
    if elements:
        E = D['ELEMS'].copy()-1
    if boundaries:
        bNodes = []
        bFaces = []
        n = D['BNODES'].T  #node nr, boundary nr pairs
        f = D['BFACES'].T   #element nr, face nr, boundary nr
        
        for i in range(1,np.max(n[1])+1):
            bNodes += [n[0][n[1]==i].copy()-1]
            bFaces += [f[:2][:,f[2]==i].copy()-1]
    if precision == 'single':
        return (float32(N),E,bNodes,bFaces)
    else:
        return (float64(N),E,bNodes,bFaces)
    
    
    
def fetchVar(VAR, N = 500, outputDir = "D:\\00PLASMA_DATA\\trellesTorch", precision = 'single'):
    tName = 'trellesTorch_'
    with open(os.path.join(outputDir,tName+VAR+'.pkl'),'rb') as f:
        D = pkl.load(f)
    y = D['y'][:N].copy()
    t = D['t'][:N].copy()
    if precision == 'single':
        return (t,float32(y))
    else:
        return (t,float64(y))
    
def fetchSVD(VAR, N = 500, outputDir = "D:\\00PLASMA_DATA\\trellesTorch", compute = False):
    tName = 'trellesTorch_'
    if compute or not os.path.isfile(os.path.join(outputDir,tName+VAR+'_svd_%d.pkl'%N)):
        t,y = fetchVar(VAR, N, outputDir)
        yBar = np.mean(y,axis = 0)
        y-= yBar
        U,s,Vh = linalg.svd(y, full_matrices = False)
        with open(os.path.join(outputDir,tName+VAR+'_svd_%d.pkl'%N),'wb') as f:
            pkl.dump(dict(zip(['U','s','Vh','columnMeans'],[U,s,Vh,yBar])), f, -1)
    else:
        with open(os.path.join(outputDir,tName+VAR+'_svd_%d.pkl'%N),'rb') as f:
            D = pkl.load(f)
            U,s,Vh,yBar = map(float32, (D['U'].copy(),D['s'].copy(),D['Vh'].copy(),D['columnMeans'].copy()))
    return U,s,Vh,yBar

if __name__=='__main__':
    print 'fetching mesh...'
    x,el,bx,bel = fetchMesh(elements = True, boundaries = True)
    print 'done. Beginning SVD computation for all variables!'
    vars = ['A1', 'A2', 'A3', 'E', 'P', 'Te', 'Th', 'U1','U2','U3']
    
    # close('all')
    # fig()
    # for V in vars:
        # print 'Now computing SVD of '+V
        # U,s,Vh,yBar = fetchSVD(V,N = 200, compute = False)
        # print 'Complete!'
        # plot(20*log10(s/s[0]), '-',linewidth = 1.5, label = V)
    # ylim(-100, 1)
    # legend()
    # xlab('rank')
    # ylab('singular values s$_i$, dB ref s$_0$')
    # tit('summary of singular values for all variables')
    # axisFontSize()
    # sf('SingularValues')
    t,y = fetchVar('U3',N = 200)
    U,s,Vh,yBar = fetchSVD('U3',N = 200, compute = False)
    
    n = 3
    A = array(matrix(U[:n,:n])*matrix(diag(s[:n]))*matrix(Vh[:n]))
    idx = (x[2]==np.max(x[2]))
    # for i in range(4):
        # mlm.scat(r_['0,2',x[0][idx],x[1][idx],w[i]*Vh[i][idx],w[i]*Vh[i][idx]], f = mlm.fig('mode %d'%i))
    
#    mlm.scat(r_['0,2',x[0][idx],x[1][idx],y[0][idx],y[0][idx]], f = mlm.fig('full solution'))
#    mlm.scat(r_['0,2',x[0][idx],x[1][idx],yBar[idx]-y[0][idx],yBar[idx]-y[0][idx]], f = mlm.fig('mean solution error'))
#    mlm.scat(r_['0,2',x[0][idx],x[1][idx],yBar[idx]+A[0][idx]-y[0][idx],yBar[idx]+A[0][idx]-y[0][idx]], f = mlm.fig('rank %d approximation error'%n))
#    mlm.ml.show()
    
    
    