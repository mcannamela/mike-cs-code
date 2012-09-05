import os 
import scipy
scipy.pkgload('io')
from numpy import *
import cPickle as pkl
import time

inputDir = "D:\\00PLASMA_DATA\\juansMatlabData\\frames"
outputDir = "D:\\00PLASMA_DATA\\trellesTorch"
outputNm = "trellesTorch"

print inputDir
print outputDir
print outputNm
print "opening mesh file to determine output size..."
meshNm = 'meshfe_torch5.mat'
M = scipy.io.loadmat(os.path.join(inputDir, os.pardir, meshNm) )
X = M['NODES'].T.copy()
del M
n = X.shape[1]

print "obtaining second output dimension from walk of inputDir"
fnames = os.walk(inputDir).next()[-1]
m = len(fnames)

print "allocating output dictionary"
D = dict(zip(           ['t', 'y'],[zeros(m), zeros((m,n)) ]    ))

componentDict = dict(zip(range(10),['P', 'U1','U2','U3', 'Th','Te', 'E', 'A1','A2','A3'] ))

print "beginning conversion...have patience"
start = time.time()
for i in range(10):
    print "now working on "+componentDict[i]
    for j in range(m):
        if j%20==0:
            print "completed  %d of %d frames for %d of %d components in %d s"%(j, m, i, 10, time.time()-start)
        idx = int32(fnames[j][fnames[j].rfind('_')+1:-4])
        Din = scipy.io.loadmat(os.path.join(inputDir,fnames[j]))
        D['y'][idx] = squeeze(Din['Y'][slice(i,None,10)].copy())
        D['t'][idx] = squeeze(Din['t'].copy())
        del Din
    print "pickling "+ componentDict[i]
    f = open(os.path.join(outputDir, outputNm+'_'+componentDict[i]+'.pkl'), 'wb')
    pkl.dump(D, f,-1)
    f.close()
