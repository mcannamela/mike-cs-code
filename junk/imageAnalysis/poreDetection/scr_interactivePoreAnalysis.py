from ridgeFinder import *
import os
import gc

#['meanStdDev', 'meanBrightness', 
#                            'minBrightness','medianBrightness',
#                            'nBlobPoints','nFootprintPoints',
#                            'meanCurvature']



##enable the following sliders
#sliderKeys = ['meanStdDev', 'meanBrightness', 
#              'minBrightness','medianBrightness',
#              'meanCurvature']
#
##whether slider defines upper or lower bound 
##-1 for upper bound, +1 for lower bound
#sliderSense = [-1, -1, 
#                -1, -1, 
#                1]
#                
#enable the following sliders
sliderKeys = ['meanStdDev', 
              'medianBrightness',
              'meanCurvature']

#whether slider defines upper or lower bound 
#-1 for upper bound, +1 for lower bound
sliderSense = [-1, 
                -1,
                1]

#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#||||||||||||||make your choices here! |||||||||||||||||||||||
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#pick the file you want to work on                
idx = 6

#force fresh computation of blob points?
computeOn = True

#force fresh computation of bob network?
forceBuildNetwork = True

#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

#number of scales to use in scale space
nscales = 75
#largest feature we'd like to detect in terms of fraction  of image size
maxFeatureSize = .30

#function that defines connectivity in scale space
neighFun = neigh3d_124connected

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::build up file names:::::::::::::::::::
fnames = []
for fn in os.walk('.').next()[2]:
    if fn[-3:]!='TIF' and fn[-3:]!='tif':
        continue
    fnames += [fn]

pklNames = [fn[:-4]+'BF.pkl' for fn in fnames]
netNames = [fn[:-4]+'INW.pkl' for fn in fnames]
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::
print 'the file is '+fnames[idx]

#load the images
IM = [sum(atleast_3d(flipud(float64(imread(fn)))), axis = 2) for fn in fnames]


print 'now working on '+fnames[idx]

#//////////////////////////////////////////////////////
#////////////////////// compute blobFinder ////////////
if computeOn or not os.path.exists(pklNames[idx]):
    
    print "initializing scale space..."
    start = time.time()
    ss = scaleSpace2d(IM[idx], nscales = nscales, maxFeatureSize = maxFeatureSize)
    print "done. %d"%(time.time()-start)
    
    print "initializing blob finder..."
    start = time.time()
    bf = blobFinder2d(ss)
    print "done. %d"%(time.time()-start)
    
    print "finding blobs..."
    start = time.time()
    r = bf()
    print "done. %d"%(time.time()-start)
    
    print 'bouncing to pickle...'
    start = time.time()
    with open(pklNames[idx], 'wb') as f:
        pkl.dump(bf, f, -1)
    print "done. %d"%(time.time()-start)
#//////////////////////////////////////////////////////    

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#\\\\\\\\\\\\\\  load blob finder \\\\\\\\\\\\\\\\\\\\\\
else:
    if forceBuildNetwork or not os.path.exists(netNames[idx]):
        print "force compute flag not thrown, loading pickle..."
        start = time.time()
        with open(pklNames[idx], 'rb') as f:
            bf = pkl.load(f)
            ss = bf.ss
            im = bf.ss.ss[0]
        print "done. %d"%(time.time()-start)
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

#--------------------------------------------------------
#----------------build the blob network------------------
if forceBuildNetwork or not os.path.exists(netNames[idx]):
    print 'building the network...'
    start = time.time()
    BN = interactiveBlobNetwork(bf,neighFun, fname = fnames[idx])
    print "done. %d"%(time.time()-start)
    
    try:
        print 'bouncing to pickle...'
        start = time.time()
        with open(netNames[idx], 'wb') as f:
            pkl.dump(BN, f, -1)
        print "done. %d"%(time.time()-start)
    except:
        print 'pickling failed'
#--------------------------------------------------------

#========================================================
#=====================load a blobNetwork=================
else:
    print "force build network flag not thrown, loading pickle..."
    start = time.time()
    with open(netNames[idx], 'rb') as f:
        BN = pkl.load(f)
    print "done. %d"%(time.time()-start)
#========================================================    

BN(sliderKeys, sliderSense)