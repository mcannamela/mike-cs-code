from ridgeFinder import *
import os
import gc

idx = 0

computeOn = True
forceBuildNetwork = True
nscales = 75
maxFeatureSize = .30
cutoff = 5.5
nKeep = 3000
pctKeep = 1.
saliencyPctKeep = .4

neighFun = neigh3d_124connected


fnames = []
for fn in os.walk('.').next()[2]:
    if fn[-3:]!='TIF' and fn[-3:]!='tif':
        continue
    fnames += [fn]

pklNames = [fn[:-4]+'BF.pkl' for fn in fnames]
netNames = [fn[:-4]+'NW.pkl' for fn in fnames]

print 'the file is '+fnames[idx]

IM = [sum(atleast_3d(flipud(float64(imread(fn)))), axis = 2) for fn in fnames]


for i in range(len(fnames)):
    idx = i
    print 'now working on '+fnames[idx]
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
        
    else:
        if forceBuildNetwork or not os.path.exists(netNames[idx]):
            print "force compute flag not thrown, loading pickle..."
            start = time.time()
            with open(pklNames[idx], 'rb') as f:
                bf = pkl.load(f)
                ss = bf.ss
                im = bf.ss.ss[0]
            print "done. %d"%(time.time()-start)
    
    if forceBuildNetwork or not os.path.exists(netNames[idx]):
        print 'building the network...'
        start = time.time()
        BN = blobNetwork(bf, sizeCutoff= 5.5, 
                             nKeep = nKeep, 
                             pctKeep = pctKeep , 
                             saliencyPctKeep = saliencyPctKeep )
        print "done. %d"%(time.time()-start)
        
        try:
            print 'bouncing to pickle...'
            start = time.time()
            with open(netNames[idx], 'wb') as f:
                pkl.dump(BN, f, -1)
            print "done. %d"%(time.time()-start)
        except:
            print 'pickling failed'
    else:
        print "force build network flag not thrown, loading pickle..."
        start = time.time()
        with open(netNames[idx], 'rb') as f:
            BN = pkl.load(f)
        print "done. %d"%(time.time()-start)
    
    BN.retrim(cutoff)
    BN.resnip((nKeep, pctKeep, saliencyPctKeep))
    
    print int(BN.snidx.shape[0])
    
    subplot(2,3,i)
    hist(BN.meanSig[BN.snidx], 20)
    gc.collect()

( ellipseFig, histFig, scatFig) = BN.plot()

 
show()