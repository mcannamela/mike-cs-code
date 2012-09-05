#this script draws many samples from a gauss-distributed ensemble of particles and forms their images
#the images are inverted to obtain an output ensemble, and the statistics of the input and output are compared

import adHocInverseSentinel as invSens
import particleModelEnsemble as pmEns
import numpy as np
from pylab import *
from plotMacros import *
import mlabMacros as mlm
import time
import os
import scipy
sf = lambda figName_:  savefig(figName_+'.png', format = 'png', transparent = True)

def superHist(InArrs, OutArrs, forceBins = False, N = 100):
    figure(figsize = (20,20))
    cnt = 1
    bins = empty(InArrs.shape, dtype = 'object')
    for idx, A in ndenumerate(InArrs):
        suplot(idx[i], idx[j], cnt)
        [n, bins[idx], patches] = hist(A.ravel(), N,normed = True,alpha = .85, facecolor = 'k', edgecolor = 'w',label = 'input')
        cnt +=1
    cnt = 1
    for idx, A in ndenumerate(OutArrs):
        suplot(idx[i], idx[j], cnt)
        if forceBins:
            hist(A, bins[idx], normed = True,alpha = .45, facecolor = 'r',edgecolor = 'w',label = 'output')
        else:
            hist(A, N, normed = True,alpha = .45, facecolor = 'r',edgecolor = 'w',label = 'output')
        cnt+=1

def inHist(arr, N = 100, bins = None):
    if bins==None:
        [n, bins, patches] = hist(arr.ravel(), N,normed = True,alpha = .85, facecolor = 'k', edgecolor = 'w',label = 'input')
        
    else:
       [n, bins, patches] =  hist(arr.ravel(), bins,normed = True,alpha = .85, facecolor = 'k', edgecolor = 'w',label = 'input')
    return (n, bins, patches)
    
def outHist(arr, N = 100, bins = None, color = 'r', label = 'observed'):
    if bins==None:
        [n, bins, patches] = hist(arr.ravel(), N,normed = True,alpha = .45, facecolor = color,edgecolor = 'w',label = label)
        
    else:
        [n, bins, patches] = hist(arr.ravel(), bins,normed = True,alpha = .45, facecolor =color,edgecolor = 'w',label = label)

    return (n, bins, patches)
#//////////////////////////////    SCRIPT START      \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\    
if __name__=='__main__':
    #simulation tag
    simTag =  'recognitionProb'
    
    figFolder = os.path.join(os.curdir, simTag)
    if not os.path.isdir(figFolder):
        os.mkdir(figFolder)
    close('all')
    ################# NR FRAMES AND NR PARTICLES###################    
    nFrames = 1000
    nPart = 100
    #################################################################
    
    ################### INITIALIZE OBJECTS ##################################
    PE = pmEns.emittingParticleEnsemble(pmEns.marginalRandomEmittingParticleGenerator)
    IS = invSens.adHocInverseSentinel_2dLUT(highPassOn = True, lowPassOn = True)
    
    print IS.FS.gains
    IS.FS.gains*= 0
    IS.FS.gains+= 6*(3e13)
    print IS.FS.gains
    
    IS.highPassConstant = .99
    #######################################################################
    
    ##################  INITIALIZE DATA ARRAYS ###########################
    # input
    Tin = zeros((nFrames,nPart))
    Din = zeros((nFrames,nPart))
    Jin = zeros((nFrames,nPart)) #intensity
    Zin = zeros((nFrames,nPart))
    
    #particle wise output
    Tout = zeros((nFrames,nPart))
    Dout = zeros((nFrames,nPart))
    Jout = zeros((nFrames,nPart)) #intensity

    
    #plume output
    Tobs = array([])
    Dobs = array([])
    Jobs = array([]) 
    
    ###################################################################
    
    ####################    SIMULATE ####################################
    print "simulating and inverting %d frames with %d particles in each"%(nFrames, nPart)
    start = time.time()
    for i in range(nFrames):
        print "completed %d of %d in %d s"%(i,nFrames,time.time()-start )
        
        #generate the ith ensemble
        PE.generate(tStats = [3000, 200],zStats = [0, 4e-3], dStats = [50e-6, 15e-6], N = nPart)
        
        #update input arrays
        Tin[i] = PE.temperature().flatten()
        Din[i] = PE.diameter().flatten()*1e6
        Zin[i] = PE.lateralPosition().flatten()
        
        #invert particle images individually
        IS.particleWiseInvert(PE)
       
       #update input intensities
        for j, P in ndenumerate(IS.FS.sensorOutput['particleImages'].ravel()):
            Jin[i,j] = sum(P[0])
       
       #update particle wise output arrays
        Tout[i] = IS.temperature().flatten()
        Dout[i] = IS.diameter().flatten()*1e6
        Jout[i] = IS.intensity.flatten()
        
        #invert whole plume image
        IS.invert(PE)
        
        print (np.mean(IS.sensorImage[0]), np.max(IS.sensorImage[0]))
        
        #update observed output arrays
        Tobs = r_[atleast_1d(Tobs), IS.temperature()]
        Jobs = r_[atleast_1d(Jobs), IS.intensity]
        Dobs = r_[atleast_1d(Dobs), IS.diameter()*1e6]
    
    ######################################################################
    
    #////////////////////////////     SUMMARY AND FRAME PLOTTING /////////////////////////////////////
    print "detected %d of %d, rate = %f"%(Tobs.ravel().shape[0],Tin.ravel().shape[0], float64(Tobs.ravel().shape[0])/Tin.ravel().shape[0])
#    IS.show()
#    tit(simTag)
#    sf(os.path.join(figFolder, 'filteredFrame_%d_part_%d_frames'%(nPart,nFrames)))
    
    figure(figsize = (20,20))
    plot(IS.FS.sensorOutput['sensorImage'][0])
    ylim(0,256)
    axisFontSize()
    tit(simTag)
    sf(os.path.join(figFolder,'rawFrame_%d_part_%d_frames'%(nPart,nFrames)))
    #////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
   
    #::::::::::::::::::::::::::::::::::::::   HISTOGRAMS ::::::::::::::::::::::::::::::::::
    
    ############ MARGINALS#######################
    figure(figsize = (20,20)); subplot(2,2,1)
    #temprature
    #vvvvvvvvvvvvvvvvvvvvvv
    h, b, p = inHist(      Tin     )
    outHist(      Tout   ,color = 'y', label = 'per particle output')
    outHist(      Tobs   )
    #^^^^^^^^^^^^^^^^^^^
    xlab('temperature, K');          axisFontSize()
    tit(simTag)
    
    #diameter
    subplot(2,2,2)
    #vvvvvvvvvvvvvvvvvvvvvv
    inHist(         Din       )
    outHist(      Dout      , color = 'y', label = 'per particle output'  )
    outHist(      Dobs     )
    #^^^^^^^^^^^^^^^^^^^
    tit('auto binning')
    xlab('diameter, $\mu$m');legend();axisFontSize()
    
    #diameter autocorrelation
    subplot(2,2,3)
     #vvvvvvvvvvvvvvvvvvvvvv
    inHist(         log(Jin)     )
    outHist(      log(Jout),    color = 'y', label = 'per particle output')
    outHist(      log(Jobs)     )
    #^^^^^^^^^^^^^^^^^^^
    xlab('log(total intensity)')
    tit(simTag)
    axisFontSize()
       
    #diameter forced bins
    subplot(2,2,4)
    #vvvvvvvvvvvvvvvvvvvvvv
    [nDout, b, patches]  = inHist(      Din   )
    outHist(      Dout        , bins = b, color = 'y', label = 'per particle output')
    outHist(      Dobs        , bins = b)
    #^^^^^^^^^^^^^^^^^^^
    xlab('diameter, $\mu$m');axisFontSize()
    tit('forced binning')
    sf(os.path.join(figFolder,'T-DHist_%d_part_%d_frames'%(nPart,nFrames)))
    #######################################################
    
    dBins = linspace(10, 500, 400)
	
	
  
    
    