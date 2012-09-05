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
    
def outHist(arr, N = 100, bins = None):
    if bins==None:
        [n, bins, patches] = hist(arr.ravel(), N,normed = True,alpha = .45, facecolor = 'r',edgecolor = 'w',label = 'output')
        
    else:
        [n, bins, patches] = hist(arr.ravel(), bins,normed = True,alpha = .45, facecolor = 'r',edgecolor = 'w',label = 'output')

    return (n, bins, patches)
#//////////////////////////////    SCRIPT START      \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\    
if __name__=='__main__':
    #simulation tag
    simTag =  'widePlume'
    
    figFolder = os.path.join(os.curdir, simTag)
    if not os.path.isdir(figFolder):
        os.mkdir(figFolder)
    close('all')
    ################# NR FRAMES AND NR PARTICLES###################    
    nFrames = 100
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
    Vin = zeros((nFrames,nPart))
    Din = zeros((nFrames,nPart))
    Jin = zeros((nFrames,nPart)) #intensity
    JmaxIn = zeros((nFrames,nPart))
    
    #output
    Tout = array([])
    ToutPert = array([]) #dithering array for temperature
    Dout = array([])
    Jout = array([]) 
    JmaxOut = array([])
    ###################################################################
    
    ####################    SIMULATE ####################################
    print "simulating and inverting %d frames with %d particles in each"%(nFrames, nPart)
    start = time.time()
    for i in range(nFrames):
        print "completed %d of %d in %d s"%(i,nFrames,time.time()-start )
        
        #generate the ith ensemble
        PE.generate(tStats = [3000, 200],zStats = [0, 4e-3], dStats = [50e-6, 15e-6], N = nPart)
        
        #update input arrays
        Vin[i] = PE.parameterVectors[...,6].flatten()
        Tin[i] = PE.parameterVectors[...,2].flatten()
        Din[i] = PE.parameterVectors[...,1].flatten()*1e6
       
       # print sort(PE.parameterVectors[...,1].flatten()*1e6)
       
       #do the inversion
        IS.invert(PE)
        
        #update input intensities
        for j, P in ndenumerate(IS.FS.sensorOutput['particleImages'].ravel()):
            Jin[i,j] = sum(P[0])
            JmaxIn[i,j] = np.max(P[0])
        
        print (np.mean(IS.sensorImage[0]), np.max(IS.sensorImage[0]))
        
        #update output arrays
        Tout = r_[atleast_1d(Tout), IS.temperature()]
        ToutPert = r_[atleast_1d(ToutPert), (21*random(IS.temperature().shape)- 10.5)]
        Jout = r_[atleast_1d(Jout), IS.intensity]
        JmaxOut = r_[atleast_1d(JmaxOut), IS.maxIntensity]
        Dout = r_[atleast_1d(Dout), IS.diameter()*1e6]
    
    #compute velocity proxy numerator
    vProxyIn = Din**2*Tin**4
    vProxyOut = Dout**2*Tout**4
    ######################################################################
    
    #////////////////////////////     SUMMARY AND FRAME PLOTTING /////////////////////////////////////
    print "mean, std T in: %d , %d"%(mean(Tin), std(Tin))
    print "mean, std T out: %d , %d"%(mean(Tout), std(Tout))
    print "mean, std  D in: %d , %d"%(mean(Din), std(Din))
    print "mean, std D out: %d , %d"%(mean(Dout),std(Dout))
    print "detected %d of %d, rate = %f"%(Tout.ravel().shape[0],Tin.ravel().shape[0], float64(Tout.ravel().shape[0])/Tin.ravel().shape[0])
    IS.show()
    tit(simTag)
    sf(os.path.join(figFolder, 'filteredFrame_%d_part_%d_frames'%(nPart,nFrames)))
    
    figure(figsize = (20,20))
    plot(IS.FS.sensorOutput['sensorImage'][0])
    ylim(0,256)
    axisFontSize()
    tit(simTag)
    sf(os.path.join(figFolder,'rawFrame_%d_part_%d_frames'%(nPart,nFrames)))
    #////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
   
    #::::::::::::::::::::::::::::::::::::::   HISTOGRAMS ::::::::::::::::::::::::::::::::::
    
    ############ TEMP, DIAMETER#######################
    figure(figsize = (20,20)); subplot(2,2,1)
    #temprature
    #vvvvvvvvvvvvvvvvvvvvvv
    inHist(      Tin     )
    outHist(      Tout   )
    #^^^^^^^^^^^^^^^^^^^
    xlab('temperature, K');    legend();        axisFontSize()
    tit(simTag)
    
    #diameter
    subplot(2,2,2)
    #vvvvvvvvvvvvvvvvvvvvvv
    [nDin, b, patches] = inHist(      Din   )
    [nDout, b, patches] = outHist(      Dout )
    #^^^^^^^^^^^^^^^^^^^
    tit('auto binning')
    xlab('diameter, $\mu$m');legend();axisFontSize()
    
    #diameter autocorrelation
    subplot(2,2,3)
    inin = scipy.signal.correlate(float64(nDin), float64(nDin)); inin/=sum(inin)
    outout = scipy.signal.correlate(float64(nDout), float64(nDout)); outout/= sum(outout)
    #vvvvvvvvvvvvvvvvvvvvvv
    plot(inin, color = 'k', linewidth = 2, label = 'input')
    plot(outout, color = 'r', linewidth = 2, label = 'output')
    #^^^^^^^^^^^^^^^^^^^
    
    xlab('autocorrelation of diameter histogram');legend(); axisFontSize()
       
    #diameter forced bins
    subplot(2,2,4)
    #vvvvvvvvvvvvvvvvvvvvvv
    [nDout, b, patches]  = inHist(      Din   )
    outHist(      Dout , bins = b)
    #^^^^^^^^^^^^^^^^^^^
    xlab('diameter, $\mu$m');legend();axisFontSize()
    tit('forced binning')
    sf(os.path.join(figFolder,'T-DHist_%d_part_%d_frames'%(nPart,nFrames)))
    #######################################################
    
    #################### LOG INTENSITY #####################
    figure(figsize = (20,20))
    subplot(1,2,1)
    #vvvvvvvvvvvvvvvvvvvvvv
    inHist(      log(Jin)     )
    outHist(      log(Jout)  )
    #^^^^^^^^^^^^^^^^^^^
    xlab('log(total intensity)')
    tit(simTag)
    legend();axisFontSize()
    
    subplot(1,2,2)
    #vvvvvvvvvvvvvvvvvvvvvv
    inHist(      log(JmaxIn)     )
    outHist(      log(JmaxOut)  )
    #^^^^^^^^^^^^^^^^^^^
    xlab('log(max intensity)')
    legend();axisFontSize()
    sf(os.path.join(figFolder,'intensityHist_%d_part_%d_frames'%(nPart,nFrames)))
    #######################################################
    
    #################### VELOCITY #####################
    figure(figsize = (20,20))
    subplot(2,2,1)
    #vvvvvvvvvvvvvvvvvvvvvv
    [n, b, patches] = inHist(      vProxyIn/Jin    )
    outHist(      vProxyOut/Jout   , bins = b)
    #^^^^^^^^^^^^^^^^^^^
    xlab('D$^2$T$^4$/I')
    tit(simTag)
    legend();axisFontSize()
    subplot(2,2,2)
    #vvvvvvvvvvvvvvvvvvvvvv
    [n, b, patches] = inHist(      vProxyIn/(Din*Jin)    )
    outHist(      vProxyOut/(Dout*Jout) , bins = b)
    #^^^^^^^^^^^^^^^^^^^
    xlab('DT$^4$/I$_{max}$')
    legend();axisFontSize()
    subplot(2,2,3)
    #vvvvvvvvvvvvvvvvvvvvvv
    [n, b, patches] = inHist(      log(vProxyIn/Jin)    )
    outHist(      log(vProxyOut/Jout )  , bins = b)
    #^^^^^^^^^^^^^^^^^^^
    xlab('log(D$^2$T$^4$/I)')
    legend();axisFontSize()
    subplot(2,2,4)
    #vvvvvvvvvvvvvvvvvvvvvv
    [n, b, patches] = inHist(      log(vProxyIn/(Din*Jin))    )
    outHist(      log(vProxyOut/(Dout*Jout)) , bins = b)
    #^^^^^^^^^^^^^^^^^^^
    xlab('log(DT$^4$/I$_{max}$)')
    legend();axisFontSize()
    sf(os.path.join(figFolder,'velocityProxyHist_%d_part_%d_frames'%(nPart,nFrames)))
    ######################################################################
    #::::::::::::::::::::::::::::::::::::::   END HISTOGRAMS ::::::::::::::::::::::::::::::::::
#    try:
#		sly = slice(0,None, 100)
#		mlm.scat(r_['0,2', Tin.ravel()[sly], vProxyIn.ravel()[sly]/Vin.ravel()[sly], Jin.ravel()[sly], Tin.ravel()[sly]], axeLab = ['T', 'BB Intensity', 'Intensity'], f = mlm.fig('Intensity Scatter')) 
#		mlm.scat(r_['0,2', Tin.ravel()[sly], Din.ravel()[sly], Jin.ravel()[sly], vProxyIn.ravel()/Jin.ravel()[sly]], axeLab = ['T', 'D', 'I', 'Vproxy'], f = mlm.fig('Input Scatter')) 
#		mlm.scat(r_['0,2', Tout[sly],Dout[sly], Jout[sly], vProxyOut[sly]/Jout[sly]], axeLab = ['T', 'D', 'I', 'Vproxy'], f = mlm.fig('Output Scatter')) 
#		mlm.ml.show()
#    except:
#		show()
