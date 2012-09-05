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
    simTag =  'meshedPlume'
    
    figFolder = os.path.join(os.curdir, simTag)
    if not os.path.isdir(figFolder):
        os.mkdir(figFolder)
    close('all')
    
    #####################   NR MESH POINTS ######################
    n = (10,10,10,1) # (D, T, z, V)
    ############################################################
    
    ################### INITIALIZE OBJECTS ##################################
    PE = pmEns.emittingParticleEnsemble(pmEns.basicGriddedEmittingParticleGenerator)
    IS = invSens.adHocInverseSentinel_2dLUT(highPassOn = True, lowPassOn = True)
    
    print IS.FS.gains
    IS.FS.gains*= 0
    IS.FS.gains+= 3*(3e13)
    print IS.FS.gains
    
    IS.highPassConstant = .99
    #######################################################################
    
     ##################  INITIALIZE DATA ARRAYS ###########################
    # input
    Jin = zeros(n) #intensity
    JmaxIn = zeros(n)
    
    ###################################################################
    
    ####################    SIMULATE ####################################
    
    start = time.time()
  
    #generate the ith ensemble
    PE.generate(n=n)
         
   #do the inversion
    IS.particleWiseInvert(PE)
    
    #update input intensities
    for j, P in ndenumerate(IS.FS.sensorOutput['particleImages']):
        Jin[j] = sum(P[0])
        JmaxIn[j] = np.max(P[0])
        
    #compute velocity proxy numerator
    vProxyIn = PE.diameter()**2*PE.temperature()**2
    vProxyOut = IS.diameter()**2*IS.temperature()**4
    ######################################################################
    
    #////////////////////////////     SUMMARY AND FRAME PLOTTING /////////////////////////////////////
#    mlm.scat(r_['0,4,-1', squeeze(PE.temperature()), squeeze(1e6*PE.diameter()),squeeze(PE.lateralPosition())], axeLab = ['T', 'D', 'z'], f = mlm.fig('input mesh') )
#    mlm.scat(r_['0,4,-1', squeeze(IS.temperature()), squeeze(1e6*IS.diameter()),squeeze(PE.lateralPosition())], axeLab = ['T', 'D', 'z'], f = mlm.fig('output mesh') )
#    mlm.scat(r_['0,4,-1', squeeze(PE.temperature()), squeeze(1e6*PE.diameter()),squeeze(PE.lateralPosition()), 1e6*squeeze(IS.diameter()-PE.diameter())], axeLab = ['T', 'D', 'z'], f = mlm.fig('diameter error') )
    mlm.surf3(1e6*PE.diameter()[:,0,:,0], PE.lateralPosition()[:,0,:,0],1e6*(IS.diameter()[:,0,:,0]-PE.diameter()[:,0,:,0]), axeLab = ['D', 'z', 'error'], f = mlm.fig('error in diameter') )
    
    
    #::::::::::::::::::::::::::::::::::::::   HISTOGRAMS ::::::::::::::::::::::::::::::::::
    
    ############ TEMP, DIAMETER#######################
    figure(figsize = (20,20)); subplot(1,2,1)
    #temprature
    #vvvvvvvvvvvvvvvvvvvvvv
    inHist(      IS.temperature().ravel() -PE.temperature().ravel()    )
    #^^^^^^^^^^^^^^^^^^^
    xlab('error in temperature, K');    legend();        axisFontSize()
    tit(simTag)
    
    #diameter
    subplot(1,2,2)
    #vvvvvvvvvvvvvvvvvvvvvv
    inHist(      1e6*IS.diameter().ravel()-1e6*PE.diameter().ravel()   )
    #^^^^^^^^^^^^^^^^^^^
    xlab('error in diameter, $\mu$m');legend();axisFontSize()
    
    figure(figsize = (20,20)); subplot(2,2,1)
    #temprature
    #vvvvvvvvvvvvvvvvvvvvvv
    inHist(      PE.temperature().ravel()    )
    outHist(      IS.temperature().ravel()   )
    #^^^^^^^^^^^^^^^^^^^
    xlab('temperature, K');    legend();        axisFontSize()
    tit(simTag)
    
    #diameter
    subplot(2,2,2)
    #vvvvvvvvvvvvvvvvvvvvvv
    [nDin, b, patches] = inHist(      1e6*PE.diameter().ravel()   )
    [nDout, b, patches] = outHist(      1e6*IS.diameter().ravel() )
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
    [nDout, b, patches]  = inHist(      1e6*PE.diameter().ravel()   )
    outHist(       1e6*IS.diameter().ravel() , bins = b)
    #^^^^^^^^^^^^^^^^^^^
    xlab('diameter, $\mu$m');legend();axisFontSize()
    tit('forced binning')
    sf(os.path.join(figFolder,'T-DHist'+str(n)))
    #######################################################
    
    #################### LOG INTENSITY #####################
    figure(figsize = (20,20))
    subplot(1,2,1)
    #vvvvvvvvvvvvvvvvvvvvvv
    inHist(      log(Jin.ravel())     )
    outHist(      log(IS.intensity.ravel())  )
    #^^^^^^^^^^^^^^^^^^^
    xlab('log(total intensity)')
    tit(simTag)
    legend();axisFontSize()
    
    subplot(1,2,2)
    #vvvvvvvvvvvvvvvvvvvvvv
    inHist(      log(JmaxIn.ravel())     )
    outHist(      log(IS.maxIntensity.ravel())  )
    #^^^^^^^^^^^^^^^^^^^
    xlab('log(max intensity)')
    legend();axisFontSize()
    sf(os.path.join(figFolder,'intensityHist'+str(n)))
    #######################################################
    
    #################### VELOCITY #####################
    figure(figsize = (20,20))
    subplot(2,2,1)
    #vvvvvvvvvvvvvvvvvvvvvv
    [h, b, patches] = inHist(      vProxyIn.ravel()/Jin.ravel()    )
    outHist(      vProxyOut.ravel()/IS.intensity.ravel()   , bins = b)
    #^^^^^^^^^^^^^^^^^^^
    xlab('D$^2$T$^4$/I')
    tit(simTag)
    legend();axisFontSize()
    subplot(2,2,2)
    #vvvvvvvvvvvvvvvvvvvvvv
    [h, b, patches] = inHist(      vProxyIn.ravel()/(PE.diameter().ravel()*Jin.ravel())    )
    outHist(      vProxyOut.ravel()/(IS.diameter().ravel()*IS.maxIntensity.ravel()) , bins = b)
    #^^^^^^^^^^^^^^^^^^^
    xlab('DT$^4$/I$_{max}$')
    legend();axisFontSize()
    subplot(2,2,3)
    #vvvvvvvvvvvvvvvvvvvvvv
    [h, b, patches] = inHist(      log(vProxyIn.ravel()/Jin.ravel())    )
    outHist(      log(vProxyOut.ravel()/IS.intensity.ravel() )  , bins = b)
    #^^^^^^^^^^^^^^^^^^^
    xlab('log(D$^2$T$^4$/I)')
    legend();axisFontSize()
    subplot(2,2,4)
    #vvvvvvvvvvvvvvvvvvvvvv
    [h, b, patches] = inHist(      log(vProxyIn.ravel()/(PE.diameter().ravel()*Jin.ravel()))    )
    outHist(      log(vProxyOut.ravel()/(IS.diameter()*IS.intensity.ravel())) , bins = b)
    #^^^^^^^^^^^^^^^^^^^
    xlab('log(DT$^4$/I$_{max}$)')
    legend();axisFontSize()
    sf(os.path.join(figFolder,'velocityProxyHist'+str(n)))
    ######################################################################
    #////////////////////////////////////////////////////////////////////////////////////////////////////////////
    mlm.ml.show()
   
    
