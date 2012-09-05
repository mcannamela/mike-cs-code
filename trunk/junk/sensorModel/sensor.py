#this model simulates the output of a sensor 
print 'importing utility modules'
import pdb
import os
import sys
import cPickle
import time 

print 'importing helpers'
import helpers as hp

print 'importing numpy and scipy'
from numpy import *
import scipy

print'importing supporting user modules'
import particle as part
import particleModelEnsemble as pmEns
import modelEnsemble as ens

print 'importing lens'
import Lens 


print 'importing pylab'
from pylab import *

print 'import complete'


class sensor(object):
    def __init__(self):
        """
        initialize the sensor
        """
        hp.virtualWarning()
        
    def sense(self, modEns):
        """
        compute sensor output for the ensemble passed, calling getTFArray and transfer
        """
        #get sensor input by observing an ensemble of models
        self.observe(modEns)
        #get an array of transfer functions, one for each model in the array
        self.getTFArray()
        #apply the transfer functions to the sensor inputs, generate sensor output
        self.transfer()
       
        return self.sensorOutput
        
    def observe(self, modEns):
        
        S = modEns.modelArray.shape
        
        self.sensorInput = dict(zip(['axialPosition','verticalPosition','lateralPosition', 'intensityProfile'],
                         [zeros(S),zeros(S),zeros(S), zeros_like(modEns.modelArray) ]))
        for idx, M in ndenumerate(modEns.modelArray):
            self.sensorInput['axialPosition'][idx] = M.PS.p.axialPosition()
            self.sensorInput['verticalPosition'][idx] = M.PS.p.verticalPosition()
            self.sensorInput['lateralPosition'][idx] = M.PS.p.lateralPosition()
            
            self.sensorInput['intensityProfile'][idx] = M.PS.p.intensityProfile((self.magnification*M.PS.p.diameter()/self.pixelSize),
                                                                                         self.wavelength*1.e-9)
                
            
       
        
    def getTFArray(self):
        """
        compute sensor transfer functions for each member of ensemble
        """
        hp.virtualWarning()
    
    def transfer(self):
        """
        apply the transfer function to all sensor inputs
        """
        hp.virtualWarning()
        
    
class opticalParticleSensor(sensor):
    def __init__(self, lens  = None, psnr = Inf, spectralSensitivityFile = 'spectralSensitivity.txt', gains = atleast_1d(1)):
        """
        setup lens, signal noise ratio, and spectral sensitivity
        """
        self.lens = lens
        self.psnr = psnr
        self.magnification = lens.magnification
        
        #pixel size(s) in meters
        self.pixelSize = (1e-6,)
        
        #the area of each pixel in m**2
        self.pixelArea = (1,)
        
        #spacing between pixel arrays, if any
        self.pixelSpacing = (0,)
        
        #shape of the pixel array
        self.sensorShape = (256,)
        
        #gain to apply to the transferred images
        self.gains = gains
        
        #need to know the spectral sensitivity of the pixels
        self.setSpectralSensitivity(spectralSensitivityFile)
        
    def setSpectralSensitivity(self, f):
        """
        read the spectral sensitivites from a file
        """
        try:
            ssf = open(f)
        except IOError:
            ssf = open(os.path.join(os.pardir, 'sensorModel', f))
            
        ss = fromfile(ssf,'double',-1,' ').reshape(5,-1)
        self.wavelength          = ss[0]
        self.spectralSensitivity = ss[1:]
        ssf.close()    
    
class basicFluxSentinel(opticalParticleSensor):
    def __init__(self, **kwargs):
        try:
            if kwargs['lens']==None:
                kwargs['lens'] == Lens.rtmFSLens()
        except KeyError:
            kwargs = dict(zip(kwargs.keys()+['lens'],kwargs.values()+[Lens.rtmFSLens()]))
        
        opticalParticleSensor.__init__(self, **kwargs)
        self.sensorShape = (4,8196)
        self.pixelSize = 5e-6
        self.gains = 3e13*ones(self.spectralSensitivity.shape[0]) 
        self.axialPos = .1
        self.sensorOutput = dict(zip(['particleImages', 'imageSlices', 'sensorImage' ],
                                                             [None,
                                                             None,
                                                             (zeros(self.sensorShape))]))
		
    def getTFArray(self):
        """
        lens will provide the transfer functions
        """
        self.tfArray = self.lens.tfFactory(self.sensorInput, self.wavelength, self.pixelSize)
    
    def transfer(self):
        """
        create the sensorOutput dictionary of color images for each particle, plus the accumulated image
        """
        #sensor outputs individual particle images, the acutal sensed image, and slices for each particle image into the sensor image
        self.sensorOutput = dict(zip(['particleImages', 'imageSlices', 'sensorImage' ],
                                                             [empty(self.sensorInput['intensityProfile'].shape, dtype='object' ),
                                                             empty(self.sensorInput['intensityProfile'].shape, dtype ='object'),
                                                             (zeros(self.sensorShape))]))
        #temporaries for readability    
        pImg = self.sensorOutput['particleImages']
        sImg = self.sensorOutput['sensorImage']
        sly = self.sensorOutput['imageSlices']
        
        #this distance divided by the velocity (already embedded in the intensity profiles) gives the time the particle light is incident on the sensor
        timeFac = self.pixelSize/self.lens.magnification
        
        #loopt over all particle input images
        for idx, img in ndenumerate(self.sensorInput['intensityProfile']):
           
           #apply the lens transfer function to each input image
           pImg[idx] = timeFac*self.gains[:,newaxis]*sum(self.tfArray[idx](img)[newaxis,...]*self.spectralSensitivity[:,newaxis,...], axis = 2)
           
           #locate the center of each image. center of the sensor array is 0, positive direction is up
           imgLoc = round(self.magnification*self.sensorInput['verticalPosition'][idx]/self.pixelSize+sImg.shape[1]/2)
           
           #cases for odd and even length particle images, build slices that don't go off the edge
           if mod(pImg[idx].shape[1],2)==0:
               sly[idx] = int32(arange( max([0,imgLoc - pImg[idx].shape[1]/2]), min([sImg.shape[1],imgLoc + pImg[idx].shape[1]/2])))
           else:
               sly[idx] = int32(arange( max([0,imgLoc - pImg[idx].shape[1]/2]), min([sImg.shape[1],imgLoc + 1+pImg[idx].shape[1]/2])))
           
           #add to the sensor image at the proper location, slicing out only the parts of the image that didn't go off the edges
           sImg[:,sly[idx]]+= pImg[idx][:,int32(sly[idx]-imgLoc+pImg[idx].shape[1]/2)]
           
          
        
    def setNoiseField(self):
        """
        use psnr to set up the noise field
        """
        pass

if __name__=='__main__':
    
    #gridded PE
    #PE = pmEns.emittingParticleEnsemble(pmEns.basicGriddedEmittingParticleGenerator)
    
    #random PE
    PE = pmEns.emittingParticleEnsemble(pmEns.marginalRandomEmittingParticleGenerator)
    
    start = time.time()
    PE.generate(N = 1e3)
    #PE.generate(dMinMax = [20e-6, 60e-6],tMinMax = [2700, 4000], zMinMax = [0,1e-3], n = (3,3,1,1),nPixels =8196 )
    elapsed = time.time()-start
    print 'generated %d models in %f seconds: %.2f models/s'%(PE.modelArray.size, elapsed,PE.modelArray.size/elapsed )
    
    
    start = time.time()
    FS = basicFluxSentinel()
    elapsed = time.time()-start
    #print 'initialized flux sentinel in %f seconds'%elapsed
    
    start = time.time()
    D = FS.sense(PE)
    elapsed = time.time()-start
    print 'produced the image of  %d models in %f seconds: %.2f models/s'%(PE.modelArray.size, elapsed,PE.modelArray.size/elapsed )
    
    close('all')
    #plot the sensor image 
    figure()
    c = ['k', 'r', 'g', 'b']
    [plot((FS.pixelSize*arange(FS.sensorShape[1])/FS.magnification)*1e6, D['sensorImage'][i], 
                                            color = c[i], linewidth = 2) for i in range(D['sensorImage'].shape[0])]

    ylim([0,255])
    xlabel('vertical coordinate, $\mu$m')
    ylabel('intensity')
    title('sensor image')
    
    
    #trends of intensity with diameter and temperature
    figure(figsize = (20,20))
    x = zeros(PE.modelArray.size)
    y = zeros_like(x)
    for i in range(PE.modelArray.size):
        x[i] = 1e6*PE.modelArray.ravel()[i].PS.p.diameter()
        y[i] = sqrt(sum(D['particleImages'].ravel()[i][0]))
    subplot(1,2,1)
    plot(x,y , 'o')
    xlabel('diameter, $\mu$m')
    ylabel('$\sqrt{total Intensity}$')
    title('trend of square root of intensity with diameter')
        
    
    for i in range(PE.modelArray.size):
        x[i] = PE.modelArray.ravel()[i].PS.p.surfaceTemperature
        y[i] = sum(D['particleImages'].ravel()[i][0])**.25
    subplot(1,2,2)
    plot(x,y , 'o')
    xlabel('surface temperature, K')
    ylabel('totalIntensity $^{.25}$')
    title('trend of fourth root of intensity with temperature')
	
    figure()
    for i in range(PE.modelArray.size):
	    M = PE.modelArray.ravel()[i].PS.p
	    x[i] = (M.diameter()**2)*(M.surfaceTemperature**4)/M.axialVelocity()
	    y[i] = sum(D['particleImages'].ravel()[i][0])
    plot(x,y , 'o')
    xlabel('D$^2$T$^4$/V')
    ylabel('totalIntensity')
    title('trend of intensity with D$^2$T$^4$/V')
	
    figure()
    for i in range(PE.modelArray.size):
	    M = PE.modelArray.ravel()[i].PS.p
	    x[i] = (M.diameter()**2)*(M.surfaceTemperature**4)/M.axialVelocity()
	    y[i] = np.max(D['particleImages'].ravel()[i][0])
    plot(x,y , 'o')
    xlabel('D$^2$T$^4$/V')
    ylabel('maxIntensity')
    title('trend of max intensity with D$^2$T$^4$/V')
    
    
    #transfer function kernel diagnostics
    figure(figsize = (20,20))
    y2 = zeros_like(y)
    for i in range(PE.modelArray.size):
        x[i] = PE.modelArray.ravel()[i].PS.p.lateralPosition()*1e6
        y[i] = sum(FS.tfArray.ravel()[i].kernel)
        y2[i] = FS.tfArray.ravel()[i].kernel.shape[0]*FS.pixelSize*1e6/FS.magnification
		
        
    subplot(1,2,1)
    plot(x,y,'o')
    xlabel('lateral position, $\mu$m')
    ylabel('kernel area')
    
    subplot(1,2,2)
    plot(x,y2,'o')
    xlabel('lateral position, $\mu$m')
    ylabel('kernel width, $\mu$m')
    title('size of lens transfer function kernel')
    show()
	
	
