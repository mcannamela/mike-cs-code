#here we implement classes for simple raytracing lenses
import pdb

from helpers import *
from numpy import *
import scipy
scipy.pkgload()
from pylab import *

class lens(object):
    """
    bas class for all lens objects
    """
    def __init__(self, focalLength = 1, fNumber = 1.5625, imagePlane = 1):
        #lens parameters in meters
        self.fNumber =fNumber                            #aperture of the system
        self.focalLength = focalLength
        self.diameter = focalLength/fNumber                        #physical diameter of the aperture
        self.imagePlane = imagePlane                         #distance from lens to image plane
        self.objectPlane = (1/focalLength-1/imagePlane)**-1     #distance from lens to object plane, thin lens approximation
        self.magnification = self.imagePlane/self.objectPlane #magnification of the system
        
    def tfFactory(self, positionDict, wavelength, pixelSize):
        """
        compute the transfer function for every location in locationArray
        locationArray - array of 3-vectors giving
        """
        pass
        
            
     
class transferFunction(object):
    """
    make a tiny class to be the transfer functions
    """
    def __init__(self, kernel, dcGain = 1,wavelength = None ):
        """
        initialize with wavelengths and convolution kernel
        """
        self.wavelength = wavelength
        self.kernel = atleast_1d(float64(kernel)).copy()
        self.kernel/=sum(self.kernel)
        self.kernel*=dcGain
        
        
    def __call__(self, img):
        """
        convolve the img with the kernel and return the output image
        """
        pass
        
class bwTransferFunction(transferFunction):
    def __call__(self, img):
        idx = argmax(sum(img,axis = 0))
        try:
            outputImgPrototype = scipy.ndimage.convolve1d(r_['0,1',zeros(self.kernel.shape[0]-1), atleast_1d(squeeze(img[:,idx])),zeros(self.kernel.shape[0]-1)], self.kernel, mode = 'constant')
        except ValueError:
            pdb.set_trace()
            
        scales = squeeze(sum(img,axis = 0)/sum(img[:,idx]))
        return outputImgPrototype[:,newaxis]*scales
  
class rayTracingLens(lens):
    """
    base class for all ray transfer matrix lens objects, it is fully functional but uses 
    the simplest possible assumptions:
    
    1) shift invariant, rectangular point spread function (psf) 
    2) no diffraction in psf
    3) exit pupil as large as entrance pupil
    4)thin lens
    
    as input, a lens object should take a list of source points (y,z) and their intensities I,
    and a sample length delta
    where
    
    y = list of distances off the optical axis
    z = list of distances along the optical axis with z=0 the object plane
    I = list of relative intensity of light radiated from each (y,z)
    
    the output will be a 1d array y representing the image of the source points
    """
    
    def __init__(s,**kwargs ):
        #lens parameters in meters
        lens.__init__(s,**kwargs)
        #can precompute the last two transfer matrices
        lensMat = mat(array([1, 0, -1/s.focalLength, 1]).reshape(2,2)) #thin lens transfer matrix
        tubeMat =mat(s.freeSpaceTransferMatrix(s.imagePlane) )#free space propagation matrix: lens to image plane
        
        #successive transforms left multiply one another
        s.fixedTransferMatrix = array(tubeMat*lensMat)
    
    def freeSpaceTransferMatrix(self,d):
        """
        the transfer matrix for propagation through free space
        """
        d = atleast_1d(d)
            
        s = (1,)+d.shape+(1,)
        return r_[r_[str(-1), ones(s), d[newaxis,...,newaxis]], r_[str(-1),zeros(s), ones(s)]]
        
    def boundingRayAngles(s,y=0,z=0):
        """
        compute minimum and maximum ray angles for given off axis distance y (positive up)
        and out of focal plane distance z (positive away from lens)
        """
        thetaMax = abs(arctan(  (s.diameter/2-y)/(s.objectPlane-z)  ))
        thetaMin = -abs(arctan(  (s.diameter/2+y)/(s.objectPlane-z)  ))
        return r_[thetaMax[newaxis,...], thetaMin[newaxis,...]]
        
    def raytrace(self,y,z,theta):
        """
        trace the rays defined by (y,theta) where
        y,z,theta are all arrays of the same shape
        """
        n = y.ndim+2
        s = y.shape
        
        raysIn = r_[str(-1), y[newaxis,...,newaxis],theta[newaxis,...,newaxis]]
        raysAtLens = sum(self.freeSpaceTransferMatrix(self.objectPlane-z)*raysIn, axis = -1)
        TM = self.fixedTransferMatrix.copy()
        for i in range(y.ndim):
            TM = TM[:,newaxis,...]
        
        
        raysOut = sum(TM*transpose(raysAtLens[...,newaxis], (n-1,)+tuple(range(1,n-1))+(0,)), axis = -1)
        
        raysIn = transpose(raysIn, (raysIn.ndim-1,)+tuple(range(1,raysIn.ndim-1))+(0,))
        raysIn = reshape(raysIn, raysIn.shape[:-1])
        return (raysIn, raysAtLens, raysOut)
    
    def rayplot(self, y, z, theta = None):
        """
        plot the rays specified by y, z, and theta
        if theta is None, used the boundingRayAngles 
        """
        if theta==None:
            theta = self.boundingRayAngles(y.flatten(),z.flatten())
        
        rays = []
        if theta.ndim>y.ndim:
            for i in range(theta.shape[0]):
                rays+=[self.raytrace(y.flatten(),z.flatten(),theta[i])]
        else:
            rays+=[self.raytrace(y.flatten(),z.flatten(),theta)]
            
        on = ones_like(z.flatten())
        x = r_['0,2',z.flatten(), self.objectPlane*on, (self.objectPlane+self.imagePlane)*on]
        
        for r in rays:
            plot(1e3*x, 1e3*r_['0,2', r[0][0], r[1][0], r[2][0]])
        
        return rays
        

class rtmFSLens(rayTracingLens):
    """
    the flux sentinel's lens
    """
    def __init__(self, f = .075, fNr = 1.5625, ip = .188):
        """
        use flux sentinel values for defaults
        """
        rayTracingLens.__init__(self,focalLength = f, fNumber = fNr, imagePlane = ip)
        
    def tfFactory(self, positionDict, wavelength, pixelSize):
        """
        compute the transfer function for every location in locationArray
        locationArray - array of 3-vectors ()
        """
        #temporaries to hold position for convenient reading
        y = positionDict['verticalPosition']
        z = positionDict['lateralPosition']
        
        #viewFac is the fractional area of the lens in a sphere centered on the particle, approximate for off axis particles
        viewFac = .5*self.diameter/sqrt(self.diameter**2+(self.objectPlane-z)**2)
        #viewFac = ones_like(z)
        
		
        #angles for the rays
        theta = self.boundingRayAngles(y,z)
        
        #do the raytracing
        upperRays = self.raytrace(y,z,theta[0])
        lowerRays = self.raytrace(y,z,theta[1])
        
       	#compute  the kernels from the difference between upper and lower rays
        kernelWidth = abs(upperRays[-1][0]-lowerRays[-1][0])/pixelSize
        kernelMeans = kernelWidth/2
        kernelStdDev = kernelWidth/6
        kernelSupport = float64(arange(ceil(np.max(kernelWidth))))
        
        #allocate array for the transfer functions, loop to create them
        tfArray = empty(y.shape, dtype = 'object')
        for idx, tf in ndenumerate(tfArray):
            
            if kernelStdDev[idx]<=1./6.:
                tfArray[idx] = bwTransferFunction(atleast_1d(1), dcGain = viewFac[idx])
            else:
                tfArray[idx] = bwTransferFunction(scipy.stats.distributions.norm(loc = kernelMeans[idx], 
                                                                                                                                     scale = kernelStdDev[idx]).pdf(
                                                                                                                                         kernelSupport[:ceil(kernelWidth[idx])]  ), 
                                                                                                                                    dcGain = viewFac[idx])
                                                                                            
                                                                                
        return tfArray
        
        
if __name__== "__main__":
    close("all")
    L = rtmFSLens()
    z = linspace(0, 3e-3, 3)
    y = zeros_like(z)
    figure(1)
    L.rayplot(y,z)
    
    xlabel('distance from lens plane, mm')
    ylabel('distance off optical axis, mm')
    title('maximum and minimum rays from defocused locations')
    axis('equal')
    
    figure(2)
    r = L.raytrace(zeros(1000),1e-3+zeros(1000), linspace(0,L.boundingRayAngles(0,1e-3)[0], 1000))
    hist(r[-1][0]*1e3, 100)
    xlabel('final distance off optical axis, mm')
    ylabel('frequency')
    show()   