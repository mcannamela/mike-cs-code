from frameProcessing import *
from numpy import *
from scipy.interpolate import interp1d
from scipy.integrate import cumtrapz

from pylab import *
#import mlabMacros as mlm
import time


class particleImageMachine(object):
    """
    generates sharp particle images based on diameter
    """
    def __init__(self):
        #array of pixel x-locations        
        self.x = float64(arange(100))
        
        #pixel locations in normalized units on a fine mesh
        self.xHat   = linspace(-1.1,1.1,2200)
        
        #the shape function computed on xHat
        self.f_xHat = (1-self.xHat**2)
        
        #force shape function to be positive
        self.f_xHat[self.f_xHat<0] = 0
        
        #now we can take the square root
        self.f_xHat = self.f_xHat**.5
        
        #compute the image by integration of the shape function
        self.F_xHat = cumtrapz(self.f_xHat, self.xHat) 
        
        #set up an interpolator so we can make any sized particle we want
        self.F = interp1d(self.xHat[1:], self.F_xHat)
        
        #preallocate for the image
        self.J = zeros_like(self.x)
        
    def __call__(self, D):
        """
        D - diameter in pixels
        
        returns:
            the sharp image of D
        """
        
        #catch cases where diameter is less than one or two pixels, as these 
        #images are trivial and require padding with zeros
        if D<=1:
            return float64([0,1,0])
        elif D<=2:
            return float64([0,1,1,0])
        
        #get meshpoints for this particle
        xHat = self.get_xHat(D*.5)
        
        #any xHat outside the radius will have the value zero. this 
        #prevents bounds errors on the interpolator
        xHat[xHat<-1.01] = -1.01
        xHat[xHat>1.01] = 1.01
        
        #image is formed by differentiating F evaluated at xHat and padding with 
        #zeros on either end
        return float32(r_['0,1',zeros(1),diff(self.F(xHat)),zeros(1)])
        
    def get_xHat(self, R):
        """
        compute gridpoints (in pixels) for a particle of diameter R, 
        ensuring that we include an extra pixel for the endpoints
        """
        return linspace(-ceil(R)/R, ceil(R)/R, ceil(2*R))
    
    
class particleGrid(object):
    """
    represents a mesh over the particle state space (Diameter, distanceFromFocalPlane)
    z is the coordinate in the lens axis direction:
        z = 0, focal plane
        z<0 closer to camera
        z>0 further from camera
    """
    def __init__(self, D_tup = (1, 100, 50), Z_tup = (-1, 1, 50)):
        """
        D_tup - linspace style specifier for the grid in diameter, (dmin, dmax, nSteps)
        Z_tup - linspace style specifier for the grid in z-location (mm), (zmin, zmax, nSteps)
        """
        self.D_tup = D_tup
        self.Z_tup = Z_tup
        
        #preallocate for the sharp image        
        self.sharpIm = zeros(D_tup[-1], dtype = 'object')
        
        #create arrays of z and D values
        self.z_ = linspace(Z_tup[0], Z_tup[1], Z_tup[2])
        self.D_ = linspace(D_tup[0], D_tup[1], D_tup[2])
        
        #initialize the sharp image generator
        self.PIM = particleImageMachine()
        
        #create sharp images for every diameter
        for i,d in enumerate(self.D_):
            self.sharpIm[i] = float32(self.PIM(d))
        
        #grid up the (D,z) space
        self.G = mgrid[slice(D_tup[0], D_tup[1], 1j*D_tup[2]), 
                       slice(Z_tup[0], Z_tup[1], 1j*Z_tup[2])]
        
        #this class will provide an iterator over the particles in the grid ,
        #so we must initialize it              
        self.resetIter()
    
    def resetIter(self):
        #we'll use the ndenumerate iterator to get our indices when iterating
        self.it = ndenumerate(self.G[1])
            
    def __iter__(self):
        """
        standard iterator function
        """
        return self
        
    def next(self):
        """
        iterator, when called serves up the index (i, j), axial position (z),
        and sharp image of that particle (sharpIm) all together in a tuple
        
        """
        try:        
            (i,j),z = self.it.next()
            return (i,j,z,self.sharpIm[i])
            
        except StopIteration:
            self.resetIter()
            raise StopIteration
            
    def D(self):
        """
        getter for the two-d grid of diamter
        """
        return self.G[0]
    def Z(self):
        """
        getter for the two-d grid of diamter
        """
        return self.G[1]
            
    def zeroImArr(self):
        """
        return an array the same size as the particle grid having dtype object
        """
        return zeros(self.G[0].shape, dtype = 'object')
        
    def zeroArr(self):
        """
        return a zero array the same size as the particle grid
        """
        return zeros_like(self.G[0])
        
                    

class blurrifier(object):
    """
    make blurry the image of a particle at a specified depth and diameter
    """
    def __init__(self, K, z):
        """
        initialize with a list of 
            3 x nKern[i] kernels K 
        and an array of 
            z-locations z
        """
        self.K = K
        self.z = z
        
    def __call__(self, im, z):
        """
        choose the appropriate kernel from K based on z, and blur the image
        im
        
        returns:
            blim - 3 x nImage array holding the blurred image
            M - len 3 array holding the observed diameters for r,g,b
            K - len nKernel  array holding the (possibly interpolated) used for
                this z 
        """
        
        #interpolate between kernels for this particular z location
        k = self.blendKernels(z)
        
        #allocate for the blurry image and convolve with the kernel to obtain it
        blim = zeros((k.shape[0], k.shape[1]+im.shape[0]-1))
        for i in range(k.shape[0]):
            blim[i] = convolve(im, k[i])
        
        #obtain the observed diameters from the blurry image
        M = float32(self.measure(blim))
        
        #return the kernel as well
        K = float32(k)
    
        return (float32(blim), M, K)
                
    def closestKernel(self, z):
        """
        find the index into the list of kernels K that is closest to z
        """
        return argmin(abs(self.z-z))
    
    def bracketingKernels(self, z):
        """
        find the kernels for z locations bracketing the given value of z
        """
        
        #index below the given value of z
        LIdx = argmax(flatnonzero(self.z<=z))
        
        #index above the given value of z
        UIdx = LIdx+1
        
        #make lists of bracketing z values and kernels
        Zb = [self.z[LIdx], self.z[UIdx]]
        Kb = [self.K[LIdx], self.K[UIdx]]
        
        
        return (Zb, Kb)
        
        
    def blendKernels(self, z):
        """
        for a given value of z, find the kernels at bracketing z locations 
        and interpolate them to get a more precise kernel at z
        """
        try:
            #obtain the neighboring kernels         
            Zb, Kb = self.bracketingKernels(z)
            
            #make sure the narrowest kernel is first  in the lists K, Z        
            if Kb[0].shape[1]<=Kb[1].shape[1]:
                K = Kb
                Z = Zb
            else:
                K = [Kb[1], Kb[0]]
                Z = [Zb[1], Zb[0]]
            
            #lengths of the kernels
            N = (K[0].shape[1], K[1].shape[1])
            
            #peak locations of the kernels        
            p = [argmax(K[0][0]), argmax(K[1][0])]
            
            
            #we must pad each kernel with the appropriate number of zeros 
            #such that the peaks will line up when we do the interpolation between them
            
            #pad0 will be a list of left and right padding number of zeros for the 
            #smaller kernel
            pad0 = [max([p[1]-p[0], 0]), 
                   max([(N[1]-p[1]) - (N[0]-p[0]), 0])]
            
            #as pad0 for the larger kernel
            pad1 = [max([p[0]-p[1], 0]), 
                   max([(N[0]-p[0]) - (N[1]-p[1]), 0])]
            
            #padded kernels will have the peaks at the same pixel        
            padded_k0 = r_['1,2', zeros((3, pad0[0])), K[0], zeros((3, pad0[1]))]
            padded_k1 = r_['1,2', zeros((3, pad1[0])), K[1], zeros((3, pad1[1]))]
            
            #f is the fraction of the second kernel to use
            f = (Z[1]-z)/(Z[1]-Z[0])
            assert (f<=1 and f>=0), "kernel blending failed because interpolation factor is outside the interval (0,1) \n check the z values!"
            
            #interpolate linearly between kernels
            blendedK = (1-f)*padded_k0+f*padded_k1
    
            return  blendedK
        except IndexError:
            idx = self.closestKernel(z)
            return self.K[idx]
        

        
    def measure(self, im, frac = .5):
        
        #indices of maxes for each channel, and maxes for each channel        
        mxIdx = argmax(im, axis = 1)        
        mx = amax(im, axis = 1)[...,newaxis]
        
        #left and right bounding indicies for frac 
        L = int32(zeros((im.shape[0], 2)))
        R = int32(zeros((im.shape[0], 2)))
        
        #for each color, find the largest index on the left half 
        #that is less than the fraction of the max, 
        #and the smallest index on the right half that is greater than the fraction 
        #of the max; add one to get the bracketing indices
        for i in range(im.shape[0]):
            splitIdx = mxIdx[i]
            leftHalf = im[:,:splitIdx]
            rightHalf = im[:,splitIdx:]
        
            L[i, 0] = int32(amax(flatnonzero(leftHalf[i]<=mx[i]*frac)))
            L[i, 1] = L[i,0]+1
            R[i, 0] = int32(amax(flatnonzero(rightHalf[i]>mx[i]*frac)))+splitIdx
            R[i, 1] = R[i,0]+1
        
        #needed to form the slicing tuples
        a3 = arange(3)
        
        #L0, R0, L1, R1 will slice out from each color the bracketing indices
        #for the left and right points bracketing frac
        L0 = (a3, L[:,0])
        L1 = (a3, L[:,1])
        
        R0 = (a3, R[:,0])
        R1 = (a3, R[:,1])
        
        #interpolate linearly to get the exact coordinate where the image is 
        #equal to frac
        mL = (im[L1]-im[L0])
        xL = L[:,0]+(frac*squeeze(mx)-im[L0])/mL
        
        mR = (im[R1]-im[R0])
        xR = R[:,0]+(frac*squeeze(mx)-im[R0])/mR
                
        return xR-xL
        
        
if __name__ == "__main__":
    print "reading z-scans at one array location"
    scanFolder = os.path.join(os.path.curdir, "Initial Scans", "00")
    SR = scanFolderReader()
    x,Z,Y = SR(scanFolder)
    
    KE = kernelExtractor()
    
    print "extracting kernels"
    KK = [KE(y)[-1] for y in Y]    
    
    def kernelVariance(K):
        v = zeros((3, len(K)))
        for j, k in enumerate(K):
            for i in range(K[0].shape[0]):
                y = k[i]/trapz(k[i]) 
                x = linspace(-k[i].shape[0]/2., k[i].shape[0]/2., k[i].shape[0])
                mu = trapz(x*y)
                v[i,j] = 6*sqrt(trapz(((x-mu)**2)*y))
        return v
        
    v = kernelVariance(KK)
    z0Idx = argmin(v[0])
    Z -= Z[z0Idx]
    
    PG = particleGrid(D_tup = (1, 100, 200), Z_tup = (-1.4, 3.5, 300))
    BL = blurrifier(KK, Z)
    
    IM = PG.zeroImArr()
    blendedKernels = PG.zeroImArr()
    D_obs = zeros(PG.zeroArr().shape+(3,))
    
    print "blurrifying"
    start = time.time()
    cnt = 0
    for i,j,z,im in PG:
        
#        if mod(cnt, 3) == 0:
#            plot(linspace(-im.shape[0]*.5, im.shape[0]*.5, im.shape[0]), im/amax(im))    
        IM[i,j], D_obs[i,j],blendedKernels[i,j] = BL(im, z)
#        if (z<.1) and (z>-.1) and (PG.D_[i]>70):
#            pdb.set_trace()
        cnt+=1
    elapsed = time.time()-start
    print "done. time was %d s for %d images"%(elapsed, cnt)
    print "that's %.2e s/image!"%(elapsed/double(cnt))
    

#    mlm.surf3(PG.G[0], PG.G[1], D_obs[:,:,0], f = mlm.fig('red diameter'))
#    mlm.ml.show()
        
    
        
        
        
               