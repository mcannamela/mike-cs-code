############################################################
############################################################
#here we have all the classes needed to extract a kernel from a frame
############################################################

from lensFolderIO import *
from numpy import *
from scipy import integrate
import pdb

class framePreProcessor(object):
    """
    perform basic frame operations: smoothing and rough extraction of ROI 
    """
    def __init__(self):
        """
        set up the smoothing kernel, 
        ROI halfwidth (roughChopHW),
        and the fraction of maximum frame value we should identify as the step
        """
        self.roughChopHW = 100        
        self.buildSmoothingKernel()
        self.stepFrac = .5

    def __call__(self, frame):
        """
        perform the sequence of operations and return the finished ROI
        """
        self.y = float32(frame)
        
        #divide by maximum frame value
        self.normalize()
        
        #make a smooth version before extracting the ROI
        self.smooth()
        
        #find the step region
        self.locateStep()
        
        #extract the step region
        self.roughChop()
        
        #make the smooth version of the step region
        self.smooth()
        
        return self.ySmooth.copy()

    def normalize(self):
        self.y/=amax(self.y, axis = 1)[...,newaxis]
        
    def locateStep(self):
        """
        find the first two pixels that bracket the 
        stepFrac value
        """
        LT = self.ySmooth[0]<=self.stepFrac
        GT = self.ySmooth[0]>=self.stepFrac
        b = LT[:-1]*GT[1:]
        self.stepPixel = amin(flatnonzero(b))
        
    def roughChop(self):
        """
        extract a region of radius roughChopHW around the step pixel
        """
        self.b = (max([0,self.stepPixel-self.roughChopHW]),
             min([self.y.shape[1]-1, self.stepPixel+self.roughChopHW]))
        self.y = self.y[:, self.b[0]:self.b[1]]   
        self.y-= self.y[:,0][...,newaxis]
        
        
    def smooth(self):
        """
        make a smooth version of the frame
        """
        y = []
        for i in range(self.y.shape[0]):
            y += [convolve(self.y[i], self.kern, mode ='valid')]
        self.ySmooth = zeros((len(y), y[0].shape[0]))
        for i in range(len(y)):
            self.ySmooth[i] = y[i]-y[i][0]
            
        self.ySmooth/=amax(self.ySmooth, axis = 1)[...,newaxis]

    def buildSmoothingKernel(self):
        """
        make the smoothing kernel that corresponds to the one used 
        in the flux sentinel algo. if that changes, this function must 
        be changed too!
        """
        k = float32(ones(3))
        self.kern = k.copy()
        for i in range(4):
            self.kern = float32(convolve(self.kern, k))
            
class kernelExtractor(object):
    """
    this class will extract the kernel from a raw frame. it uses the frame
    preprocessor class to extract the smooth ROI, then differentiates to find
    the kernel. 
    """
    def __init__(self):
        """
        set threshold for further ROI extraction, initialize the preprocessor
        """
        #only keep the top 1-fineChopThreshold fraction of the frame
        self.fineChopThreshold = .05

        self.PP = framePreProcessor()
        
    def __call__(self, frame):
        """
        extract a kernel from frame
        frame - 3 x N array of rgb values
        
        return:
                tuple of 
                    (locationOfKernelCenter, kernelWidth, kernel)
        """

        #pre-process the frame            
        self.y = self.PP(frame)

        #take derivative of the smooth frame
        self.differentiate()
        
        #make the peak of the differentiated frame 1           
        self.normalize()
        
        #pin down the kernel center by its peak
        self.locateStep()
        
        #reduce the ROI, retaining the kernel down to fineChopThreshold
        #of it's height
        try:
            self.fineChop()
        except ValueError:
            print "Error extracting central portion of kernel. Try setting a breakpoint here and viewing the frame and it's derivative"
            raise
            
        fwhm = self.fullWidthHalfMax()

        retval = (self.stepPixel, 
                  self.fullWidthHalfMax(),
                  float32(self.yPrime.copy()))
        return retval
        
    def locateStep(self):
        #step pixel is the start of the ROI plus the distance in the ROI to the 
        #kernel peak
        self.stepPixel = self.PP.b[0]+ round(mean(argmax(self.yPrime, axis = 1)))
        
    def normalize(self):
        self.yPrime/=amax(self.yPrime, axis = 1)[...,newaxis]
        
    def differentiate(self):
        #differentiate and force positive
        self.yPrime = diff(self.y, axis = 1)
        self.yPrime[self.yPrime<0] = 0
        
    def fineChop(self):
        #left and right indices for r, g, and b
        L,R = self.argFractionPeakHeight(self.fineChopThreshold)
        
        #slice from the leftmost to rightmost index 
        sly = slice(amin(L), amax(R)+1)
        
        #refine the ROI
        self.y = self.y[:, sly]
        self.yPrime = self.yPrime[:, sly]
        
        
    def argFractionPeakHeight(self, frac = .5):
        """
        moving outward from the kernel peak, find the pixels where the value of 
        the kernel falls below frac*100 percent of the peak height
        
        returns:
            L,R: length 3 arrays, left and right indices respectively for rgb
        """
        #indices of peaks in rgb, array of length 3
        mid = floor(self.yPrime.shape[1]/2.0)
        offset = min([20, mid-1, self.yPrime.shape[1]-mid-1])
        midSly = slice(mid-offset, mid+offset)
        mxIdx = argmax(self.yPrime[slice(None, None), midSly], axis = 1)+mid-offset     
        
        # length 3 array, value of the peak maximum
        mx = array([self.yPrime[i,mxIdx[i]] for i in range(3)])[...,newaxis]
        
        #pre-allocate for left and right indices
        L = int32(zeros(self.yPrime.shape[0]))
        R = int32(zeros(self.yPrime.shape[0]))
        
        #split the array into left and right halves, search each half separately
        #for the first pixel that falls below frac fraction of the max
        leftHalf = self.yPrime[:,:mxIdx[0]]
        rightHalf = self.yPrime[:,mxIdx[0]:]
        for i in range(self.yPrime.shape[0]):
            L[i] = int32(amax(flatnonzero(leftHalf[i]<=mx[i]*frac)))
            R[i] = int32(amin(flatnonzero(rightHalf[i]<=mx[i]*frac)))+mxIdx[0]
            
        return (L, R)
       
        
    def fullWidthHalfMax(self):
        """
        use argFractionPeakHeigth(.5) to find the width of the kernel at 
        half of its maximum value
        """
        L,R = self.argFractionPeakHeight()
        return R-L
        
    
        
if __name__ == '__main__':
    print "reading z-scans at one array location"
    scanFolder = os.path.join(os.path.curdir, "Initial Scans", "00")
    SR = scanFolderReader()
    x,z,Y = SR(scanFolder)
    
    
    KE = kernelExtractor()  
    
    print "extracting kernels"
    KK = [KE(y) for y in Y]
    
    for k in KK:
        plot(k[-1].T)
        
    
    
#    w = zeros((3, len(KK)))    
#    v = zeros((3, len(KK)))    
#    for j, k in enumerate(KK):
#        for i in range(3):
#            w[i, j] = k[1][i]
#            
#            y = k[-1][i]/trapz(k[-1][i])
#            
#            x = linspace(-k[-1][i].shape[0]/2., k[-1][i].shape[0]/2., k[-1][i].shape[0])
#            mu = trapz(x*y)
#            v[i,j] = 6*sqrt(trapz(((x-mu)**2)*y))
#    
#    figure()        
#    plot(w[0], 'r')
#    plot(w[1], 'g')
#    plot(w[2], 'b')
#    
#    figure()        
#    plot(v[0], 'r')
#    plot(v[1], 'g')
#    plot(v[2], 'b')
#    show()
        
    