from numpy import *
from helpers import virtualWarning
import scipy
scipy.pkgload('ndimage')
import pdb
from scipy.interpolate import interp1d

def regularize2d(x,y,z,n, nCut = (-1,-1)):
    """
    where x is has z.shape[0], y has z.shape[1] and form an an irregular mesh of z, 
    resample z on a regular mesh having dimensions given in n
    """
    assert z.ndim==2, 'z must have only 2 dimensions!' 
    xReg = linspace(x[0], x[nCut[0]-1], n[0])
    yReg = linspace(y[0],y[nCut[1]-1], n[1])
    
    zRegx = zeros((z[:nCut[0]].shape[0], n[1]))
    zReg = zeros(n)
    
    pol = lambda x_, z_, xi:  interp1d(x_,z_, kind = 'cubic')(xi)
    
    for i in range(zRegx.shape[0]):
        zRegx[i] = pol(y, z[i], yReg)

    

    for j in range(zRegx.shape[1]):
        zReg[:,j] = pol(x[:nCut[0]],zRegx[:,j], xReg)


    return (xReg, yReg, zReg)
    
        
    


    
class abstractInterpolator(object):
    def __init__(self, axesMin, axesMax, z, precision = 'single'):
        """
        set up the scaling for the axes 
        axesMin = sequence of lower bounds on the axes
        axesMax = sequence of upper bounds on the axes
        z - n-d data on a uniform grid
        """
        
        if axesMin == None:
            axesMin = zeros(z.ndim)
        if axesMax == None:
            axesMax = z.shape
            
        if precision == 'single':
            self.cast = self.toSingle
        elif precision =='double':
            self.cast = self.toDouble
        elif precision =='noCast':
            self.cast = self.noCast
        else:
            raise TypeError, 'unknown precision specified, pick "single" or "double"'
            
        self.mn = self.cast(axesMin)
        self.mx = self.cast(axesMax)
        self.z = self.cast(z)
        self.delta = self.cast((self.mx-self.mn)/(array(self.z.shape)-1))
        
        
    def toSingle(self,c):
        return float32(c)
        
    def toDouble(self,c):
        return float64(c)
        
    def noCast(self,c):
        return c
        
    def scale(self, X):
        """
        rescale to fractional indices
            X - ndarray with shape (z.ndim, ...)
            
        """
        if self.mn.ndim<X.ndim:
            
            self.mn = self.cast(resize(self.mn,
                        [self.z.ndim]+ones(X.ndim-1).tolist()))
            self.mx = self.cast(resize(self.mx,
                        [self.z.ndim]+ones(X.ndim-1).tolist()))
            self.delta = self.cast(resize(self.delta,
                        [self.z.ndim]+ones(X.ndim-1).tolist()))
        elif self.mn.ndim>X.ndim:
            self.mn.resize([self.z.ndim]+ones(X.ndim-1).tolist())
            self.mx.resize([self.z.ndim]+ones(X.ndim-1).tolist())
            self.delta.resize([self.z.ndim]+ones(X.ndim-1).tolist())
        
        self.Y = (X[:self.z.ndim]-self.mn)/self.delta
            

    def interp(self):
        """
        perform the interpolation with the current coords
        """
        virtualWarning()
        
    def __call__(self, X):
        """
        use the list of coordinate arrays to perform interpolation
        X - ndarray with shape (z.ndim, ...)
        """
        self.scale(X)
        return self.interp()

class cubicInterpolator(abstractInterpolator):
    def __init__(self,*args, **kwargs):
        """
        set up the scaling for the axes 
        axesMin = sequence of lower bounds on the axes
        axesMax = sequence of upper bounds on the axes
        z - n-d data on a uniform grid
        """
        abstractInterpolator.__init__(self,*args, **kwargs)
        self.coeffs = self.cast(scipy.ndimage.spline_filter(self.z))
        
        
    def interp(self):
        """
        use map_coordinates to apply cubic spline interpolation
        """
        return squeeze(scipy.ndimage.map_coordinates(self.coeffs, 
                                            self.Y[...,newaxis], 
                                            mode = 'nearest', 
                                            prefilter = False))
    
                                            
class zerothOrderInterpolator(abstractInterpolator):
     
    def interp(self):
        """
        interpolate by indexing...fast!
        """
        self.Y = int32(self.Y)
        
        return self.z[tuple(int32(self.Y))]
    
        
        
class linearInterpolator(abstractInterpolator):
    def __init__(self,*args, **kwargs):
        """
        compute the step size in each axis
        """
        abstractInterpolator.__init__(self,*args, **kwargs)
        
        tern = {True:0, False:1}
        
        self.D = array([r_[str(k)+','+str(self.z.ndim), 
                            diff(self.z, axis = k),
                            zeros(tuple(array(self.z.shape)-
                                (self.z.shape[k]-1)*int32(arange(self.z.ndim)==k))) ] 
                                for k in range(self.z.ndim)])
                            
        self.DPrime = array([diff(self.z, axis = k)[tuple([slice(tern[i==k],None) 
                            for i in range(self.z.ndim)])]
                            for k in range(self.z.ndim)])
        
        self.zPrime = self.z[tuple([slice(1,None) for k in range(self.z.ndim)])]-sum(self.DPrime,0)
    
    def scale(self, *args):
        abstractInterpolator.scale(self,*args)
        self.idx = int32(self.Y)
        self.Y-=self.idx
        self.idxTup = tuple(self.idx)
        
    def interp(self):
        zInt = zeros_like(self.Y)
        if all(self.Y)<.5:
            print 'forward case'
            return (self.z[self.idxTup]+
                sum(self.D[(slice(0,None),)+self.idxTup]*self.Y,axis = 0))
        else:
            print 'reverse case'
            return (self.zPrime[self.idxTup]+
                sum(self.DPrime[(slice(0,None),)+self.idxTup]*self.Y,axis = 0))
        
class linear1DGridInterpolator(object):
    def __init__(self, x, y, boundaryErrorOn = False):
        self.xMin = x[0]
        self.xMax = x[-1]
        self.y = y
        self.dy = diff(y)
        self.n = len(y)
        self.boundaryErrorOn = boundaryErrorOn
        
    def __call__(self, xi):
        I,f = self.idx(xi)
        try:
            yi = self.y[I]+f*self.dy[I]
        except:
            pdb.set_trace()
        
        return yi
    
    def idx(self, xi):
        x = float64((self.n-1)*(xi-self.xMin)/(self.xMax-self.xMin))
        I = int32(x)
 
        if self.boundaryErrorOn:
            if any(I<xMin) or any(I>xMax):
                raise ValueError, "interpolation point outside the domain!"
                
        f = x-I
        f[I>=(self.n-1)] = 1-1e-10
      
        I[I<0] = 0
        I[I>=(self.n-1)] = self.n-2
        
        return (I,f)
        
        
if __name__=='__main__':
    a = arange(12.).reshape(3,4)
    ax = r_['0,2', zeros(2), ones(2)]
    I = linearInterpolator(ax[0], ax[1], a)
    x = r_['0,2', array([0,.6,1]),array([0,.6,1])]
    print I(x)
    
    
    
    
        