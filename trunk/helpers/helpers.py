#define some helpers for python 
import numpy as np
from numpy import *
rfft = fft.rfft
irfft = fft.irfft
fftshift = fft.fftshift
import scipy
scipy.pkgload()
import pdb

def ngrid(X):
    """
    like meshgrid for n dimensions
    X is an n-tuple of e.g. (x,y,z) points to arrange in mesh
    """
    G = ones((len(X),)+tuple([x.shape[0] for x in X]))
    for i in range(len(X)):
        sly = []
        for j in range(len(X)):
            if j==i:
                sly+= [slice(None)]
            else:
                sly+=[newaxis]    
                
        G[i] *=  X[i][sly]
    return G

def grid2idx(g):
    idx = zeros_like(g[0])
    rg = array(g.shape)[1:]
    gPrime = array([g[i]-amin(g[i]) for i in range(g.shape[0])])
    
    
    rgIdx = flipud(arange(g.shape[0])).tolist()
    
    for i in range(g.shape[0]):
        rgIdx.pop()
        if rgIdx==[]:
            factor = 1
        else:
            factor = prod(rg[rgIdx])
        
        idx+= factor*gPrime[i]
        
    return idx
    

    
def x2gamma(x):
    X = x.T    
    rg = array([1+amax(X[i])-amin(X[i]) for i in range(X.shape[0])])
    rgIdx = flipud(arange(X.shape[0])).tolist()
    
    XPrime = array([(X[i]-amin(X[i]))/rg[i] for i in range(X.shape[0])])
    
    gamma = zeros_like(X[0])
    for i in range(X.shape[0]):
        rgIdx.pop()
        if rgIdx==[]:
            factor = 1
        else:
            factor = prod(rg[rgIdx])
        
        gamma+= factor*XPrime[i]
        
    return gamma
    
    
    
def vectorQuantize(x, downsample = 1, nClusters = 100):
    """
    use clustering to quantize x, treating vectors along the 0 axis as the observations
    """
    obs = zeros((prod(x.shape[1:]), x.shape[0]))
    y = zeros_like(x)
    c = zeros_like(x[0])
    d = zeros_like(x[0])
    cnt = 0
    for idx, a in ndenumerate(x[0]):
        obs[cnt] = x[(slice(None),)+idx]
        cnt+=1
    
    codebook , dist= scipy.cluster.vq.kmeans(obs[slice(0,None, downsample),:], nClusters)
    codes, dist = scipy.cluster.vq.vq(obs, codebook)
    
    quantized = codebook[codes,:]
    
    cnt = 0
    
    for idx, a in ndenumerate(x[0]):
        y[(slice(None),)+idx] = quantized[cnt]
        c[idx] = codes[cnt]
        d[idx] = dist[cnt]
        cnt+=1
    return (y, c,d,codebook)

def matchImageSize(master, slave):
    """
    use map_coordinates to resample slave along axis to be the same size 
    as master
    """
    sly = []
    for i in range(slave.ndim):
        sly+=[slice(0,slave.shape[i],master.shape[i]*1j)]
        
    G = mgrid[sly]
    g = []
    for i in range(G.shape[0]):
        g+=[G[i].flatten()]
        
    return scipy.ndimage.map_coordinates(slave, g, mode = 'nearest').reshape(master.shape)
    
    

        
class oneD_LUT(object):
    """
    lookup-table for a one-d array
    """
    def __init__(self,x,y):
        self.x = x
        self.y = y
    def __call__(self, xi):
        virtualWarning()
        return None
    
class indicialOneD_LUT(oneD_LUT):
    """
    zeroth order interpolation for indicial coordinates in a one d array 
    """
    def __init__(self,x,y):
        """
        interpolate input arrays on indicial grid
        """
        self.x = arange(round(np.max(x)))
        self.y = float32(scipy.interpolate.interp1d(x,y)(self.x))
        
    def __call__(self,xi):
        """
        lookup yi from the integer part of xi
        """
        return self.y[int16(xi)]
        


    
def pnorm(x, dim=0, p=2.0):
    """
    return the p-norm of an array along dimension dim
    """
    return  pow(sum(pow(abs(x),p) ,dim), 1/double(p))

def harmonicMean(x, w = None, dim = 0):
    """
    weighted harmonic mean along the specified axis
    """
    if w == None:
        w = ones_like(x)
        
    m = sum(w, axis = dim)/sum(w/x,axis = dim)
    return m



def excise(x,r,idx):
    """
    cut a cube of side length 2r out of an array x around a point idx
    """
    #slice out within bounds 
    
    b = [  [max([   idx[i]-ceil(r),     0                     ]),  
           min([    idx[i]+ceil(r),    x.shape[i]-1   ])]
                                for i in range(x.ndim) ]
    
    sly = [slice(b[i][0], b[i][1], 1) for i in range(x.ndim)]                           
    tup = tuple([int32(arange(b[i][0],b[i][1])) 
                                for i in range(x.ndim) ])
#    print b
#    print sly
#    print tup 
#    
    
    return (sly, x[sly], tup)

def virtualWarning():
    print "this is a virtual function! subclass and overwrite"
    
def uniqueRowsIdx(x):
    """
    return indices of the unique rows of the sorted 2d array x
    """
    d = diff(x,axis = 0)
    idx = flatnonzero(pnorm(d,1)!=0)
    return idx
    
def oneAxisConvolve(x, y, dim = 0):
    """
    convolve 2 arrays of the same size along dimension dim using the fft
    """
    z = fftshift(irfft(rfft(x,x.shape[dim],dim )*rfft(y,y.shape[dim],dim),y.shape[dim],dim),[dim])
    return z

def diffsForIntegral(x):
    """
    compute the widths for nonuniform rectangle approximation to an integral
    """
    xd = diff(concatenate((x[0]-diff(x[:2]), x, x[-1]+diff(x[-2:]))))
    return (xd[:-1]+xd[1:])/2
    

def decimate(x,factor):
    """
    downsample x by an integer factor, applying an antialiasing filter at 1/factor
    """
    #make the bands, allowing 9 percent of the cutoff for a transition band
    cutoff = .499/factor
    passBand = [0, .9*cutoff]
    stopBand = [1.05*cutoff, .499]
    
    #get the taps using remez exchange
    taps = scipy.signal.remez(32, array(passBand+stopBand),[1, 0])
    
    #filter once each way for zero phase
    y = scipy.signal.lfilter(flipud(taps), 1,scipy.signal.lfilter(taps,1,x))
    
    #downsample
    return y[slice(0,y.shape[0],factor)]
    

def lognormalParameters(mean, stdDev):
    #typecast to double
    variance = double(stdDev)**2
    mean = double(mean)
    #compute lognormal parameters
    mu = log(mean)-.5*log(1+variance/mean**2)
    sig = sqrt(log(variance/mean**2 +1))
    return (mu, sig)
    
def randomLognormal(N,mean = 1., stdDev = 1.):
    """
    generate N samples of a lognormal distribution with expected value mean and standard deviation stdDev
    """
    #typecast to double
    variance = double(stdDev)**2
    mean = double(mean)
    
    #compute lognormal parameters
    mu = log(mean)-.5*log(1+variance/mean**2)
    sig = sqrt(log(variance/mean**2 +1))
    
    #exponentiate the normally distributed samples to get a lognormal distribution
    return exp(random.standard_normal(N)*sig+mu)

def lognormalPdf(x,mean = 1., stdDev = 1.):
    """
    return the value of the lognormal pdf at x
    """
    
    #typecast to double
    variance = double(stdDev)**2
    mean = double(mean)
    
    #compute lognormal parameters
    mu = log(mean)-.5*log(1+variance/mean**2)
    sig = sqrt(log(variance/mean**2 +1))
    
    return exp(-.5*((log(x)-mu)/sig)**2)/(sig*x*sqrt(2*pi))

def normalPdf(x,m=0,s=1):
    return exp(-.5*( (x-m) /s)**2)/(s*sqrt(2*pi))
    
# def peaks(x, n = 1, dim = 0):
    # """
    # along dimension dim of array x, locate peaks in x that are greater than n neighbors in each direction, 
    # returning a vector of indices idx into dimension dim of x corresponding to said peaks
    # """
    # S = [slice(None,None) for i in range(x.ndim)]
    
    # b = 0*x;
    # b = [ x[i]>x[i-n]
    