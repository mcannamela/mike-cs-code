from numpy import *
from scipy.ndimage import gaussian_filter
from scipy.ndimage.filters import sobel
from plotMacros import * 
import matplotlib as mpl
#import mlabMacros as mlm
import pdb
from scipy.stats import norm as normRand 
import cPickle as pkl
import time
import networkx
from helpers import matchImageSize as match
from matplotlib.collections import EllipseCollection
import gc

def smoothg(x, t, ord = 0, output = None):
    return gaussian_filter(x, sigma = t**.5, order = ord, mode = 'nearest', output = output)
    
def gauss(x,mu,sig):
    return exp(-(x-mu)**2/(2*sig**2))/sqrt(2*pi*sig**2)

def multiGauss32(x,mu,sig):
    """
    x - nDim x nX, float32
    mu - nDim x nGauss
    sig - nGauss
    
    returns: sum of gaussians with centers mu and standard devs sig
    """
    mu32 = float32(mu)[:, newaxis, ...]
    sig32 = float32(sig)
    xarg = sum((x[...,newaxis]-mu32)**2, axis = 0)
    
    supp = unique(flatnonzero(any(xarg**.5<sig32, axis = -1)))
 
    G = sum(exp(-.5*xarg/sig32**2)/sqrt(2*pi*sig32**2), axis = -1)
    G[supp]/=sum(G[supp])

    return (supp, G[supp])
    
def symeig2x2(a):
    """
    a - symmetric 2x2 matrix
    
    returns:
        e - sorted eigenvalues
        v - corresponding right eigenvectors
    """
    n = a.shape[2:]
    
    e = zeros((2,)+n)
    v = zeros((2,2)+n)
    
    det = ((a[0,0]-a[1,1])**2+4*a[0,1]**2)**.5
    trc = a[0,0]+a[1,1]
    
    
    e[0] = .5*(trc+det)
    e[1] = .5*(trc-det)
    
    e = sort(e, axis = 0)
    
    denom = (a[0,1]**2+(e-a[0,0])**2)**.5
    
    
    v[0] = e-a[0,1]/denom
    v[1] = e-a[0,0]/denom
    
    return (e,v)
    
    
class mask2graph(object):
    def __init__(self, mask, nodevals = None):
        self.m = mask.copy()
        self.g = networkx.Graph()
        
        
        if nodevals == None:
            self.nodevals = ones_like(self.m)
        else:
            self.nodevals = nodevals.copy()
            
        #add nodes to the network
        for i in flatnonzero(self.m):
            self.g.add_node(unravel_index(i, self.m.shape), {'val': self.nodevals.ravel()[i]})
            
        
    def __call__(self, neighborFun = None):
        idxs = self.g.nodes()
        
        if neighborFun == None:
            if array(idxs).shape[1]==2:
                nidxs = self.neigh2d(array(idxs))
            if array(idxs).shape[1]==3:
                nidxs = self.neigh3d(array(idxs))
        else:
            nidxs = neighborFun(self, array(idxs))
        
        for i,x in ndenumerate(zeros(nidxs.shape[:-1])):
            nid = tuple(nidxs[i+(slice(None, None),)])
            id = idxs[i[0]]
            
            try:
                if self.m[nid]:
                    self.g.add_edge(nid, id, {'weight':sum((float64(nid)-float64(id))**2)**.5})
            except IndexError:
                continue
                
        return self.g
    
    def neigh2d(self, idxs):
        """
        idxs - #indices x 2, array of 2-d indices to retrieve the neighbors for
        
        returns: 
            nidxs - #indices x #neighbors x2 array of 2d indices of the neighbors, nidxs[i,j] gives 
                a 2d index of point i's j'th neighbor
        """
        #build stencils to get each neighboring index
        stencils = array([])
        for i in range(-1,2):
            for j in range(-1,2):
                if i==0 and j==0:
                    continue 
                try:
                    stencils = r_['0,2', stencils,  [i,j] ]
                except ValueError:
                    stencils = r_['0,1', stencils,  [i,j] ]
        
        #in 2d use 8 connected pixels
        assert stencils.shape[0]==8
        
        #add stencils to the indices to get the neighbors
        
        nidxs = int64(idxs[:,newaxis,...]+stencils)
        
        return nidxs
        
    def neigh3d(self, idxs):
        """
        idxs - #indices x 3, array of 3-d indices to retrieve the neighbors for
        
        returns: 
            nidxs - #indices x #neighbors x 3 array of 3d indices of the neighbors, nidxs[i,j] gives 
                a 3d index of point i's j'th neighbor
        """
        #build stencils to get each neighboring index
        stencils = array([])
        for i in range(-1,2):
            for j in range(-1,2):
                for k in range(-1,2):
                    if i==0 and j==0 and k==0:
                        continue 
                    try:
                        stencils = r_['0,2', stencils,  [i,j,k] ]
                    except ValueError:
                        stencils = r_['0,1', stencils,  [i,j,k] ]
        
        #in 2d use 8 connected pixels
        assert stencils.shape[0]==26
        
        #add stencils to the indices to get the neighbors
        
        nidxs = int64(idxs[:,newaxis,...]+stencils)
        
        return nidxs
        
def neigh3d_124connected(self, idxs):
        """
        idxs - #indices x 3, array of 3-d indices to retrieve the neighbors for
        
        returns: 
            nidxs - #indices x #neighbors x 3 array of 3d indices of the neighbors, nidxs[i,j] gives 
                a 3d index of point i's j'th neighbor
        """
        #build stencils to get each neighboring index
        stencils = array([])
        for i in range(-2,3):
            for j in range(-2,3):
                for k in range(-2,3):
                    if i==0 and j==0 and k==0:
                        continue 
                    try:
                        stencils = r_['0,2', stencils,  [i,j,k] ]
                    except ValueError:
                        stencils = r_['0,1', stencils,  [i,j,k] ]
        
        #in 2d use 8 connected pixels
        try:
            assert stencils.shape[0]==124
        except AssertionError:
            pdb.set_trace()
            
        
        #add stencils to the indices to get the neighbors
        #pdb.set_trace()
        nidxs = int64(idxs[:,newaxis,...]+stencils)
        
        return nidxs
def neigh3d_26connected(self, idxs):
        """
        idxs - #indices x 3, array of 3-d indices to retrieve the neighbors for
        
        returns: 
            nidxs - #indices x #neighbors x 3 array of 3d indices of the neighbors, nidxs[i,j] gives 
                a 3d index of point i's j'th neighbor
        """
        #build stencils to get each neighboring index
        stencils = array([])
        for i in range(-1,2):
            for j in range(-1,2):
                for k in range(-1,2):
                    if i==0 and j==0 and k==0:
                        continue 
                    try:
                        stencils = r_['0,2', stencils,  [i,j,k] ]
                    except ValueError:
                        stencils = r_['0,1', stencils,  [i,j,k] ]
        
        #in 2d use 8 connected pixels
        try:
            assert stencils.shape[0]==26
        except AssertionError:
            pdb.set_trace()
            
        
        #add stencils to the indices to get the neighbors
        #pdb.set_trace()
        nidxs = int64(idxs[:,newaxis,...]+stencils)
        
        return nidxs        

class abstractScaleSpaceRepresentation(object):
    """
    represent an image in scale space 
    """
    def __init__(self, im, #image to operate on
                 nscales = 50, # number of scales to compute
                 smoothFun = smoothg, #smoothing function must accept arguments (image, scaleParameter, derivativeOrderTuple)
                 maxFeatureSize = .25, #largest scale in fraction of image size
                 ):
        """
        set the image im, and compute the scaled versions
        nscales - number of scales to compute
        """
        self.nscales = nscales
        self.maxFeatureSize = maxFeatureSize
        self.smooth = smoothFun
        self.im = im
        
        #container for the derivatives
        self.D = [None, None, None]
        self.scaleSpace()
    
    def updateDerivatives(self):
        self.gradient()
        self.hessian()
        self.thirdOrderDerivatives()
    def scaleNorm(self, t, gamma = .5):
        """
        return the normalizing factor for a derivative at the scale t
        """
        return t**gamma
    
    def scaleSpace(self, im = None):
        if im !=None:
            self.im = im
        
        self.n = self.im.shape
        self.scales = linspace(.75, (self.maxFeatureSize*max(self.n)), self.nscales)
        self.ss = zeros((self.nscales,)+self.n  )
        
        for i in range(self.nscales):
            self.ss[i] = self.smooth(self.im, self.scales[i])
            
        return self.ss
    
    def gradient(self):
        """
        compute the gradient of the scale space images
        """
        pass
        
    def hessian(self):
        """
        compute the hessian for the scale space images
        """
        pass
        
    def thirdOrderDerivatives(self):
        """
        compute the third order partial derivatives of the scale space images
        """
        pass
        
    def dget(self, idx = [1,0]):
        """
        return the derivative of order len(idx),
        which has been differentiated idx[i] times in the ith direction
        """
        pass
            
            
            
            
        
class scaleSpace2d(abstractScaleSpaceRepresentation):
    def __init__(self,im, **kwargs):
        assert len(im.shape)==2, "im must be 2d!"
        abstractScaleSpaceRepresentation.__init__(self, im,**kwargs)
    
    def updateDerivatives(self):
        self.hessian()
    
    
    def gradient(self, normed = True):
        k = ones_like(self.scales)
        if normed:
            k = self.scaleNorm(self.scales)
        
        
        self.D[0] = zeros((2,)+self.ss.shape)
        order = ([1,0],[0,1])
        cnt = 0
        for ord in order:
            for i in range(self.nscales):
                self.D[0][cnt,i] = k[i]*self.smooth(self.im, self.scales[i], ord = ord)
            cnt+=1
    
    def hessian(self, normed = True):
        k = ones_like(self.scales)
        if normed:
            k = self.scaleNorm(self.scales)**2
        self.D[1] = zeros((3,)+self.ss.shape)
        order = ([2,0],[1,1],[0,2])
        cnt = 0
        for ord in order:
            #print "ord is"+ str(ord)+ "cnt is %d"%cnt
            for i in range(self.nscales):
                self.D[1][cnt,i] = k[i]*self.smooth(self.im, self.scales[i], ord = ord)
 
            cnt+=1
            
    
    def thirdOrderDerivatives(self, normed = True):
        k = ones_like(self.scales)
        if normed:
            k = self.scaleNorm(self.scales)**3
            
        self.D[2] = zeros((4,)+self.ss.shape)
        order = ([3,0],[2,1],[1,2],[0,3])
        cnt = 0
        for ord in order:
            for i in range(self.nscales):
                self.D[2][cnt,i] = k[i]*self.smooth(self.im, self.scales[i], ord = ord)
            cnt+=1
            
            
    def dget(self, idx = [1,0]):
        n = sum(array(idx))       
        if n==1:
            order = (str([1,0]),str([0,1]))
            d = dict(zip(order,range(2)))
        if n==2:
            order = (str([2,0]),str([1,1]),str([0,2]))
            d = dict(zip(order,range(3)))
        if n==3:
            order = (str([3,0]),str([2,1]),str([1,2]),str([0,3]))
            d = dict(zip(order,range(4)))
            
        return self.D[n-1][d[str(idx)]].copy()
    
    def show(self, figname = 'contours in scale space'):
        f = mlm.cont4(s = self.ss, axeLab = ['t','x','y'], f = mlm.fig(figname))
        return f

class bareBonesScaleSpace2d(scaleSpace2d):
    """
    a stripped down version of scaleSpace2d containing only the essentials
    """
    def __init__(self, ss2d):
        """
        copy the appropriate data members from ss2d, leaving the rest behind
        """
        self.ss = float32(ss2d.ss.copy())
        self.scales = ss2d.scales.copy()
        self.nscales = self.scales.shape[0]
        self.maxFeatureSize = 0+ss2d.maxFeatureSize
        self.im = ss2d.im.copy()
        

    