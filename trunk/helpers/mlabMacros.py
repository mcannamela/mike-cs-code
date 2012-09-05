import numpy as np
from numpy import *
try:
    from mayavi import mlab as ml
except ImportError:
    from enthought.mayavi import mlab as ml
import pdb
#some macros to aid plotting, range, scale, quiver, scatter, figure

#support functions
def rg(x):
    r = float64(array([np.min(x.ravel()), np.max(x.ravel())]))
    #assert all(r!=0), 'dimension %d with zero range detected'%(find(r==0)[0])
    return r
    
def rgs(x,y,z):
    return array([rg(x), rg(y), rg(z)]).ravel()
def sc(x):
    return float64(x)/diff(rg(x.ravel()))
def scLike(x,master):
   return float64(x)/diff(rg(master.ravel()))
   
ex = (-1.,1.,  -1.,1.,  -1.,1.)

def scaleExtent(slaveRg, masterRg, masterEx = ex):
    """
    return the extent needed to reconcile rg into extent masterEx if masterEx spans masterRg
    """
    ev = slice(0,None, 2)
    odd = slice(1,None,2)
    ME = array(masterEx)
    deltaEx = ME[odd]-ME[ev]
    deltaRg = masterRg[odd]-masterRg[ev]
    
    slaveEx = zeros_like(ME)
    slaveEx[ev] = (slaveRg[ev]/deltaRg)*deltaEx + ME[ev]
    slaveEx[odd] = (slaveRg[odd]/deltaRg)*deltaEx + ME[ev]
    
    print slaveEx
    return tuple(slaveEx)
    
gridSlice = lambda x_:slice(x_[0],x_[-1],1j*x_.shape[0])

#creation and detailing funcitons
out = lambda:                                     ml.outline(extent = ex)
def fig(name = 'foo'):
    try:
        f = ml.figure(figure = name, size = (1024, 720))
    except TypeError:
        f = ml.figure(name = name, size = (1024, 720))
    return f
    
def axe(labs, rgs, f):
    """
    add an axes object to the plot
    labs - list of axis labels
    rgs - array of axis ranges, may be produced by rgs
    f - figure to which to add the axes
    """
    

    try:
        ml.axes(xlabel = labs[0], ylabel = labs[1], zlabel = labs[2],
                ranges =rgs, extent = ex, nb_labels = 5, figure = f)
    except:
        ml.axes(xlabel = labs[0], ylabel = labs[1], zlabel = labs[2],
                ranges =rgs, extent = ex, figure = f)
    
    
        


#plotting functions

def quiv(X,U,f=None ,axeLab = ['', '', ''], clf = True):
    if clf:
        ml.clf(f)
    
    ml.quiver3d(sc(X[0]), sc(X[1]), sc(X[2]),
                           sc(U[0]), sc(U[1]), sc(U[2]),
                           figure = f,  scale_factor = .01, scale_mode = 'none',mode = 'cone')
                                
    
    
    axe(axeLab, rgs(X[0],X[1],X[2]),f)
    out()

def scat(X, f =None,axeLab = ['', '', ''], clf = True, scale = 'none'):
    """
    X.shape[0] = 3 or 4
    """
            
    if clf:
        ml.clf(f)
    
    if X.shape[0]==3:
        ml.points3d(sc(X[0].ravel()), sc(X[1].ravel()), sc(X[2].ravel()), figure = f, scale_mode =scale, scale_factor = .02, extent = ex)
    elif X.shape[0] ==4:
        p = ml.points3d(sc(X[0].ravel()), sc(X[1].ravel()), sc(X[2].ravel()),sc(X[3].ravel()), figure = f, scale_mode =scale, scale_factor = .02, extent = ex)
        ml.colorbar(p)
    
    axe(axeLab, rgs(X[0],X[1],X[2]),f)
    out()
    return f
    
def plot3(X, f =None,axeLab = ['', '', ''], clf = True):
    """
    X.shape[0] = 3 or 4
    """
            
    if clf:
        ml.clf(f)
    
    for i in range(X.shape[0]):
        assert diff(rg)!=0, 'range of axis %d is zero!'%(i)
            
            
    if X.shape[0]==3:
        ml.plot3d(sc(X[0]), sc(X[1]), sc(X[2]), figure = f,  extent = ex)
    elif X.shape[0] ==4:
        p = ml.plot3d(sc(X[0]), sc(X[1]), sc(X[2]),sc(X[3]), figure = f,  extent = ex)
        ml.colorbar(p)
    
    axe(axeLab, rgs(X[0],X[1],X[2]),f)
    out()
    return f
        
def surf3(x=None,y=None,z=None, f =None,axeLab = ['', '', ''], clf = True):
    """
    3 argument surface plot
    """
    if x==None:
        x = arange(z.shape[0])
    if y==None:
        y = arange(z.shape[1])
    if z==None:
        z = arange(x.shape[0]*y.shape[0])
    if clf:
        ml.clf(f)
    if z.ndim==3:
        for i in range(z.shape[-1]):
            ml.surf(sc(x),sc(y),scLike(z[...,i], z[...,0]), extent = ex, figure =f)
            axe(axeLab, rgs(x,y,z[...,0]),f)
    else:
        ml.surf(sc(x),sc(y),sc(z), extent = ex, figure =f)
        axe(axeLab, rgs(x,y,z),f)
    out()
    return f

def surfLike(x = None,y = None,z = None,X = None, f = None):
    """
    plot z(x,y) with the same scaling as X[2](X[0], X[1])
    """
    try:
        switch = (X.ndim ==2)
    except AttributeError:
        switch = False
        
        if switch:
            if x==None:
                x = arange(X.shape[0])
            if y==None:
                y = arange(X.shape[1])
            
            ml.surf(scLike(x,x), scLike(y,y),scLike(z, X), extent = scaleExtent(rgs(x,y,z), rgs(x,y,X)),  figure = f)
        else:
            print "switch is false"
            if x==None:
                x = X[0]
            if y==None:
                y = X[1]
            if z==None:
                z = arange(x.shape[0]*y.shape[0])
            ml.surf(scLike(x,X[0]), scLike(y,X[1]), scLike(z,X[2]),extent = scaleExtent(rgs(x,y,z), rgs(X[0],X[1],X[2])), figure = f)
    
        
        
    return f
    
    
def mesh3(x=None,y=None,z=None, f =None, axeLab = ['', '', ''],clf = True):
    """
    3 argument surface plot
    """
    if x==None:
        x = arange(z.shape[0])
    if y==None:
        y = arange(z.shape[1])
    if z==None:
        z = arange(x.shape[0]*y.shape[0])
    
    
    if x.ndim==1:
        x = tile(x[...,newaxis], (1,z.shape[1]))
    if y.ndim==1:
        y = tile(y[newaxis,...], (z.shape[0], 1) )
        
    if clf:
        ml.clf(f)
    ml.mesh(sc(x),sc(y),sc(z), extent = ex, figure =f)
    axe(axeLab, rgs(x,y,z),f)
    out()
    return f

def cont3(x=None,y=None,z=None, contours = 20, f =None,axeLab = ['', '', ''], clf = True):
    """
    3 argument surface plot
    """
    if x==None:
        x = arange(z.shape[0])
    if y==None:
        y = arange(z.shape[1])
    if z==None:
        z = arange(x.shape[0]*y.shape[0])
    if clf:
        ml.clf(f)
    if z.ndim==3:
        for i in range(z.shape[-1]):
            ml.contour_surf(sc(x),sc(y),scLike(z[...,i], z[...,0]), 
                contours = contours, extent = ex, figure =f)
            axe(axeLab, rgs(x,y,z[...,0]),f)
    else:
        ml.contour_surf(sc(x),sc(y),sc(z), 
                contours = contours, extent = ex, figure =f)
        axe(axeLab, rgs(x,y,z),f)
    out()
    return f 

def cont4(x=None,y=None,z=None,s = None,c=10, opacity = .5,axeLab = ['', '', ''], f =None):
    """
    4 argument contour plot
    x, y, z should be arrays
    """
    if x==None:
        x = arange(s.shape[0])
    if y==None:
        y = arange(s.shape[1])
    if z==None:
        z = arange(s.shape[2])
    if s==None:
        s = arange(x.shape[0]*y.shape[0]*z.shape[0]).reshape((x.shape[0],y.shape[0],z.shape[0]))
    
    g = mgrid[gridSlice(sc(x)),gridSlice(sc(y)),gridSlice(sc(z))]
    ml.contour3d(g[0],g[1],g[2],s, contours = c,opacity = opacity, transparent = True, extent = ex, figure =f)
    axe(axeLab, rgs(x,y,z),f)
    out()
    ml.colorbar()
    return f

