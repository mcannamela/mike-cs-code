#here we implement a class that models a spherical particle in a fluid field
import helpers as hp
import interpolators as pol
import solvers as sol

import gasMixture as gm
import cPickle as pkl

import scipy
from numpy import *
from pylab import *
#import mlabMacros as mlm
import time
import pdb 
import cProfile
import pstats
import os
from blockMeshTools import blockMesh
from gasModel import gas, binaryGasMixture
import gc

class plasmaField(object):
    """
    class to represent a finished jet object
    provides spatial interpolation for 
    """
    def __init__(self, simName = "cur-700_flo-50_gas-argon68helium32_make-SG100", 
                    fastOn = True, 
                    velocityScaling = 1, 
                    temperatureScaling = 1):
        """
        initialize with a pickled summary file
        """
        try:
            f=open(simName+".pkl",'rb')
        except IOError:
            try:
                f = open(os.path.join(os.pardir, 'pyrticle',simName+'.pkl'),'rb')
            except IOError:
                try:
                    f = open(os.path.join('D:\\CODE', 'pyrticle',simName+'.pkl'),'rb')
                except IOError: 
                    f = open(os.path.join('/media/ubuntuData/CODE', 'pyrticle',simName+'.pkl'),'rb')
                    
        #data dictionary
        self.D = pkl.load(f)
        f.close()
        self.reconstructGas()
        self.D['velocity'][self.D['velocity']<0]=0
        
        self.velocityScaling = velocityScaling
        self.temperatureScaling = temperatureScaling
        
        self.setCallbacks(fastOn)
        
    
    def __getstate__(self):
        """
        need to exclude all lambda objects
        """
        black =['polDict', 'dTz', 'dVz', 'dTr', 'dVr']
        S = []
        K = []
        for k in self.__dict__.keys():
            if k in black:
                continue
            K+=[k]
            S+=[self.__dict__[k]]
            
        return dict(zip(K,S))
        
    def __setstate__(self, D):
        """
        gonna have to call the methods that reconstruct the lambdas, which cannot be pickled
        """
        self.__dict__=D.copy()
        self.setCallbacks(fastOn = True)

        
    def reconstructGas(self):
        """
        parse out the gas names and volume fractions and reconstruct the gas for use in property fetching
        """
        validGasNames = ['argon', 'hydrogen', 'helium', 'nitrogen']
        gasNames = []
        volumeFractions = []
        for k in validGasNames:
            #if it fails then that gas was not in the mix
            try:
                volumeFractions+= [self.D[k]]
                gasNames+= [k]
            except:
                continue
        self.gas = gm.gasMixture(gasNames,float32(volumeFractions))

        
    def setCallbacks(self, fastOn = True):
        """
        set up interpolating callbacks that give the state of the plasma field at a given x,y,z
        """
        extendOne = lambda x_: hstack((  transpose(atleast_2d(x_[:,0]))  ,x_))
        if fastOn:
            #zeroth order interpolation will be fast, but need a fine mesh for accuracy
            master = zeros((100,100))
            T = float32(hp.matchImageSize(master, extendOne(self.D['temperature'])))
            V =  float32(hp.matchImageSize(master,extendOne(self.D['velocity'])))
        else:
            T = extendOne(self.D['temperature'])
            V = extendOne(self.D['velocity'])
        
        #temp variables for improved readability
        z = self.D['axialCoord']
        if self.D['radialCoord'][0]>0:
            r = hstack((zeros(1), self.D['radialCoord']))
        else:
            r = self.D['radialCoord']
                 
        T[T<0]=0
        V[V<0]=0
        
        #convenient and arbitrary scaling 
        T*= self.temperatureScaling
        V*= self.velocityScaling
        
        #need the gradient to see what's ahead!
        self.dT = float32(gradient(T, diff(z[-2:]),diff(r[-2:])) )
        self.dV =  float32(gradient(V, diff(z[-2:]),diff(r[-2:])) )
        
        #axis limits
        self.mn = float32([np.min(z), np.min(r)])
        self.mx = float32([np.max(z), np.max(r)])
        
        #set fast and slow interpolators here
        if fastOn:
            I = pol.zerothOrderInterpolator
        else: 
            I = pol.cubicInterpolator
        
        #initialize interpolators for the gradients 
        self.dTr = I(self.mn, self.mx,self.dT[1])
        self.dTz = I(self.mn, self.mx,self.dT[0])
        
        self.dVr = I(self.mn, self.mx,self.dV[1])
        self.dVz = I(self.mn, self.mx,self.dV[0])
        
        
        
        #build a dictionary of interpolators for the field variables and properties
        self.polDict = dict(zip(['temperature', 'velocity', 'conductivity', 
                                    'density', 'enthalpy', 'specificHeat', 'viscosity'],
                                    [I(self.mn, self.mx,T),
                                    I(self.mn, self.mx,V),
                                    I(self.mn, self.mx,self.gas.Conductivity(T)),
                                    I(self.mn, self.mx,self.gas.Density(T)),
                                    I(self.mn, self.mx,self.gas.Enthalpy(T)),
                                    I(self.mn, self.mx,self.gas.SpecificHeat(T)),
                                    I(self.mn, self.mx,self.gas.Viscosity(T))]))
        
        
    
    def temperature(self,x):
        return self.polDict['temperature'](self.cart2cyl(x))
    def velocity(self,x):
        return self.polDict['velocity'](self.cart2cyl(x))
    def conductivity(self,x):
        return self.polDict['conductivity'](self.cart2cyl(x))
    def density(self,x):
        return self.polDict['density'](self.cart2cyl(x))
    def enthalpy(self,x):
        return self.polDict['enthalpy'](self.cart2cyl(x))
    def specificHeat(self,x):
        return self.polDict['specificHeat'](self.cart2cyl(x))
    def viscosity(self,x):
        return self.polDict['viscosity'](self.cart2cyl(x))
     
        
    def temperatureGradient(self, x):
        """
        return a 2 vector of the temperature gradient at the point x
        """
        y = self.cart2cyl(x)
        
        return hstack((self.dTr(y), self.dTz(y)))
    
    def velocityGradient(self, x):
        """
        return a 2 vector of the temperature gradient at the point x
        """
        y = self.cart2cyl(x)
        
        return hstack((self.dVr(y), self.dVz(y)))
        
    def cart2cyl(self, X):    
        """
        convert the cartesian torch frame a radial coordinate
        X - rows are
            x, positive out the torch axis
            y, positive up
            z, x cross y
            
        output has rows
            x, positive out the torch axis
            r, radially away from x
            theta, ccw angle from old z axis
        """
        Y = zeros_like(X)
        
        Y[0] = X[0].copy()        
        Y[1] = sqrt(X[1]**2+X[2]**2)
        Y[2] = arctan2(X[1],X[2])
        
        #enforce bounds of the array
        
        if Y.ndim==1:
            Y[0] = min(Y[0],.99*self.mx[0])
            Y[1] = min(Y[1],.99*self.mx[1])
            Y[1] = max(Y[1], self.mn[1])
        else:
            Y[0][Y[0]>self.mx[0]] = .99*self.mx[0]
            Y[1][Y[1]>self.mx[1]] = .99*self.mx[1]
            Y[1][Y[1]<self.mn[1]] = self.mn[1]
        
        return Y 
        
    def showField(self, temperature = True, velocity = True):
        """
        make contour plots of the field
        """
        z = self.D['axialCoord']
        r = hstack((-flipud(self.D['radialCoord']),self.D['radialCoord']))
        T = vstack((flipud(transpose(self.D['temperature'])), transpose(self.D['temperature'])))
        try:
            V = vstack((flipud(transpose(self.D['velocity'])), transpose(self.D['velocity'])))
        except:
            V = vstack((flipud(transpose(self.D['axialVelocity'])), transpose(self.D['axialVelocity'])))
        figure(figsize = (20, 8))
        if temperature and velocity:
            subplot(2,1,1)
            contourf(z,r, T,20)
            colorbar(format = '%d')
            subplot(2,1,2)
            contourf(z,r,V,20)
            colorbar(format = '%d')
        elif temperature:
            contourf(z,r, T,20)
            colorbar(format = '%d')
        elif velocity:
            contourf(z,r,V,20)
            colorbar(format = '%d')
        
        
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::        
class dummyField(plasmaField):
    """
    make a silly plasma field for testing having linear velocity and temperature fields with a simple ficticious gas 
    """
    def __init__(self, vscale = 0, tscale = 1, fastOn = True):
        """
        """
        r = linspace(0,1,10)
        z = linspace(0,1,20)
        v = linspace(1,0,10)[:,newaxis]*ones_like(z)*vscale
        T = linspace(1,0,10)[:,newaxis]*ones_like(z)*tscale
        self.velocityScaling = 1
        self.temperatureScaling = 1
        
        
        self.gas = gm.gasMixture(['simpleGas'])
        self.D = dict(zip(['axialCoord', 'radialCoord', 'temperature', 'velocity'],[z,r,T,v]))
        self.setCallbacks(fastOn)
        

class lavaField(plasmaField):
    """
    make a plasma field from LAVA data
    """
    def __init__(self, simName = "lavaJet.pkl", 
                    simDict = None,
                    fastOn = True, 
                    velocityScaling = 1, 
                    temperatureScaling = 1,
                    gasPickle = None):
        self.simName = simName
        
        #use existing gas pickle if it exists
        if gasPickle==None and os.path.exists(self.pickleName()):
            gasPickle = self.pickleName()
        
        #if we are not passed a dictionary to use we'll ahve to unpickle one
        if simDict==None:
            
            try:
                
                if os.path.exists(simName):
                    f=open(simName,'rb')
                else:
                    f=open(simName+".pkl",'rb')
            except IOError:
                f = open(os.path.join(os.pardir, 'pyrticle',simName),'rb')
            
            #data dictionary
            self.D = pkl.load(f)
            f.close()
        else:
            self.D = simDict

        self.D['axialVelocity'][self.D['axialVelocity']<0]=0
        
        if np.max(self.D['axialCoord'])>1:
            self.D['axialCoord']*=1e-2
            self.D['radialCoord']*=1e-2
        self.reconstructGas(gasPickle)
                
        self.velocityScaling = velocityScaling
        self.temperatureScaling = temperatureScaling
        
        self.setCallbacks(fastOn)
    
    def pickleName(self):
        return self.simName[:-4]+'_gas.pkl'
    def reconstructGas(self,  gasPickle = 'lavaJet_gas.pkl'):
        """
        parse out the gas names and volume fractions and reconstruct the gas for use in property fetching
        """
        
        self.gasNames = ['argon', 'hydrogen', 'air']
        self.molFrac = r_['0,3', self.D['ArMols'], self.D['HMols'], self.D['airMols']]
        self.gas = gm.gasMixture(self.gasNames,None, molFractions = squeeze(self.molFrac[:,0,0]) )
        
        if gasPickle==None:
            print 'mixing gasses in all cells, this could take a while...'
            self.Conductivity = zeros_like(self.D['airMols'])
            self.Density= zeros_like(self.D['airMols'])
            self.Enthalpy= zeros_like(self.D['airMols'])
            self.SpecificHeat= zeros_like(self.D['airMols'])
            self.Viscosity = zeros_like(self.D['airMols'])
            
            T = self.D['temperature']
            print 'quantizing...'
            mf, codes, distortion, codebook = hp.vectorQuantize(self.molFrac, downsample = 50)
            print '...complete'
            cnt = 0
            for i in range(codebook.shape[0]):
                print 'completed %d of %d'%(i,codebook.shape[0])
                self.gas.blend(self.gasNames, None, squeeze(codebook[i]))
                idx = (codes == i)
                self.Conductivity[idx] = self.gas.Conductivity(T[idx])
                self.Density[idx] = self.gas.Density(T[idx])
                self.Enthalpy[idx] = self.gas.Enthalpy(T[idx])
                self.SpecificHeat[idx] = self.gas.SpecificHeat(T[idx])
                self.Viscosity[idx] = self.gas.Viscosity(T[idx])
                cnt+=1
    
            self.gas.blend(self.gasNames, None, molFractions = squeeze(self.molFrac[:,0,0]))
            with open(self.pickleName(), 'wb') as f:
                pkl.dump(dict(zip(['conductivity', 'density', 'enthalpy', 'specificHeat', 'viscosity'],
                                [self.Conductivity, self.Density, self.Enthalpy, self.SpecificHeat, self.Viscosity])),
                                f,-1)
            print 'mixing complete, gas has been pickled!'
            
        else:
            with open(gasPickle, 'rb') as f:
                gp = pkl.load(f)
            self.Conductivity = gp['conductivity']
            self.Density = gp['density']
            self.Enthalpy = gp['enthalpy']
            self.SpecificHeat = gp['specificHeat']
            self.Viscosity = gp['viscosity']
            
        

        
    def setCallbacks(self, fastOn):
        """
        set up interpolating callbacks that give the state of the plasma field at a given x,y,z
        """
        
        if fastOn:
            #zeroth order interpolation will be fast, but need a fine mesh for accuracy
            master = zeros((1000,1000))
            T = float32(hp.matchImageSize(master, (self.D['temperature'])))
            V =  float32(hp.matchImageSize(master,(self.D['axialVelocity'])))
        else:
            
            T = (self.D['temperature'])
            V = (self.D['axialVelocity'])
        
        #temp variables for improved readability
        z = self.D['axialCoord']
        r = self.D['radialCoord']
                 
        T[T<0]=0
        V[V<0]=0
        
        #convenient and arbitrary scaling 
        T*= self.temperatureScaling
        V*= self.velocityScaling
        
        #need the gradient to see what's ahead!
        self.dT = float32(gradient(T, diff(z[:2]),diff(r[:2])) )
        self.dV =  float32(gradient(V, diff(z[:2]),diff(r[:2])) )
        
        #axis limits
        self.mn = float32([np.min(z), np.min(r)])
        self.mx = float32([np.max(z), np.max(r)])
        
        #set fast and slow interpolators here
        if fastOn:
            I = pol.zerothOrderInterpolator
        else: 
            I = pol.cubicInterpolator
        
        #initialize interpolators for the gradients 
        self.dTr = I(self.mn, self.mx,self.dT[1])
        self.dTz = I(self.mn, self.mx,self.dT[0])
        
        self.dVr = I(self.mn, self.mx,self.dV[1])
        self.dVz = I(self.mn, self.mx,self.dV[0])
        
        
        #build a dictionary of interpolators for the field variables and properties
        self.polDict = dict(zip(['temperature', 'velocity', 'conductivity', 
                                    'density', 'enthalpy', 'specificHeat', 'viscosity'],
                                    [I(self.mn, self.mx,T),
                                    I(self.mn, self.mx,V),
                                    I(self.mn, self.mx,self.Conductivity),
                                    I(self.mn, self.mx,self.Density),
                                    I(self.mn, self.mx,self.Enthalpy),
                                    I(self.mn, self.mx,self.SpecificHeat),
                                    I(self.mn, self.mx,self.Viscosity)]))
                                    
class openFoamPlasma(object):
    def __init__(self, runFolder, tMin, tMax, xMin,xMax, G):
        """
        initialize with 
        
        runFolder  - path to the openFoam run
        tMin       - minimum time to include
        tMax       - maximum time to include
        xMin       - len 3 list of x,y,z minimum values to read between
        xMax       - len 3 list of x,y,z maximum values to read between
        G          - binaryGasMixture object representing the gas for the 
                        openFoam run
        """
        self.gas = G
        
        R = openFoamJetReader(runFolder)
        self.T, self.U, self.f, self.X , self.t = R(tMin, tMax, xMin, xMax)

        self.t-=self.t[0]
        
        self.gmT = float32([float32(unevenGradientMagnitude(self.X, T)) for T in self.T])
        
        self.polNames = ['temperature', 'xVelocity', 'yVelocity', 'zVelocity', 
                    'plasmaFraction', 'temperatureGradient']
        polList = [self.T]+ [self.U[:,i,...] for i in range(3)]+[self.f, self.gmT]
        
        self.pols = [unevenSpaceTimeInterpolator(p, self.X, self.t) 
                           for p in polList ]
        
        self.polDict = dict(zip(self.polNames, self.pols))
        del R
        
        
    def __call__(self, x, t):
        """
        return plasma properties at the specified place and time
        """
        U = zeros(3)
        T, U[0], U[1], U[2], f, gmT = [p(x,t) for p in self.pols] 
        
        rho, mu, k, Cp, M = self.gas(T,f)

        return (T, U, f,rho, mu, k, Cp, M)
    def setCallbacks(self, fastOn = True):
        pass
        
    def temperatureGradient(self, x, t):
        """
        return the gradient of the temperature field at the specified x and t
        """
        return self.polDict['temperatureGradient'](x,t)
        
    def temperature(self,x,t):
        return self.polDict['temperature'](x,t)
    def velocity(self,x,t):
        return self.polDict['velocity'](x,t)
    def conductivity(self,x,t):
        return self.polDict['conductivity'](x,t)
    def density(self,x,t):
        return self.polDict['density'](x,t)
    def enthalpy(self,x,t):
        return self.polDict['enthalpy'](x,t)
    def specificHeat(self,x,t):
        return self.polDict['specificHeat'](x,t)
    def viscosity(self,x,t):
        return self.polDict['viscosity'](x,t)
        
    def showField(self, temperature = True, velocity = True, planeNormal = 'y', t = None, withColor = True):
        """
        make contour plots of the field
        """
        if t==None:
            t = self.t[-1]/2            
        
        idx = argmin(abs(self.t-t))

        if planeNormal == 'y':
            r = self.X[1]
            T = self.T[idx][:,len(self.X[1])/2, :]
            V = self.U[idx,0][:,len(self.X[1])/2, :]
        elif planeNormal == 'z':
            r = self.X[2]
            r = self.X[1]
            T = self.T[idx][:, :, len(self.X[2])/2,]
            V = self.U[idx,0][:, :, len(self.X[2])/2]
            
        z = self.X[0]
        

        figure(figsize = (20, 8))
        
        if withColor:
            if temperature and velocity:
                subplot(2,1,1)
                contourf(z,r, T.T,15, cmap = cm.Reds )
                colorbar(format = '%d', cmap = cm.Reds)
                subplot(2,1,2)
                contourf(z,r,V.T,15, cmap = cm.Reds)
                colorbar(format = '%d')
            elif temperature:
                C = contourf(z,r, T.T,15, cmap = cm.Reds, linestyles = 'solid')
                #clabel(C, fontsize = 12, fmt = '%d', inline = 1, colors = 'k')
                colorbar(format = '%d')
                
            elif velocity:
                contourf(z,r,V.T,15, cmap = cm.Reds)
                colorbar(format = '%d')
        else:
            
            if temperature and velocity:
                subplot(2,1,1)
                C1 = contour(z,r, T.T,10, colors = 'k')
                clabel(C, fontsize = 12, fmt = '%d', inline = 1)
                
                subplot(2,1,2)
                C2 = contour(z,r,V.T,10, colors = 'k')
                clabel(C, fontsize = 12, fmt = '%d', inline = 1)
                
            elif temperature:
                C = contour(z,r, T.T,10, colors = 'k')
                clabel(C, fontsize = 12, fmt = '%d', inline = 1)
                
            elif velocity:
                C = contour(z,r,V.T,10, colors = 'k')
                clabel(C, fontsize = 12, fmt = '%d', inline = 1)
           
                
            
    def showFieldC(self, temperature = True, velocity = True, planeNormal = 'y', t = None):
        """
        make contour plots of the field
        """
        
        if planeNormal == 'y':
            r = self.X[1]
            T = self.T[0][:,len(self.X[1])/2, :]
            V = self.U[0,0][:,len(self.X[1])/2, :]
        elif planeNormal == 'z':
            r = self.X[2]
            r = self.X[1]
            T = self.T[0][:, :, len(self.X[2])/2,]
            V = self.U[0,0][:, :, len(self.X[2])/2]
            
        z = self.X[0]
        

        figure(figsize = (20, 8))
        
        if temperature and velocity:
            subplot(2,1,1)
            contourf(z,r, T.T,20)
            colorbar(format = '%d')
            subplot(2,1,2)
            contourf(z,r,V.T,20)
            colorbar(format = '%d')
        elif temperature:
            contourf(z,r, T.T,20)
            colorbar(format = '%d')
        elif velocity:
            contourf(z,r,V.T,20)
            colorbar(format = '%d')
            
class clonedOpenFoamPlasma(openFoamPlasma):
    def __init__(self, T, U, f, X , t , G):
        """
        initialize with 
        
        runFolder  - path to the openFoam run
        tMin       - minimum time to include
        tMax       - maximum time to include
        xMin       - len 3 list of x,y,z minimum values to read between
        xMax       - len 3 list of x,y,z maximum values to read between
        G          - binaryGasMixture object representing the gas for the 
                        openFoam run
        """
        self.gas = G
        
        self.T, self.U, self.f, self.X , self.t = (T,U,f,X,t) 

        self.t-=self.t[0]
        
        self.gmT = float32([float32(unevenGradientMagnitude(self.X, T)) for T in self.T])
        
        self.polNames = ['temperature', 'xVelocity', 'yVelocity', 'zVelocity', 
                    'plasmaFraction', 'temperatureGradient']
        polList = [self.T]+ [self.U[:,i,...] for i in range(3)]+[self.f, self.gmT]
        
        self.pols = [unevenSpaceTimeInterpolator(p, self.X, self.t) 
                           for p in polList ]
        
        self.polDict = dict(zip(self.polNames, self.pols))
        
    
class openFoamJetReader(object):
    def __init__(self, runFolder):
        self.runFolder  = runFolder
        self.BM = blockMesh(self.meshFolder())
        self.arrayNames = 'Temperature U f'.split()
        
    def __call__(self, tMin, tMax, xMin, xMax): 
        t,F = self.getTimes(tMin, tMax)
        iMin = self.BM.cc2sub(atleast_2d(xMin).T)
        iMax = self.BM.cc2sub(atleast_2d(xMax).T)
        
        slyTup = tuple([slice(int(iMin[i]), int(iMax[i])) for i in range(3)])
        
        bigArrs = [float32(zeros((len(t),)+tuple(iMax-iMin)))  
                    for nm in self.arrayNames]
        bigArrs[1] = float32(zeros((len(t),3)+tuple(iMax-iMin)))
        print 'reading timesteps'
        for i,f in enumerate(F):
            AA = self.readTime(f)
            for j in range(len(self.arrayNames)):
                if j==1:
                    bigArrs[j][i] = copy(AA[j][(slice(None, None),)+slyTup])
                else:
                    bigArrs[j][i] = AA[j][slyTup]
        return bigArrs+[[copy(self.BM.cc[i][st]) 
            for i,st in enumerate(slyTup)], copy(t)]
            
    
    def readTime(self, F):
        AA = []
        for nm in self.arrayNames:
            self.BM.setField(os.path.join(self.runFolder, F, nm))
            AA+=[float32(self.BM.fieldDict[nm])]
        
        return AA
            
    def meshFolder(self):
        return os.path.join(self.runFolder, 'constant', 'polyMesh')
    def timeStepFolders(self):
        t = []
        F = []
        for f in os.walk(self.runFolder).next()[1]:
            try:
                t+=[double(f)]
                F+=[f]
            except ValueError:
                continue
        a = argsort(array(t))
        
        return (array(t)[a], array(F,dtype = 'object')[a])
    def getTimes(self, tMin, tMax):
        t, F = self.timeStepFolders()
        idx = flatnonzero((t>=tMin)*(t<=tMax))
        return(t[idx],F[idx])

class unevenSpaceTimeInterpolator(object):
    def __init__(self, U, X, t):
        """
        U - nt x nx x ny x nz array of function values to interpolate
        X - len 3 list of x, y, z cell centers in resp. nx, ny, nz length arrays
        t - nt array of times
        """
        #set arguments as members        
        self.U = U
        self.X = X
        self.t = t
        
        
        #pre-cache limits on the bounds of space/time coordinates
        self.N = int32([len(x)*100 for x in self.X])
        self.Nminus = self.N-1
        self.I_ = int32(r_['0,2', zeros(3), int32([len(x)-1 for x in self.X])])
        self.idx__ = int32(zeros((2,3)))
        self.idx__[1] = int32(zeros(3))
        self.idx_ = int32(zeros((2,3)))
        self.idx_[1] = int32(self.Nminus)
        self.I = int32(zeros((2,3)))
        self.F = float32(zeros((2,3)))
        
        #pre-compute some indexing to select the 8 points for interpolation
        self.perms = [[(0, 0, 0), tuple(range(3))],
                      [(0, 0, 1), tuple(range(3))],
                      [(0, 1, 0), tuple(range(3))],
                      [(0, 1, 1), tuple(range(3))],
                      [(1, 0, 0), tuple(range(3))],
                      [(1, 0, 1), tuple(range(3))],
                      [(1, 1, 0), tuple(range(3))],
                      [(1, 1, 1), tuple(range(3))],
                       ]
        
        self.nt = len(t)        
        self.ntMinus = self.nt-1
        self.rt = t[-1]-t[0]
        self.tFac = self.ntMinus/self.rt
        
        
        self.rX = [x[-1]-x[0] for x in self.X]
        self.X0 = float64([x[0] for x in self.X])
        
        
        self.xFac = float64([ double(self.Nminus[i])/r 
                    for i, r in enumerate(self.rX)])
        
        
        self.Xn = [linspace(x[0], x[-1], len(x)*100) for x in self.X]
        self.iXn = [scipy.interpolate.interp1d(x, 
                    float64(arange(len(x))) )(self.Xn[i]) 
                        for i,x in enumerate(self.X)]
        
        self.fXn = [ixn - int32(ixn) for ixn in self.iXn]
        self.iXn = [int32(ixn) for ixn in self.iXn]
    
    def __call__(self, X,t):
        """
        interpolated the function at point X, time t
        X - len 3 array of x,y,z values of the point to interpolate
        t - scalar, time at which to interpolate
        """
        I0,fI0 = self.getIdx(X)
        i0,fi0 = self.i_t(t)
        
        i = int32([min([max([i0, 0]), self.ntMinus] ),
             min([max([i0+1, 0]), self.ntMinus] )])
             
        self.I[0] = I0
        self.I[1] = self.boundI(int32(I0)+1)
        
        self.F[1] = fI0
        self.F[0] = 1-float64(fI0)
        
        f  = float32(zeros(8))
        U0 = float32(zeros(8))
        U1 = float32(zeros(8))
        
        for j,p in enumerate(self.perms):
            f[j]  = float32(prod(self.F[p]))
            U0[j] = self.U[i[0]][tuple(self.I[p])]
            U1[j] = self.U[i[1]][tuple(self.I[p])]
        
        U = sum(U0*f)*(1-fi0)+sum(U1*f)*fi0

        return U
    
    def boundIdx(self, idx):
        self.idx__[0] = idx
        self.idx_[0] = amax(self.idx__, axis = 0)
        return amin(self.idx_, axis = 0)
    def boundI(self, I):
        self.I_[0] = I
        return amin(self.I_, axis = 0)
   
    def getIdx(self, X):
        """
        convert a point in space to a cell index in the array
        X - len 3 array giving a location in x,y,z space
        """
        #idx = [(i,j) 
        #    for i, j in enumerate(round( (X-self.X0)*self.xFac ))] 
        
        idx = self.boundIdx(int32( (X-self.X0)*self.xFac ))
        
        iXn = [ixn[idx[i]] for i,ixn in enumerate(self.iXn)]
        fXn = [fxn[idx[i]] for i,fxn in enumerate(self.fXn)]
                          
        return [iXn, fXn]
    
    def i_t(self, t):
        """
        compute the index for the time axis, which is regularly spaced
        t - double, the time at which to interpolate
        """
        f = (t-self.t[0])*self.tFac
        i  = floor(f)
        return (i, f-i)

def unevenGradientMagnitude(X, F):
    """
    return the approximate gradient magnitude of the 3d array F
    """
    
    sly = slice(None, -1)
    slyp = slice(1, None)
    
    Gz =  ( (F[..., slyp]-F[..., sly])/diff(X[2])[newaxis, newaxis, ...] )
    Gy =  ( (F[..., slyp,...]-F[..., sly,...])/diff(X[1])[newaxis, ..., newaxis] )
    Gx =  ( (F[slyp,...]-F[sly,...])/diff(X[0])[..., newaxis, newaxis] )
    
    Gy = r_['1,3', Gy, Gy[:,-1,:][:,newaxis,:]]
    Gz = r_['2,3', Gz, Gz[...,-1][...,newaxis]]
    Gx = r_['0,3', Gx, Gx[-1,...][newaxis,...]]
         
    return (Gz**2+Gy**2+Gx**2)**.5

def defaultOpenFoamField():
    runFolder = "/media/bigBoy/openFoamRuns/jet_sg100_ArPilot_40_500_DGMesh"
        
    with open('/media/raidArray/CODE/gasModel/argonPropertiesDict.pkl', 'rb') as f:
        DAr = pkl.load(f)
        DAr['enthalpy']-=DAr['enthalpy'][0]
    Ar = gas(40.0, DAr, 10000)


    with open('/media/raidArray/CODE/gasModel/airPropertiesDict.pkl', 'rb') as f:
        DA = pkl.load(f)
        DA['enthalpy']-=DA['enthalpy'][0]
    Air = gas(29.0, DA, 10000)
    
    G = binaryGasMixture(Ar, Air)
    
    tMin = .001
    tMax = .003
    xMin = array([0.0, -.02, -.02])
    xMax = array([.1, .02, .02])
    
    P = openFoamPlasma(runFolder, tMin, tMax, xMin,xMax, G)
    return P

if __name__== "__main__":
    fieldTestOn = False
    lavaFieldTestOn = False
    unevenInterpolatorTestOn = False
    openFoamJetReaderTestOn = False    
    openFoamPlasmaTestOn = False
    
    
    if lavaFieldTestOn :
        pf = lavaField(simName = "lavaJet.pkl", gasPickle = None)
        #pf = lavaField(gasPickle = 'lavaJet_gas.pkl')
#        mlm.surf3(z = pf.Density, f = mlm.fig('density') )
        
        
        
    
    if fieldTestOn:
        pf = plasmaField('..\\pyrticle\\cur-700_flo-50_gas-argon68helium32_make-SG100')
        dpf = dummyField()
        
        print 'plasma density'
        print pf.density(zeros(3))
        
        print 'dummy plasma density'
        print dpf.density(zeros(3)) 
    
    if unevenInterpolatorTestOn:
        a = arange(12.).reshape(3,4)
        A = r_['0, 3', a, a+12]

        A = r_['0,4', A, A+24]
        t = array([0, 1])
        x = array([2, 3])
        y = linspace(4,5,3)
        z = linspace(6,7,4)
        
        p = unevenSpaceTimeInterpolator(A, [x,y,z], t)
        print p(array([2.1, 4.1,6.1]), .25)
        
    if openFoamJetReaderTestOn:
        runFolder = "/media/bigBoy/openFoamRuns/jet_sg100_ArPilot_40_500_DGMesh"
        R = openFoamJetReader(runFolder)
        
        tMin = .001
        tMax = .0015
        xMin = array([0.0, -.001, -.001])
        xMax = array([.005, .001, .001])
        
        T,U, f, X,t = R(tMin, tMax, xMin,xMax)
        
#        figure()
#        subplot(1,3,1)
#        contourf(X[2],X[1],T[0][1,...],30); colorbar()
#        subplot(1,3,2)
#        contourf(X[2],X[0],T[0][:,18,...],30); colorbar()
#        subplot(1,3,3)
#        contourf(X[1],X[0],T[0][...,18],30); colorbar()
#        
#        figure()
#        subplot(1,3,1)
#        contourf(X[2],X[1],U[0,0][1,...],30); colorbar()
#        subplot(1,3,2)
#        contourf(X[2],X[0],U[0,0][:,18,...],30); colorbar()
#        subplot(1,3,3)
#        contourf(X[1],X[0],U[0,0][...,18],30); colorbar()
#        
#        figure()
#        subplot(1,3,1)
#        contourf(X[2],X[1],f[0][1,...],30); colorbar()
#        subplot(1,3,2)
#        contourf(X[2],X[0],f[0][:,18,...],30); colorbar()
#        subplot(1,3,3)
#        contourf(X[1],X[0],f[0][...,18],30); colorbar()
        
    if openFoamPlasmaTestOn:
        runFolder = "/media/bigBoy/openFoamRuns/jet_sg100_ArPilot_40_500_DGMesh"
        
        with open('/media/raidArray/CODE/gasModel/argonPropertiesDict.pkl', 'rb') as f:
            DAr = pkl.load(f)
            DAr['enthalpy']-=DAr['enthalpy'][0]
        Ar = gas(40.0, DAr, 10000)

    
        with open('/media/raidArray/CODE/gasModel/airPropertiesDict.pkl', 'rb') as f:
            DA = pkl.load(f)
            DA['enthalpy']-=DA['enthalpy'][0]
        Air = gas(29.0, DA, 10000)
        
        G = binaryGasMixture(Ar, Air)
        
        tMin = .00098
        tMax = .00104
        xMin = array([0.0, -.01, -.01])
        xMax = array([.03, .01, .01])
        
        P = openFoamPlasma(runFolder, tMin, tMax, xMin,xMax, G)
        T = P.T
        X = P.X
        t = P.t
        p = unevenSpaceTimeInterpolator(T,X,t)
        p(array([0,0,0]), 0)
        
#        for k in Ar.__dict__.keys():
#            try:
#                blah = pkl.dumps(Ar.__dict__[k], -1)
#            except pkl.PicklingError:
#                print k
    
    pickleDefaultOpenFoamField = True
    if pickleDefaultOpenFoamField:
        P = defaultOpenFoamField()
        with open('defaultOpenFoamField.pkl','wb') as f:
            pkl.dump(P, f, -1)
        
        
        