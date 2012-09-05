import numpy as np
from numpy import *
import scipy

import particleSolver as ps
import particle as part
import modelEnsemble as ens
import sensor as sens
from bigPickle import bigPickle

import pdb
import time
import mlabMacros as mlm



class basicParticleModel(ens.model):
    """
    a particle model consists of a particle object and a solver
    """
    def __init__(self, pVec, material = part.gm.defaultMaterial):
        """
        initialize with a parameter vector, (dIn, dOut, x0,v0)
        """
        self.PS = ps.particleSolver(part.particle(pVec[:2],
                                    T = pVec[2],
                                    thermalTimestep = 5e-5,
                                    kinematicTimestep = 5e-5,
                                    x0 = pVec[3:6], 
                                    v0 = pVec[6:9],
                                    nNodes = 30, 
                                    nStep = 1000,
                                    material = material))
        #self.PS.p.adaptiveTimestepping = False
        self.PS.p.thermalStep   = float32(75.)
        self.hasRun = False 
        
    def forwardSolve(self, plasmaField = part.plasmaField()):
        """
        given a plasma field, solve the particle
        """
        self.PS.p.inject(plasmaField)
        if self.hasRun:
            self.PS.reset()
        self.PS.solve()
        self.hasRun = True
        
    def parameterVector(self):
        """
        p = (innerDiameter, outerDiameter, 
                injectionPosition x,y,z,
                injectionVelocity x,y,z) 
                
        """
        return r[2*self.PS.p.mesh.xInt[[0,-1]],
                    self.PS.p.X[0],
                    self.PS.p.V[0]]
    
    def parameterKeys(self):
        return array(("innerDiameter outerDiameter temperature\
                        injPosX injPosY injPosZ  \
                        injVelX injVelY injVelZ").split())
    
    def summarize(self):
        """
        default observation behavior should return the usual solution
        """
        return dict(zip(['time', 'radius', 'temperature'],
                        self.PS.getSolnCoords()+(self.PS.getSoln(),)))
        
class emittingParticleModel(basicParticleModel):
    def __init__(self, pVec, material = part.gm.defaultMaterial):
        """
        initialize with a parameter vector, (dIn, dOut, x0,v0)
        """
        self.PS = ps.particleSolver(part.emittingParticle(pVec[:2],
                                    T = pVec[2],
                                    x0 = pVec[3:6], 
                                    v0 = pVec[6:9],
                                    nNodes = 30, 
                                    nStep = 200,
                                    material = material))
        
        self.hasRun = False 
    def parameterVector(self):
        """
        p = (innerDiameter, outerDiameter, 
                injectionPosition x,y,z,
                injectionVelocity x,y,z) 
                
        """
        return r[2*self.PS.p.mesh.xInt[[0,-1]],
                    self.PS.p.surfaceTemperature,
                    self.PS.p.X[0],
                    self.PS.p.V[0]]
    
    def parameterKeys(self):
        return array(("innerDiameter outerDiameter \
                        surfaceTemp \
                        injPosX injPosY injPosZ  \
                        injVelX injVelY injVelZ").split())
    
        
def basicGriddedParticleParameterGenerator( dMinMax = [10e-6,100e-6],
                                            vMinMax = [3, 18],
                                            x0 =array([0.001,.006,0.]),
                                            n = (10,10)):
    """
    make a regular mesh of particle parameters using solid 
    particles with a straight injection velocity
    
    dMinMax - minimum and maximum diameters in meters
    vMinMax - min and max injection velocity in m/s
    x0 - injection location
    n - number of grid points in the d, v directions
    """
    
    g = mgrid[dMinMax[0]:dMinMax[1]:n[0]*1j,
            vMinMax[0]:vMinMax[1]:n[1]*1j]
    
    G= r_['2,3,0', .01*g[0],g[0],0*g[0]+300,
                x0[newaxis,newaxis,:]*ones_like(g[0])[...,newaxis],
                zeros_like(g[1]),-g[1], zeros_like(g[1])]
    return G

#def basicRandomParticleParameterGenerator()

def basicGriddedEmittingParticleGenerator( dMinMax = [10e-6,100e-6],
                                            tMinMax = [2500, 4000],
                                            zMinMax = [0,3e-3],
                                            vMinMax = [100, 200],
                                            n = (10,1,1,5) ,
                                            nPixels = 8196,
                                            pixelSize = 5e-6,
                                            magnification = 1.5):
                                                
    g = mgrid[dMinMax[0]:dMinMax[1]:n[0]*1j,
                  tMinMax[0]:tMinMax[1]:n[1]*1j,
                  zMinMax[0]:zMinMax[1]:n[2]*1j,
                  vMinMax[0]:vMinMax[1]:n[3]*1j]

    n = g.shape[0]
    y = zeros_like(g[0])
    for i in range(n):
        if np.max(g[i])!=0:
            y+= g[i]/np.max(g[i])
    
    y-= np.min(y)
    y/= np.max(y)
    y-=.5
    
    y*= (.8*nPixels*pixelSize/magnification)
    
    G= r_['4,5,0', .01*g[0],g[0],
                         g[1], 
                        zeros_like(g[2]), y, g[2],
                        g[3], zeros_like(g[3]), zeros_like(g[3]) ]
    

    return G

def marginalRandomEmittingParticleGenerator(dStats = [50e-6, 15e-6], 
                                                                                tStats = [3200, 200], 
                                                                                yStats = [0, 4e-3], 
                                                                                zStats = [0,4e-3], 
                                                                                vStats = [150, 15], 
                                                                                N = 100):
    dist = scipy.stats.distributions
    
    d = dist.norm(loc = dStats[0], scale = dStats[1]).rvs(N)
    t = dist.norm(loc = tStats[0], scale = tStats[1]).rvs(N)
    y = dist.norm(loc = yStats[0], scale = yStats[1]).rvs(N)
    z = dist.norm(loc = zStats[0], scale = zStats[1]).rvs(N)
    v = dist.norm(loc = vStats[0], scale = vStats[1]).rvs(N)
    
    d[d<5e-6] = 10e-6
    t[t<0] = 1000
    v[v<0] = 100
    
    G = r_['1,2, 0',.01*d, d,t,zeros_like(y), y ,z , v, zeros_like(v), zeros_like(v)]
    
    return G
    

    
    
def dummyEmittingParticleGenerator(pVec = None):
    """
    pass an array of parameter vectors straight through 
    """
    return pVec
    
    
class basicParticleEnsemble(ens.modelEnsemble):
    """
    a collection of non-fancy particle objects
    """
    def __init__(self, parameterGeneratorFun, ensembleName):
        """
        set up the parameter vector generator
        """
        self.parameterGenerator = parameterGeneratorFun
        self.ensembleName = ensembleName
        
    def __getstate__(self):
        if 'ensembleName' not in self.__dict__.keys():
            D = self.__dict__
        else:
            K = []
            M = []
            for k in self.__dict__.keys():
                if k!='modelArray':
                    K+= [k]
                    M+= [self.__dict__[k]]
                else:
                    bigPickle().dump(self.modelArray, self.ensembleName+'.big.pkl')
            D = dict(zip(K,M))
        
        return [D, 
                    self.modelArray.flatten()[0].PS.p.plasmaField,  
                    part.gm.gasMixture(self.modelArray.flatten()[0].PS.p.materialName)]
                        
    def __setstate__(self, L):
        self.__dict__=L[0]
        if 'ensembleName' in self.__dict__.keys():
            print 'bigLoading...'
            self.modelArray = bigPickle().load(self.ensembleName+'.big.pkl')
        for idx, m in ndenumerate(self.modelArray):
            m.PS.p.reconstitute(L[1],L[2])
        
        
    def modelConstructor(self, parameterVector):
        """
        just use the model's constructor
        """
        return basicParticleModel(parameterVector)
    
    def ensembleSolve(self,**kwargs):
        start = time.time()
        ens.modelEnsemble.ensembleSolve(self,**kwargs)
        elapsed = time.time()-start
        print 'solved %d particles in %f sec: %f s/particle'%(self.modelArray.size, elapsed, elapsed/self.modelArray.size )
    
    def diameter(self):
        return self.parameterVectors[...,1]
    def temperature(self):
        return self.parameterVectors[...,2]
    
    def lateralPosition(self):
        return self.parameterVectors[...,5]
    def axialVelocity(self):
        return self.parameterVectors[...,6]
        
    def summary(self, plotOn = True):
        """
        make a bunch of surface plots if appropriate, and display 
        ensemble statistics
        """
        pv = self.parameterVectors
        pk = self.modelArray.flatten()[0].parameterKeys()
        if self.modelArray.ndim==2:
            #assume we meshed outer diameter and injection velocity
            g = r_['0,3',squeeze(1e6*pv[...,nonzero(pk=='outerDiameter')[0]]),
                                        squeeze(-pv[...,nonzero(pk=='injVelY')[0]])]
            axlab = ['diameter, microns', 'injection velocity, m/s']
            
            arrs = ['surfaceTemperature', 
                        'centerTemperature',
                        'meanTemperature',
                        'volumeFraction',
                        'diameterFraction',
                        'finalVelocity',
                        'enthalpy',
                        'surfaceEnthalpy',
                        'moltenFraction',
                        'moltenVolume',
                        'momentum',
                        'verticalPosition',
                        'axialPosition',
                        ]
                        
            units = ['K', 'K', 'K', '','','m/s', r'$\mu$J', 'J/kg', '', r'$\mu$$m^3$', 'mg m/s', 'mm', 'mm']
            A = zeros((len(arrs),)+ self.modelArray.shape)
           
            
            for i,M in ndenumerate(self.modelArray):
               
                    A[3][i]  = sum(M.PS.p.mesh.nodeVolume)/sum(M.PS.p.meshHistory[0].nodeVolume)
                    A[4][i]  = M.PS.p.diameter()/(2*M.PS.p.meshHistory[0].xInt[-1])
                    A[1][i]  = copy(M.PS.p.T[0])
                    A[2][i]  = sum(M.PS.p.mesh.nodeVolume*M.PS.p.T/sum(M.PS.p.mesh.nodeVolume))
                    A[0][i]  = copy(M.PS.p.surfaceTemperature)
                    A[8][i]  = M.PS.p.moltenFraction()
                    A[9][i]  = M.PS.p.moltenVolume()*1e18
                    A[6][i]  = sum(M.PS.p.enthalpy()*M.PS.p.mass*1e6)
                    A[7][i]  = M.PS.p.enthalpy()[-1]
                    A[5][i]  = sqrt(sum(M.PS.p.velocity()**2))
                    A[10][i]  = A[5][i]*M.PS.p.totalMass*1e6
                    A[11][i]  = M.PS.p.verticalPosition()*1e3
                    A[12][i]  = M.PS.p.axialPosition()*1e3
                    A[12][i]  = M.PS.p.axialPosition()*1e3
               
               
                    
            A[2][logical_or(isinf(A[2]), isnan(A[2]))] = 0
            
            if plotOn:
                F = [mlm.fig(a) for a in arrs]
                for i in range(len(arrs)):
                    try:
						mlm.surf3(g[0], g[1], A[i], axeLab = axlab+[arrs[i]+', '+units[i]], f = F[i])
                    except ValueError:
						print 'nans! check '+arrs[i]
    
            return dict(zip(['diameter', 'injectionVelocity']+arrs+['units'], 
                             [g[0], g[1]]+[A[i] for i in range(A.shape[0])] +[dict(zip(arrs, units))]  ))

class emittingParticleEnsemble(basicParticleEnsemble):
    def modelConstructor(self, parameterVector):
        """
        just use the model's constructor
        """
        return emittingParticleModel(parameterVector)


if __name__=='__main__':
    PE = basicParticleEnsemble(basicGriddedParticleParameterGenerator)
    print "constructed!"
    PE.generate(dMinMax = (15e-6,150e-6),vMinMax=(2,20), n  = (10,10))
    print 'generated!'
    print PE.modelArray.flatten()[0].PS.p.X
    PE.ensembleSolve(plasmaField = part.plasmaField(fastOn = True, temperatureScaling = 1.25))
    print 'solved!'
    print PE.modelArray.flatten()[0].PS.p.X
    S = PE.summary()
    
  
        
    #mlm.ml.show()    
    
    