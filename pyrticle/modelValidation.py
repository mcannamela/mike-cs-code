from numpy import *
seterr(invalid='raise')
import particle as part
#reload(part)
import particleSolver as PS

import scipy
import pdb
#import mlabMacros as mlm
#from matplotlib.pyplot import *
#try:
#    from mayavi import mlab as ml
#except ImportError:
#    from enthought.mayavi import mlab as ml
import cPickle as pkl
from plasmaField import *
import gasMixture as gm

def firstNRoots(fn, N = 10, startVal = 0):
    """
    find the first N roots of the one-d function fn,
    starting at startVal and increasing
    """
    findZero = lambda a_: scipy.optimize.brentq(fn, a_, Inf)
    a = startVal
    da = 1e-10
    x = zeros(N)
    for i in range(N):
        X = findZero(a)
        x[i]
        a = x[i]+da
        if r.converged == False:
            print 'warning, zero not found'
    return x


class rigidParticle(part.particle):
    """
    a particle that is rigidly fixed in space
    """
    def __init__(self, **kwargs):
        """
        enforce constant zirconia properties
        """
        kwargs['materialName'] = 'zirconiaConstant'
        part.particle.__init__(self, **kwargs)

    def fly(self):
        """
        probably not necessary, but
        """


        self.refreshCache()
        self.meshHistory+= [self.mesh]
        self.thermalTimestepHistory+=[self.thermalTimestep]
        self.X[self.stepCount] = self.X[0]
        self.stepCount+=1



    def heatTransferCoefficient(self):
        return 2*.88e4

class rigidTimeDependentJetParticle(part.timeDependentJetParticle):
    def __init__(self, **kwargs):
        """
        enforce constant zirconia properties
        """
        kwargs['materialName'] = 'zirconiaConstant'
        kwargs['material']= None
        part.particle.__init__(self, **kwargs)

    def fly(self):
        """
        probably not necessary, but
        """
        self.refreshCache()
        self.meshHistory+= [self.mesh]
        self.thermalTimestepHistory+=[self.thermalTimestep]
        if self.stepCount > 0:
                    self.t[self.stepCount] = (self.t[self.stepCount-1]
                                            +self.thermalTimestep)
        self.vmfHistory+=[self.vmf]
        self.plasmaTemperatureHistory+=[self.plasmaTemperature]
        self.X[self.stepCount] = self.X[0]
        self.stepCount+=1

    def heatTransferCoefficient(self):
        """
        don't correct for surfaceTemperature
        """
        return 2*.88e4


class rigidParticleSolver(PS.particleSolver):
        """
        use this class to solve rigidParticles for
        comparison with the solutions delivered by FourierParticleSolver
        """
        def __init__(self, **kwargs):
            #assert type(kwargs['particle']) ==type(rigidParticle()), 'this solver requires a rigidParticle object!'
            PS.particleSolver.__init__(self, **kwargs)

        def dimensionlessSoln(self):
            """
            return dimensionless temperature as a function of dimensionless radius
            """
            coord = self.p.mesh.xNode/self.p.mesh.xInt[-1]

            alpha = self.p.conductivity[-1]/(self.p.density[-1]*self.p.specificHeat()[-1])
            timeConstant = self.p.mesh.xInt[-1]**2/alpha
            t = cumsum(self.p.thermalTimestepHistory)
            
            Fo = t/timeConstant
            dimensionlessTemperature = ((self.getSoln()-self.p.plasmaTemperature)
                                                                    /(self.p.T0-self.p.plasmaTemperature))
            return (Fo, coord, dimensionlessTemperature)
class explicitRigidParticleSolver(PS.explicitParticleSolver):
        """
        use this class to solve rigidParticles for
        comparison with the solutions delivered by FourierParticleSolver
        """
        def __init__(self, **kwargs):
            #assert type(kwargs['particle']) ==type(rigidParticle()), 'this solver requires a rigidParticle object!'
            PS.explicitParticleSolver.__init__(self, **kwargs)

        def dimensionlessSoln(self):
            """
            return dimensionless temperature as a function of dimensionless radius
            """
            coord = self.p.mesh.xNode/self.p.mesh.xInt[-1]

            alpha = self.p.conductivity[-1]/(self.p.density[-1]*self.p.specificHeat()[-1])
            timeConstant = self.p.mesh.xInt[-1]**2/alpha
            t = cumsum(self.p.thermalTimestepHistory)
            
            Fo = t/timeConstant
            dimensionlessTemperature = ((self.getSoln()-self.p.plasmaTemperature)
                                                                    /(self.p.T0-self.p.plasmaTemperature))
            return (Fo, coord, dimensionlessTemperature)

class FourierParticleSolver(object):
    """
    solve heat transfer to the particle using fourier series
    """
    def __init__(self,theParticle, N = 100):
        """
        compute and cache the dimensionless quantities needed for the computation
        """
        self.N = N
        self.p = theParticle
        self.timestep = self.p.thermalTimestep
        self.nStep = self.p.nStep()

        #compute constant nondimensionalized
        self.nondimensionalize()
        self.t = linspace(self.timeConstant*1e-3,
                                    self.timeConstant*1e-3+self.timestep*self.nStep,
                                    self.nStep-1)

        self.computeEigenvalues()
        self.computeModeshapes()
        self.computeTransient()

        self.dimensionlessTemperature = sum(self.modeShapes[newaxis, ...]*self.transient, axis = 2)

    def nondimensionalize(self):
        """
        certain quantities are assumed constant and are therefore cacheable
        """
        self.coord = self.p.mesh.xNode/self.p.mesh.xInt[-1]
        self.Biot(self.p)
        self.alpha = self.p.conductivity[-1]/(self.p.density[-1]*self.p.specificHeat()[-1])
        self.timeConstant = self.p.mesh.xInt[-1]**2/self.alpha
        

    def Biot(self, p):
        """
        compute the biot number of a particle object
        """
        self.Bi = p.heatTransferCoefficient()*p.mesh.xInt[-1]/(p.conductivity[-1])
        print "Fourier Solver Biot nr: %.4f"%self.Bi



    def computeEigenvalues(self):
        """
        find the first N roots of the eigenvalue equation
        """
        f = lambda L: L*cos(L)+(self.Bi-1)*sin(L)

        delta = 1e-10
        period = pi

        a = 0
        self.eigenvalues = zeros(self.N)
        self.eigenvalues[0] = scipy.optimize.brentq(f,delta, period-delta)
        for i in range(1,self.N):
            self.eigenvalues[i] = scipy.optimize.brentq(f,a+i*period+delta, a+(i+1)*period-delta)

    def modalAmplitudes(self):
        """
        compute the modal amplitudes from the eigenvalues
        """
        L = self.eigenvalues
        return 2*(sin(L)-L*cos(L))*(L-sin(L)*cos(L))**-1


    def radialProfiles(self):
        """
        compute the radial profiles for each mode using the eigenvalues
        """
        return sinc((self.coord[...,newaxis]*self.eigenvalues[newaxis,...])/pi)


    def computeModeshapes(self):
        """
        modeshapes are the radialProfiles scaled by the modalAmplitudes
        """
        self.modeShapes = self.radialProfiles()*self.modalAmplitudes()[newaxis,...]

    def computeTransient(self):
        """
        transient term is exponential decay
        """
        theArg = -(self.eigenvalues[newaxis,newaxis,...]**2)*(self.t[...,newaxis,newaxis]/self.timeConstant)
        self.transient = exp(theArg)


    def dimensionalize(self, T0, Tinf):
        """
        return the dimensional temperature obtained if
        T0 = scalar, initial temperature
        Tinf = scalar, free stream temperature
        """
        r = self.coord*self.p.mesh.xInt[-1]
        T = Tinf+ self.dimensionlessTemperature*(T0-Tinf)

        return (self.t, r, T)

    def dimensionlessSoln(self):
        """
        return the dimensionless solution
        """
        Fo = self.t/self.timeConstant
        return (Fo, self.coord, self.dimensionlessTemperature)

    def plotAmplitude(self):
        """
        plot the log of the modal amplitudes
        """
        A = self.modalAmplitudes().copy()
        A[A==0]= np.min(A[A!=0])
#        mlm.surf3(x = self.coord, z = log10(A), f = mlm.fig('modalAmplitudes'))



if __name__== "__main__":
    pass
    close('all')
    timestep =1e-5
    nSteps = 1000
    nModes = 100
    nNodes = 20
    diameter = [.1e-6, 50e-6]
    T0 = 2000
    x0 = array([0, .01, 0])

    timeDependent = False

    g = gm.gasMixture('zirconiaConstant')
    

    if not timeDependent:
        p = rigidParticle(diameter = diameter,
                            T =T0,
                            thermalTimestep = timestep,
                            kinematicTimestep = timestep,
                            nNodes = nNodes,
                            nStep = nSteps,
                            material = g)
        p.adaptiveTimestepping=False
        p.splitScales = 0
        #create the plasmaField
        try:
            PF.density(array([0,0,0]))
        except NameError:
            PF = part.pf.plasmaField(fastOn = True)
    else:
        p = rigidTimeDependentJetParticle(diameter = diameter,
                            T =T0,
                            thermalTimestep = timestep,
                            kinematicTimestep = timestep,
                            nNodes = nNodes,
                            nStep = nSteps)
        p.adaptiveTimestepping=False
        p.splitScales = 0
        #create the plasmaField
        with open('/media/raidArray/CODE/jetModel/jet_sg100_ArPilot_40_500_DGMesh.pkl') as f:
            PF = pkl.load(f)

    p.inject(PF, x0,array([0, 0, 0]))

    S = []

#    print p.specificHeat()
    PS = explicitRigidParticleSolver(particle = p)
    PS.radiationOn = False
    PS.alpha = .5
    PS.stepErrorThres = 2e-5
    PS.maxIter = 50
#    print p.specificHeat()
    PS.solve()
#    print p.specificHeat()
   # f = PS.plotSoln()
    S+= [PS.dimensionlessSoln()]
    
    

    p.reset()
    FPS = FourierParticleSolver(p, nModes)
    S+= [ FPS.dimensionlessSoln()]
    lab = ['numerical', 'analytical']
    lStyle = ['g:', 'g-.', 'k--','k-']
#    print p.specificHeat()

    #close('all')
    figure(figsize = (22,20))
    #figure(figsize = (22,20))
    cnt = 0
    for s in S:
        figure(1)
        plot(s[0],(s[-1][:,0]),lStyle[2*cnt], linewidth = 1.5, marker = 'None',  label = lab[cnt]+ ', center')
        plot(s[0],(s[-1][:,-1]),lStyle[2*cnt+1], linewidth = 1.5, marker = 'None', label = lab[cnt]+ ', surface')
#        figure(2)
#        plot(s[0],log(s[-1][:,0])/s[0],linewidth = 1.5, marker = 'None',  label = lab[cnt]+ 'center')

        cnt+=1
    figure(1)
    xlabel(r'$\frac{t}{\tau}$',fontsize = 24)
    ylabel(r'$\frac{T-T_{inf}}{T_0-T_{inf}}$', fontsize = 24)
    title('Comparison of dimensionless solutions \n between numerical and analytical methods')
    legend()
    

    figure(figsize = (22,20))
    eSurf = (scipy.interp(S[1][0][:-10], S[0][0], S[0][-1][:,-1])-S[1][-1][:-10,-1])#*100/(S[1][-1][:-10,-1]+1e-4)
    eCenter = (scipy.interp(S[1][0][:-10], S[0][0], S[0][-1][:,0])-S[1][-1][:-10,0])#*100/(S[1][-1][:-10,0]+1e-4)
    plot(S[1][0][:-10],eCenter, 'k-', linewidth = 1.5, marker = 'None',  label = 'center')
    plot(S[1][0][:-10],eSurf, 'g--', linewidth = 1.5, marker = 'None', label = 'surface')

    xlabel(r'$\frac{t}{\tau}$')
    ylabel('% error')
    title('Error between numerical and analytical methods')
    legend()
#    figure(2)
#    xlabel('t/$tau$')
#    ylabel('exponent/time')
#    title('Comparison of exponential factors')
#    legend()

#    f = mlm.fig('error')
    e = (S[1][-1][2:]-S[0][-1][:-2])
#    mlm.surf3(x = S[0][0],y = S[0][1]*1e6, z = e, f = f)
#
#    f = mlm.fig('relative error')
#    mlm.surf3(x = S[0][0],y = S[0][1]*1e6 , z = e/S[1][-1][2:], f = f)

    v = p.mesh.nodeVolume[newaxis,:]
    v/=sum(v)*e.shape[0]

    print 'volume avg error x 1e6: %f'%(1e6*sum(v*e**2))
    print 'volume avg pct error: %f'%(100*sum(v*abs(e/S[1][-1][2:])))
