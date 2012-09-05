#here we implement a class that models a spherical particle in a fluid field
import helpers as hp
import interpolators as pol
import solvers as sol
import gasMixture as gm
import cPickle as pkl
import plasmaField as pf

import scipy
from numpy import *
from pylab import *
#import mlabMacros as mlm
import time
import pdb
import cProfile
import pstats
import os



class particle(object):
    """
    a particle immersed in hot fluid

    important methods:
    reset() - clear all histories and return the particle to its initial state
    inject() - initialize the particle into the fluid bath,
                    represented as a plasmaField object
    fly() - compute drag force and integrate to move the particle along in the field
    collapse() - remesh the particle if it is hollow and melts
    vaporize() - burn off portions of the particle if it gets too hot
    """

    def __init__(self,     diameter=float32([1.e-6, 50.e-6]),
                                x0=float32([0.,0.,0.]),
                                v0=float32([100.,0.,0.]),
                                T=float32(300.),
                                nNodes = float32(20.),
                                thermalTimestep = float32(2.e-5),
                                kinematicTimestep = float32(1.e-5),
                                nStep = uint16(500),
                                solidus_liquidus = float32([2912.,  3062.]),
                                vaporizationTemperature = float32(4600.),
                                heatOfVaporization = 5227000.,
                                heatOfFusion = 706300.,
                                vaporPressureFun = lambda T_: 101325.*exp(-4.5176e4/T_+10.088),
                                materialName = 'zirconiaGauss',
                                material = None,
                                plasmaField = None,
                                t0=0.0):
        """
        this class implements a hollow spherical particle in a hot fluid.

        x0 -  the location of the particle in the cartesian torch frame (x,y,z) where
                x = distance from torch exit plane
                y = distance above torch axis
                z = x cross y, the distance to the right of the torch looking down the torch from behind
        x0 - velocity of the particle in the same frame as X
        T - surface temperature of the particle
        rInterfaces - radial location of interfaces between nodes, first is inner shell radius last is outer shell radius
        rNodes - radial node locations node i is between interfaces i and i+1
        """
        with open('vaporModelDump', 'w') as f:
            pass
        self.nNodes = uint16(nNodes)
        self.vmf = 0.
        self.filmVaporFraction = 0.
        self.vaporPressure = 1.
        self.fv = 1.
        #need to know this to cope with the phase change
        self.solidus_liquidus = float32(solidus_liquidus) #depreciated, but may come back in the future
        self.vaporizationTemperature = float32(vaporizationTemperature)
        self.heatOfVaporization = heatOfVaporization
        self.heatOfFusion = heatOfFusion

        self.boltzmann = 1.3806503e-23
        self.vaporMassFluxThreshold = 0
        self.vaporPressureFun=vaporPressureFun

        #need to initialize these to prevent first-call errors
        self.Sh = 1;
        self.radialDecrement = 0
        self.radiusLost = 0
        self.nDecrement = 0


        #enforce two diameters
        if array(diameter).size==1:
            d = float32([.01e-6, diameter])
        else:
            d = float32((diameter))

        #radii of the interfaces and nodes
        self.mesh = sol.radial1DMesh(d[0], d[-1], nNodes)

        #timestepping stuff
        self.thermalTimestep = float32(thermalTimestep)
        self.kinematicTimestep = float32(kinematicTimestep)

        self.timestepRatio = int32(max([1, round(self.thermalTimestep/self.kinematicTimestep)])) #how many kinematic timesteps in one thermal timestep
        self.kinematicTimestep = self.thermalTimestep/self.timestepRatio #enforce integral timestep ratio

        self.splitScales = True and self.timestepRatio!=1  #split scales by default
        self.adaptiveTimestepping = True#adaptive timestepping by default
        self.thermalTimestepHistory = []


        self.maxTimestep = float32(10e-5)
        self.thermalStep   = float32(100)


        #position in the plume
        self.X = float32(zeros((nStep,3)))
        self.X[0] = float32(x0)

        #velocity in the plume
        self.V = float32(zeros((nStep,3)))
        self.V[0] = float32(v0)

        self.nStep = lambda:self.X.shape[0]

        self.t = float64(zeros(nStep))
        self.t[0] = t0

        #keep track of how many steps we've taken so we know which row of X to use
        self.stepCount = 1

        #initial temperature profile of the particle
        if array(T).size ==1:
            self.T = float32(T)*ones(nNodes,dtype = 'float32')
        else:
            self.T = float32(T)



        #initial temperature
        self.T0 = self.T.copy()
        self.surfaceTemperature = self.T[-1]


        #particle material
        self.materialName = materialName
        if material ==None:
            self.material = gm.gasMixture(self.materialName)
        else:
            self.material = material
        #setup callbacks for property fetching and particle state
        self.initializeNodePropertyCallbacks(g = self.material)
        self.initializeParticleStateCallbacks()

        #initialize specific enthalpy
        self.H = float32(self.material.Enthalpy(self.T))

        #cache particle properties that don't change (or change only with phase change)
        self.particleConstantCache()

        #initialization of other sundries
        self.initialMass = self.totalMass#so we know how much we started with if we boil off some mass
        self.terminate = False #throw true to inform the simulator to abort
        self.forceHistory = [] #keep track of the drag force over time
        self.meshHistory = [self.mesh] #keep track of the mesh over time
        self.vmfHistory = []
        self.plasmaTemperatureHistory = []

        #needed at init time for the timeDependentJetParticle class
        self.plasmaFraction = 1.0

        if plasmaField!=None:
            self.inject(plasmaField)




    def __getstate__(self):
        """
        need to exclude all lambda objects
        """
        black = ['diameter', 'verticalVelocity', 'velocity', 'axialVelocity', 'axialPosition', 'lateralPosition', 'liquidus',
                         'plasmaField', 'solidus', 'lateralVelocity', 'verticalPosition', 'enthalpy', 'specificHeat',
                         'temperatureLookup', 'position',  'time','nStep', 'vaporPressureFun']
        S = []
        K = []
        for k in self.__dict__.keys():
            if k in black:
                continue
            K+=[k]
            S+=[self.__dict__[k]]

        d = dict(zip(['enthalpy', 'specificHeat'],[self.enthalpy(), self.specificHeat() ]))

        return [dict(zip(K,S)), d]




    def __setstate__(self, D):
        """
        gonna have to call the methods that reconstruct the lambdas, which cannot be pickled
        """
        self.__dict__=D[0].copy()
        self.initializeParticleStateCallbacks()
        self.enthalpy = lambda: D[1]['enthalpy']
        self.specificHeat = lambda: D[1]['specificHeat']
        self.vaporPressureFun =  lambda T_: 101325.*exp(-4.5176e4/T_+10.088)

    def setEnthalpy(self, h):
        self.H = h
        self.T = self.material.Temperature(self.H)

            

    def reconstitute(self, PF, g):
        """
        after unpickling, we will need to re-set the plasma field and some callbacks
        """
        self.plasmaField = PF
        self.initializeNodePropertyCallbacks(g)

    def reset(self, x0 = None, v0 = None, T0 = None, t0 = None):
        """
        return the particle to x0, v0, and reset appropriate fields,
        prepping to run another sim
        """
        self.terminate = False

        if x0==None:
            x0 = self.X[0].copy()
        if v0==None:
            v0 = self.V[0].copy()
        if t0==None:
            t0 = self.t[0].copy()
        if T0!=None:
            self.T0 = T0

        self.X[...]=0
        self.V[...]=0
        self.t[...]=0

        #initialize position and velocity
        self.X[0] = x0
        self.V[0] = v0
        self.t[0] = t0
        self.T = self.T0

        #start at the beginning!
        self.stepCount = 1
        self.mesh = self.meshHistory[0]
        self.meshHistory = [self.meshHistory[0]]
        self.thermalTimestepHistory = []
        self.vmfHistory = []
        self.plasmaTemperatureHistory = []
        #update all quantities relating to the plasma
        self.particleConstantCache()
        self.refreshCache()
        self.vaporMassFlux()

    def inject(self, pField = pf.dummyField(), x0 = None, v0=None, t0=None):
        """
        shoot the particle into a plasmaField with initial position x0 and velocity v0
        """
        if x0==None:
            x0 = self.X[0].copy()
        if v0==None:
            v0 = self.V[0].copy()
        if t0==None:
            t0 = self.t[0].copy()

        self.plasmaField = pField

        self.Mplas = self.plasmaField.gas.meanMolarMass()
        #self.rplas = self.plasmaField.gas.meanMolecularRadius()
        self.Mp = self.material.meanMolarMass()
        self.rp = self.material.meanMolecularRadius()
        self.reset( x0, v0, None, t0)

    def fly(self):
        """
        integrate the equations of motion to move the particle through the plasma field for timestep seconds
        """
        if self.adaptiveTimestepping:
            self.getTimestep()

        r = self.diameter()/2
        a = (4.*pi/3)
        self.radialDecrement = r - (r**3-self.vmf*self.thermalTimestep/(a*self.density[-1]))**(1./3)




        if self.timestepRatio==1:
            self.splitScales = False

        if self.splitScales:
            T = zeros(self.timestepRatio)
            h = zeros(self.timestepRatio)
            for i in range(self.timestepRatio):
                if self.stepCount > 0:
                    self.t[self.stepCount] = (self.t[self.stepCount-1]
                                            +self.kinematicTimestep)
                self.vmfHistory+=[self.vmf]
                self.plasmaTemperatureHistory+=[self.plasmaTemperature]



                #get the acceleration
                a = self.dragForce()/self.totalMass
                self.forceHistory+=[hp.pnorm(a*self.totalMass)]


                #update position and velocity
                self.X[self.stepCount] = self.X[self.stepCount-1]+(.5*a*self.kinematicTimestep+self.velocity())*self.kinematicTimestep
                self.V[self.stepCount] = self.V[self.stepCount-1]+a*self.kinematicTimestep

                if self.V[self.stepCount][0]<0:
                    self.V[self.stepCount][0] = 0
                    print 'warning: reverse flight detected, solution likely unstable'
                if self.X[self.stepCount][0]<self.X[self.stepCount-1][0]:
                    self.X[self.stepCount][0] = self.X[self.stepCount-1][0]


                T[i] = self.plasmaTemperature
                h[i] = self.h

                #update all quantities relating to the plasma
                self.refreshCache()

                self.stepCount+=1

            #replace plasma temperature and heat transfer coefficient with their mean values
            #print [T, h]

            # self.plasmaTemperature = atleast_1d(mean(T))
            # self.h = atleast_1d(mean(h))

            Th = atleast_1d(mean(T*h))
            self.plasmaTemperature = Th/mean(h)
            self.h = mean(h)
        else:
            #get the acceleration
                if self.stepCount > 0:
                    self.t[self.stepCount] = (self.t[self.stepCount-1]
                                            +self.thermalTimestep)

                self.vmfHistory+=[self.vmf]
                self.plasmaTemperatureHistory+=[self.plasmaTemperature]

                a = self.dragForce()/self.totalMass
                self.forceHistory+=[hp.pnorm(a*self.totalMass)]

                #update position and velocity
                self.X[self.stepCount] = self.X[self.stepCount-1]+(.5*a*self.thermalTimestep+self.velocity())*self.thermalTimestep
                self.V[self.stepCount] = self.V[self.stepCount-1]+a*self.thermalTimestep

                self.stepCount+=1

                #update all quantities relating to the plasma
                self.refreshCache()
        #increment time, update mesh history
        self.meshHistory+= [self.mesh]
        self.thermalTimestepHistory+=[self.thermalTimestep]



#        with open('vaporModelDump', 'a') as f:
#            f.write('%.2e  %.2e  %.2e  %.2e  %.2e  %.2e  %.2e  %.2e  %.2e\n'%(self.surfaceTemperature,
#                                                  self.plasmaTemperature, self.Sh, self.vaporPressure,
#                                                  self.filmVaporFraction, self.B, self.diffusionCoeff,
#                                                  self.vmf*1e12, self.fVap()))

    def collapse(self):
        """
        if the particle is hollow, collapse it into a solid particle and remesh
        """
        if (self.mesh.xInt[0]/self.mesh.xInt[-1])>.05:
            #get diameter from the total volume of the particle
            d= 2*(3*sum(self.mesh.nodeVolume)/(4*pi))**(1./3.)

            #temperature of the collapsed particle comes from the mean specific enthalpy
            H = sum(self.mass*self.enthalpy())/self.totalMass
            T = self.temperatureLookup(H)

            print str([sum(self.T*self.mass)/self.totalMass , T  ])

            self.T = ones_like(self.T)*T
            self.mesh = sol.radial1DMesh(.02*d, d, self.mesh.nNodes)

            #recompute mesh dependent quantities
            self.particleConstantCache()

            self.collapsedStep = self.stepCount
            print 'collapsed'

    def vaporize(self):
        """
        if outer shells have reached the vapor point they must be discarded
        """

        r = self.diameter()/2
        a = (4.*pi/3)
        self.radialDecrement = r - (r**3-self.vmf*self.thermalTimestep/(a*self.density[-1]))**(1./3)
        self.radiusLost+=self.radialDecrement
        self.nDecrement+=1

        #print '%.2e, %.3e, %d, %.3e'%(self.radialDecrement, self.radiusLost, self.nDecrement, self.mesh.xInt[-1])


        #print 'reducing radius by %.2e'%self.radialDecrement
        R = self.mesh.xInt[-1]
        Rn = self.mesh.xNode[-1]
        dR = self.radialDecrement

        if (R-dR)>(Rn):
            self.mesh.xInt[-1]-=dR
            self.mesh.interfaceArea[-1] = self.mesh.area()[-1]
            self.mesh.nodeVolume[-1] = self.mesh.volume()[-1]
            self.mesh.totalVolume = sum(self.mesh.nodeVolume)
           # print '%.3e'%self.mesh.xInt[-1]
          #  print "vapor flux is %.2e cubic mics/microsecond, reducing mesh by %.2e mics"%(self.vmf*1e12, self.radialDecrement*1e6)

        else:
            #print "(xInt[-1]-dx, xInt[-2]) = (%.2e, %.2e)"%(self.mesh.xInt[-1]-self.radialDecrement, self.mesh.xInt[-2])
#            pdb.set_trace()
            oldMesh = self.mesh
            self.mesh = sol.radial1DMesh(2*self.mesh.xInt[0], 2*(self.mesh.xInt[-1]-self.radialDecrement), self.mesh.nNodes)

            #interpolate to get T
            self.T = scipy.interpolate.interp1d(r_[oldMesh.xInt[0], oldMesh.xNode], r_[self.T[0], self.T] )(self.mesh.xNode)
            self.H = float32(self.material.Enthalpy(self.T))

            #recompute mesh dependent quantities
        self.particleConstantCache()

        #determine whether the particle is totally gone
        if self.totalMass<=.10*self.initialMass:
            self.terminate = True
            print 'disintegrated'


    def getTimestep(self):
        """
        compute the timestep from the gradient of the plasma field
        """
        gradT = self.plasmaField.temperatureGradient(self.position())
        v = array([sqrt((self.lateralVelocity()**2+self.verticalVelocity()**2)), self.axialVelocity()])

        # gradT dot velocity gives the time derivative of forcing temperature
        denom = abs(sum(gradT*v))

        #enforce a maximum timestep
        oldTS = 0
        oldTS += self.thermalTimestep
        self.thermalTimestep = min([self.thermalStep/denom, self.maxTimestep, 1.1*oldTS])

        #kinematicTimestep in a fixed ratio for now
        self.kinematicTimestep = self.thermalTimestep/self.timestepRatio

        #update thermalMassConstant
        self.thermalMassConstant = self.mass/self.thermalTimestep


    def dragForce(self):
        """
        compute the drag on the particle
        """
        #coefficient of drag is a fn of Reynold's nr
        #A.F. Mills, "Heat Transfer", p320
        #Cd = 24/self.Re*(1+(self.Re**(2./3.))/6)
        
        #chiang, numerical analysis of convecting, vaporizing fuel droplet with variable properties
        Cd = (24/self.Re)*(1+.325*self.Re**.474)
        
        #print 'Re, Cd, fProp, fVap, Tplasma-Tpart: %.1f,  %.2f, %.2f, %.2f, %d'%(self.Re, Cd, self.fProp(), self.fVap(), self.plasmaTemperature-self.surfaceTemperature)
        

        #Ono used this formula
#        Cd_old=24/self.Re+6/(1+pow(self.Re,0.5)+.4)
#        print str([Cd- Cd_old])

        Cd*=self.fProp()

        #speed relative to the plasma determines the magnitude of the drag force
        speed = hp.pnorm(self.relativeVelocity)

        #plasma density, particle area, and square of speed premultiply the drag coeff
        #prefactor = (.5*self.plasmaDensity*speed**2*pi*(self.diameter()/2.)**2)
        prefactor = (.5*self.filmDensity()*speed**2*pi*(self.diameter()/2.)**2)
        

        #drag force is directed along the relative velocity vector
        Fd = prefactor*Cd*(self.relativeVelocity/max([speed, .001]))

        #debug statements
        #print str(self.plasmaVelocity[0])+'  '+ str(self.velocity()[0])+'   '+str(Fd[0])

        return Fd

    def heatTransferCoefficient(self):
        """
        compute the heat transfer coefficient between the plasma and particle
        """
        #wan, sampath, fincke, mass transfer factor
        
        h = float32(self.h*self.fVap())
        h*= self.fProp()

        #print h
        return h

    def fProp(self):
        """
        variable properties factor for heat transfer coefficient
        """
        #A.F. Mills, "Heat Transfer", p320
        fmu =  pow( self.plasmaField.gas.Viscosity(self.surfaceTemperature)/self.plasmaViscosity, -.25 )
        fcp = (self.plasmaSpecificHeat/self.plasmaField.gas.SpecificHeat(self.surfaceTemperature))**.38
        print 'fProp %.2f, fcp  %.2f, fmu  %.2f'%(fcp*fmu,fcp, fmu)
        return fcp*fmu
#        return 1

    def fVap(self):
        """
        vaporization factor for heat transfer coeff
        """
        if self.plasmaTemperature>=self.surfaceTemperature:               
            xi=  self.vmf*self.fVapConst()
            try:
                fv =  xi/(exp(xi)-1)
            except FloatingPointError:
                fv = 1.0
            if isnan(fv) or isinf(fv):
                fv = 1.
            return max([min([fv,1.0]),0])
        else:
            return 1.0

    def fVapConst(self):
        """
        constant applied to vapor mass flux to generate argument for fVap formula
        """
        return self.filmSpecificHeat()/(pi*self.diameter()*self.filmConductivity())
        #return self.plasmaSpecificHeat/(pi*self.diameter()*self.plasmaConductivity)
    def B_h_0(self):
        return self.filmSpecificHeat()*(self.plasmaTemperature-self.surfaceTemperature)/self.heatOfFusion

    def vaporMassFlux(self):

        #for default vapor pressure relation, see: simon. mathematical modeling
        #of the melt pool during a physical vapor deposition process.
        #MS thesis, MIT 1998
        self.vaporPressure = self.vaporPressureFun(self.surfaceTemperature)
        if self.vaporPressure<= .001*self.plasmaPressure:
            self.vmf = 0.
            self.B = 0.
            self.filmVaporFraction = 0.
            self.diffusionCoeff = 0.
            return self.vmf

#        if self.vaporPressure >=.99*self.plasmaPressure and self.plasmaTemperature>.99*self.vaporizationTemperature:

        #boiling vmf

        if self.plasmaTemperature>self.vaporizationTemperature:
            a = self.fVapConst()
            A = self.mesh.interfaceArea[-1]
            dT = max([self.plasmaTemperature-self.vaporizationTemperature, 0])
            transferVMF = (a**-1)
            transferVMF *= log(1+ a*self.h*self.fProp()*A*dT/self.heatOfVaporization)
        else:
            transferVMF = inf



        self.filmVaporFraction = self.Mp/(self.Mp+
             self.Mplas*max([self.plasmaPressure/self.vaporPressure-1,0.]))

        self.B = self.filmVaporFraction/max([1-self.filmVaporFraction, 1e-4])


        self.diffusionCoeff = 1.6e-1*((2.68e-3*self.plasmaPressure/1e24)*(self.filmTemperature()**3/
                                                                        self.material.meanMolarMass())**.5*
                                                                        (self.plasmaPressure*self.material.meanMolecularRadius()**2)**-1)


        nonBoilingVMF = self.filmDensity()*self.diffusionCoeff*pi*self.diameter()*log(1+self.B)*self.Sh
        self.vmf = max([0,nonBoilingVMF])
        
#        if not isinf(transferVMF):
#            if self.T[-1]>(self.vaporizationTemperature-200):
#                self.vmf = (1./transferVMF+1./nonBoilingVMF)**-1
#            elif self.T[-1]>(self.vaporizationTemperature-200):
#                self.vmf = transferVMF
            
            
        
        
#        try:
#            self.vmf= .5*self.vmfHistory[-1]+.5*min([transferVMF, nonBoilingVMF])
#        except IndexError:
#            self.vmf= min([transferVMF, nonBoilingVMF])
#
#        if nonBoilingVMF > transferVMF:
#            self.T[-1] = self.vaporizationTemperature
#            self.surfaceTemperature = self.T[-1]
#            self.vaporMode = 1
#        else:
#            self.vaporMode = 0


        #self.vmf = self.plasmaDensity*self.diffusionCoeff*pi*self.diameter()*log(1+self.B)*self.Sh

        #print (self.filmDensity(), self.diffusionCoeff, self.filmDensity()*self.diffusionCoeff)
#        if self.vmf>2e-7:
#            pdb.set_trace()

        return self.vmf


    def filmTemperature(self):
        """
        use the current temperature of the particle and surrounding fluid to get the film temperature for the
        convective BC
        """
        return ( self.surfaceTemperature+self.plasmaTemperature )/2

    def filmDensity(self):
        return self.plasmaField.gas.Density(self.filmTemperature())

    def filmConductivity(self):
        return self.plasmaField.gas.Conductivity(self.filmTemperature())

    def filmSpecificHeat(self):
        return self.plasmaField.gas.SpecificHeat(self.filmTemperature())

    def filmViscosity(self):
        return self.plasmaField.gas.Viscosity(self.filmTemperature())


    def nodalThermalMass(self):
        """
        compute the thermal mass of the nodes using the current temperature profile
        """
        return(self.thermalMassConstant*self.specificHeat())

    def initializeNodePropertyCallbacks(self, g = None):
        """
        use the gas mixture class to set up property fetching for the particle nodes
        """
        #it's not a gas but the object doesn't know that!
        if g==None:
            g = gm.gasMixture(self.materialName)

        #dictionary of property fetchers with keys:
        #['Temperature', 'Conductivity', 'Density', 'Enthalpy', 'SpecificHeat', 'Viscosity']
        cb= lambda propName_: lambda: g.getCallbacks()[propName_](self.T)

        self.conductivity = ones_like(self.T)*g.Conductivity(300)
        self.enthalpy = cb('Enthalpy')
        self.density = ones_like(self.T)*g.Density(300)
        self.specificHeat = cb('SpecificHeat')
        self.temperatureLookup = g.getCallbacks()['Temperature']

    def initializeParticleStateCallbacks(self):
        """
        set up small functions to return basic states of the particle
        """


        #leave a backdoor for changes
        self.diameter = lambda: 2*self.mesh.xInt[-1]
        self.solidus = lambda: self.solidus_liquidus[0]
        self.liquidus = lambda: self.solidus_liquidus[1]
        self.position = lambda: self.X[self.stepCount-1]
        self.time =     lambda: self.t[self.stepCount-1]

        self.velocity = lambda: self.V[self.stepCount-1]

        #position and velocity
        v = lambda i_: lambda: self.V[self.stepCount-1].copy()[i_]
        x = lambda i_: lambda: self.X[self.stepCount-1].copy()[i_]
        self.axialVelocity    =     v(0)
        self.verticalVelocity =     v(1)
        self.lateralVelocity  =     v(2)
        self.axialPosition    =     x(0)
        self.verticalPosition =     x(1)
        self.lateralPosition  =     x(2)


    def particleConstantCache(self):
        """
        cache some constant particle properties
        """
        #cacheable particle states
        self.mass = self.mesh.nodeVolume*self.density
        self.thermalMassConstant = self.mass/self.thermalTimestep
        self.totalMass = sum(self.mass)

    def refreshCache(self):
        """
        set up cachable quantities
        """

        #properties of the plasma at this location
        self.plasmaTemperature = self.plasmaField.temperature(self.position())
        self.plasmaVelocity =  array([self.plasmaField.velocity(self.position()), 0, 0])
        self.plasmaDensity = self.plasmaField.density(self.position())
        self.plasmaViscosity = self.plasmaField.viscosity(self.position())
        self.plasmaConductivity = self.plasmaField.conductivity(self.position())
        self.plasmaSpecificHeat = self.plasmaField.specificHeat(self.position())
        self.plasmaPressure = 100000.

        #deal with any phase changes
        #if any(self.T>self.vaporizationTemperature):
#        try:
        #self.vaporMassFlux()
        if self.vmf > self.vaporMassFluxThreshold:
            self.vaporize()
        if all(self.T>self.solidus_liquidus[-1]) and not self.terminate:
            self.collapse()
#        except AttributeError:
#            print "couldn't compute vapor mass flux, probably a first call error you can ignore"

        #cacheable particle properties
        self.surfaceTemperature = self.T[-1]
        self.relativeVelocity =self.plasmaVelocity-self.velocity()

        #dimensionless quantities
        self.Re = max([self.filmDensity()
                    *hp.pnorm(self.relativeVelocity)
                    *self.diameter()/self.filmViscosity(),.001])
        #print '%d,%.2e, %.2e, %.2e, %.2e, %.2e'%(self.stepCount,self.axialPosition(), self.plasmaDensity, sign(self.plasmaVelocity[0]-self.axialVelocity())*hp.pnorm(self.relativeVelocity),self.plasmaViscosity , self.Re)

        #Pr from film temperature is very large, like 2 or 3
        #self.Pr = self.plasmaViscosity*self.plasmaSpecificHeat/self.plasmaConductivity
        self.Pr = self.filmViscosity()*self.filmSpecificHeat()/self.filmConductivity()

        #A.F. Mills, "Heat Transfer", p320
        #print (self.Pr, self.Prf)
        self.Nu = float32(2 + ( .4*pow(self.Re,.5)+.06*pow(self.Re,2./3.) )*pow(self.Pr,.4))

        #wan and fincke, a bit higher than Mills Nu
        #self.Nu = float32(2 +  .6*self.Re**.5*self.Pr**.33333)
        #print (self.Nu, self.Nuw)

        self.h = (self.filmConductivity()*self.Nu/self.diameter())

        #wan, sampath, fincke. Model and powder particle heating,
        #melting, resolidification and evaporation in plasma spraying process.
        #j. heat transfer 1999.
        self.Sc = .7
        self.Sh = 2+.6*self.Re**.5*self.Sc**(1./3)
                ##alias for readability


        #print str(self.position()[0])+ '  '+ str(self.plasmaTemperature)+ '  '+ str(self.plasmaVelocity[0])



    def simpleStateVector(self):
        """
        a summary of the state of the particle
        """
        return (self.diameter(), self.position(), self.velocity(), self.surfaceTemperature)

    def simpleStateVectorString(self):
        """
        format the state summary into a string
        """
        return str(self.simpleStateVector).replace('[','').replace(']','').replace('(','').replace(')','')
    def moltenFraction(self, T = None):
        """
        compute the molten fraction of the particle
        """
        if T==None:
            T = self.T

        mf = (T-self.solidus())/(self.liquidus()-self.solidus())
        mf[mf<0] = 0
        mf[mf>1] = 1
        MF = sum(mf*(self.mesh.nodeVolume/sum(self.mesh.nodeVolume)))
        if isnan(MF) or isinf(MF):
            MF = 0.
        return MF

    def moltenVolume(self, T = None):
        """
        compute the molten fraction of the particle
        """
        if T==None:
            T = self.T

        mf = (T-self.solidus())/(self.liquidus()-self.solidus())
        mf[mf<0] = 0
        mf[mf>1] = 1
        return sum(mf*(self.mesh.nodeVolume))

    def norm(self):
        """
        the magnitude of the particle's state
        """
        return self.totalMass





class emittingParticle(particle):
    """
    like the particle, but with the ability to compute it's own emission spectrum
    the default emitting particle is an opaque blackbody
    """
    def intensityProfile(self, nx, wavelength, T = None):
        """
        for each wavelength in the array wavelength,
        produce a function centered about x = 0 which gives the intensity profile of the particle at that wavelength
        and sample the function at nx evenly spaced points

        nx - number of samples in the output image
        wavelength - array of wavelengths
        T - array of temperatures
        """
        if T == None:
            T = self.T
        if any(wavelength>1e1):
            print 'warning: long wavelength detected, what do you think you have, a radio telescope!? \
                    \nwavelength argument  must be in meters.'
            print 'assuming wavelengths in nanometers and converting...'
            Wavelength = float64(wavelength)*1e-9
        else:
            Wavelength = float64(wavelength)



        D = self.extinctionDepth(wavelength, T)
        E = self.spectralEmissivity(wavelength, T)
        R = self.spectralRadiance(wavelength, T)



        #shape in the 0th axis, wavelength in the 1st axis
        #!!! need to multiply by the wavelength bin width!!!
        iP = self.shapeFn(nx)[...,newaxis]*sum(E*R*self.mesh.interfaceArea[1:,newaxis]*
                                exp( -(self.mesh.xInt[-1]-self.mesh.xInt[1:,newaxis])/D) , axis = 0)


        return iP

    def extinctionDepth(self, wavelength, T= None):
        """
        return the extinction depth at given wavelengths and temperatures
        """
        #optically thick - surface radiation only
        if T == None:
            T = self.T
        return zeros((self.nNodes,)+wavelength.shape)+1e-10

    def spectralEmissivity(self, wavelength, T = None):
        """
        return the spectral emissivity at given wavelengths and temperatures
        """
        #black body assumption
        if T == None:
            T = self.T
        return ones((self.nNodes,)+wavelength.shape)

    def spectralRadiance(self,wavelength, T=None):
        """
        spherical intensity of radiated light at each wavelength and temperature, in W/steradian
        if T = None, the particle's current temperature profile is used

        wavelength - array of wavelengths
        T - array of temperatures
        """
        if T == None:
            T = self.T

        #physical constants, planck, speed o' light, boltzmann
        h=6.626e-34; #J*s
        c=3e8; #m/s
        kB= 1.381e-23;#J/K
        const=8*pi*h*c;

        #T varies along axis 0, wavelength along axis 1
        radiance=const*pow(expand_dims(wavelength,0),-5)*pow(exp(h*c/(outer(T,wavelength)*kB )) -1,-1) #W/m^3/steradians
        return radiance

    def shapeFn(self, nx, fac = 10.):
        """
        compute the shape of the intensity profile from geometric considerations and
        the paticle's velocity
        """
        #oversample by a factor of fac, smooth with a window of fac, downsample to nx

        R = float64(nx)*fac/2
        r = arange(float64(ceil(nx)*fac)) - float64(R)

        y = sqrt(R**2- r**2).reshape(ceil(nx), fac)
        y[isnan(y)] = 0

        Y = float64(squeeze(sum(y, axis = 1)))


        Y/= sum(Y)*self.axialVelocity()#need to set axialVelocity nonzero!

        if any(isnan(Y)) or any(isinf(Y)):
            pdb.set_trace()

        if nx==1:
            return atleast_1d(Y)
        else:
            return Y

class timeDependentJetParticle(particle):
    def inject(self, pField = pf.dummyField(), x0 = None, v0=None, t0=None):
        """
        shoot the particle into a plasmaField with initial position x0 and velocity v0
        """
        if x0==None:
            x0 = self.X[0].copy()
        if v0==None:
            v0 = self.V[0].copy()
        if t0==None:
            t0 = self.t[0].copy()

        self.plasmaField = pField

        self.Mplas = self.plasmaField.gas.meanMolarMass(self.plasmaFraction)
        #self.rplas = self.plasmaField.gas.meanMolecularRadius()
        self.Mp = self.material.meanMolarMass()
        self.rp = self.material.meanMolecularRadius()
        self.reset( x0, v0, None, t0)
    def getTimestep(self):
        """
        compute the timestep from the gradient of the plasma field
        """
        gradT = self.plasmaField.temperatureGradient(self.position(), self.time())
        v = array([sqrt((self.lateralVelocity()**2+
            self.verticalVelocity()**2)), self.axialVelocity()])

        # gradT dot velocity gives the time derivative of forcing temperature
        denom = abs(sum(gradT*v))

        #enforce a maximum timestep
        oldTS = 0
        oldTS += self.thermalTimestep
        self.thermalTimestep = min([self.thermalStep/denom, self.maxTimestep, 1.25*oldTS])

        #kinematicTimestep in a fixed ratio for now
        self.kinematicTimestep = self.thermalTimestep/self.timestepRatio

        #update thermalMassConstant
        self.thermalMassConstant = self.mass/self.thermalTimestep

    def fProp(self):
        """
        variable properties factor for heat transfer coefficient
        """
        #A.F. Mills, "Heat Transfer", p320
        fmu =  pow( self.plasmaField.gas.Viscosity(self.surfaceTemperature, self.plasmaFraction)/self.plasmaViscosity, -.25 )
        fcp = (self.plasmaSpecificHeat/self.plasmaField.gas.SpecificHeat(self.surfaceTemperature, self.plasmaFraction))**.38
        
        #print 'fProp %.2f, fcp  %.2f, fmu  %.2f'%(fcp*fmu,fcp, fmu)
        return fcp*fmu
#        return 1
    def refreshCache(self):
        """
        set up cachable quantities
        """

        #properties of the plasma at this location
        T, V, f, rho, mu, k, Cp, M = self.plasmaField(self.position(), self.time())
        self.plasmaFraction = f
        self.plasmaTemperature = T#self.plasmaField.temperature(self.position(), self.time())
        self.plasmaVelocity =  V#self.plasmaField.velocity(self.position(), self.time())
        self.plasmaDensity = rho#self.plasmaField.density(self.position(), self.time())
        self.plasmaViscosity = mu#self.plasmaField.viscosity(self.position(), self.time())
        self.plasmaConductivity = k#self.plasmaField.conductivity(self.position(), self.time())
        self.plasmaSpecificHeat = Cp#self.plasmaField.specificHeat(self.position(), self.time())
        self.Mplas = M
        self.plasmaPressure = 100000.

        #deal with any phase changes
        #if any(self.T>self.vaporizationTemperature):
        #self.vaporMassFlux()
        if self.vmf > self.vaporMassFluxThreshold:
            self.vaporize()
        if all(self.T>self.solidus_liquidus[-1]) and not self.terminate:
            self.collapse()

        #cacheable particle properties
        self.surfaceTemperature = self.T[-1]+(self.mesh.xInt[-1]-self.mesh.xNode[-1])*diff(self.T[-2:])/diff(self.mesh.xNode[-2:])
        self.relativeVelocity =self.plasmaVelocity-self.velocity()

        #dimensionless quantities
        self.Re = max([self.filmDensity()
                    *hp.pnorm(self.relativeVelocity)
                    *self.diameter()/self.filmViscosity(),.001])
        #print '%d,%.2e, %.2e, %.2e, %.2e, %.2e'%(self.stepCount,self.axialPosition(), self.plasmaDensity, sign(self.plasmaVelocity[0]-self.axialVelocity())*hp.pnorm(self.relativeVelocity),self.plasmaViscosity , self.Re)

        
        #self.Pr = self.plasmaViscosity*self.plasmaSpecificHeat/self.plasmaConductivity
        self.Pr = self.filmViscosity()*self.filmSpecificHeat()/self.filmConductivity()
        

        #A.F. Mills, "Heat Transfer", p320
        #print (self.Pr, self.Prf)
        #self.Nu = float32(2 + ( .4*pow(self.Re,.5)+.06*pow(self.Re,2./3.) )*pow(self.Pr,.4))
        
        #from Sazhin
        self.Nu = 2+(.4*self.Re**.5-.06*self.Re**.6666)*self.Pr**.3333 

        #wan and fincke, a bit higher than Mills Nu
        #self.Nu = float32(2 +  .6*self.Re**.5*self.Pr**.33333)
        #print (self.Nu, self.Nuw)

        self.h = (self.filmConductivity()*self.Nu/self.diameter())
        #print (k, self.filmConductivity(), self.filmTemperature(), T, f)

        #wan, sampath, fincke. Model and powder particle heating,
        #melting, resolidification and evaporation in plasma spraying process.
        #j. heat transfer 1999.
        self.Sc = .7
        self.Sh = 2+.6*self.Re**.5*self.Sc**(1./3)
                ##alias for readability


        #print str(self.position()[0])+ '  '+ str(self.plasmaTemperature)+ '  '+ str(self.plasmaVelocity[0])

    def filmDensity(self):
        return self.plasmaField.gas.Density(self.filmTemperature(), self.plasmaFraction)

    def filmConductivity(self):
        return self.plasmaField.gas.Conductivity(self.filmTemperature(), self.plasmaFraction)

    def filmSpecificHeat(self):
        return self.plasmaField.gas.SpecificHeat(self.filmTemperature(), self.plasmaFraction)

    def filmViscosity(self):
        return self.plasmaField.gas.Viscosity(self.filmTemperature(), self.plasmaFraction)


if __name__== "__main__":
    close("all")

    radianceTestOn = False
    timeDependentTestOn = False
#    fieldTestOn = False
#    lavaFieldTestOn = True
#
#    if lavaFieldTestOn :
#        pf = lavaField(simName = "lavaJet.pkl", gasPickle = None)
#        #pf = lavaField(gasPickle = 'lavaJet_gas.pkl')
#        mlm.surf3(z = pf.Density, f = mlm.fig('density') )
#
#
#
#
#    if fieldTestOn:
#        pf = plasmaField('..\\pyrticle\\cur-700_flo-50_gas-argon68helium32_make-SG100')
#        dpf = dummyField()
#
#        print 'plasma density'
#        print pf.density(zeros(3))
#
#        print 'dummy plasma density'
#        print dpf.density(zeros(3))

    if radianceTestOn:
        P = emittingParticle()
        T = arange(3000,5501,500)
        Lam = arange(300, 2000,20)
        planckRadiance = P.spectralRadiance(Lam*1e-9,T)
        figure()


        for i in range(planckRadiance.shape[0]):
            plot(Lam, planckRadiance[i], label = str(T[i]))


        title('planck distribution for some temperatures')
        xlabel('wavelength')
        ylabel('radiance')
        legend()


        P.T*=0
        P.T+=4000
        figure()
        plot(Lam,P.spectralRadiance(Lam*1e-9, T = P.T[-1] ).T)
        title('spectral radiance of particle P at T = %d K'%P.T[-1] )
        xlabel('wavelength')
        ylabel('radiance')

        ip = P.intensityProfile(10, Lam)
        show()

    if timeDependentTestOn:
        runFolder = "/media/bigBoy/openFoamRuns/jet_sg100_ArPilot_40_500_DGMesh"

        with open('/media/raidArray/CODE/gasModel/argonPropertiesDict.pkl', 'rb') as f:
            DAr = pkl.load(f)
            DAr['enthalpy']-=DAr['enthalpy'][0]
        Ar = pf.gas(40.0, DAr, 10000)


        with open('/media/raidArray/CODE/gasModel/airPropertiesDict.pkl', 'rb') as f:
            DA = pkl.load(f)
            DA['enthalpy']-=DA['enthalpy'][0]
        Air = pf.gas(29.0, DA, 10000)

        G = pf.binaryGasMixture(Ar, Air)

        tMin = .00098
        tMax = .00104
        xMin = array([0.0, -.01, -.01])
        xMax = array([.03, .01, .01])

        P = pf.openFoamPlasma(runFolder, tMin, tMax, xMin,xMax, G)

        p = timeDependentJetParticle(plasmaField = P)

