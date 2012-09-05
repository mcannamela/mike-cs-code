import helpers
from helpers import *
import solvers as sol
import particle as part
import scipy
from pylab import *
from plotMacros import *
from numpy import *
import cPickle
import time
import pdb
import cProfile
import pstats
import threading
#import mlabMacros as mlm
#try:
#    from mayavi import mlab as ml
#except ImportError:
#    from enthought.mayavi import mlab as ml


class particleSolver(sol.marchingSolver):
    """
    solves the equations of heat a particle in a hot fluid
    """
    def __init__(self, particle = part.particle(thermalTimestep = 1e-6),
                        xMax = .1,
                        maxSteps = 1000,
                        verbose = False,
                        sema = None):
        """
        initialize using the marching solver's init method and the particle's number of nodes
        """
        #whether to print extra messages
        self.verbose = verbose

        #set the particle datamember
        self.p = particle

        #stop the sim if we pass xMax
        self.xMax = xMax

        #stop the sim if we

        self.maxSteps = fix(double(self.p.X.shape[0])/self.p.timestepRatio)
        sol.marchingSolver.__init__(self, self.p.T, self.p.mesh.nNodes, nVars = 1, verbose = verbose, sema = sema)#last node doesn't count since it's the boundary temperature


        self.convHistory = []
        self.radHistory = []
        self.vapHistory = []


        #debugging messages
        if self.verbose:
            print 'injection at position '+ str(self.p.position())
            print 'with velocity '+ str(self.p.velocity())
            #print 'at temperature '+ str(self.currentSoln())
            print 'ambient temperature is ' + str(self.p.plasmaTemperature)
            print 'with velocity ' + str(self.p.plasmaVelocity)
            
        self.radiationOn = True

    def __getstate__(self):
        """
        need to exclude all lambda objects
        """
        black =['']
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


    def reset(self, x0 = None, v0 = None, t0 =None):
        """
        return the particle to its initial condition, prepping to run another sim
        """
        self.p.reset(x0,v0, t0 = t0)
        sol.marchingSolver.reset(self)

    def lastSoln(self):
        return self.solvedSteps[-1][0]
    def currentSoln(self):
        return self.soln[0].copy()


    def boundaryTemperature(self):
        return self.p.plasmaTemperature
#        return self.p.filmTemperature()

    def setCurrentSolution(self,soln):
        """
        insert a guess into the current solution
        """
        self.soln[0] = soln

    def solveStep(self):
        """
        like marchingSolver.solveStep, but
        update the particle's temperature mesh at every step
        """
#        sol.marchingSolver.solveStep(self)

        i = 0
        while self.stepError()>self.stepErrorThres and i<self.maxIter or i==0:
            self.updateCoeffs()
            self.prevIterSoln = self.soln.copy()
            for j in arange(self.nVars):
                deltaT = (self.tdmaSolve(self.coeffs[j], self.rhs[j])-self.soln[j])
                self.soln[j] += self.alpha*deltaT
                #print str(mean(deltaT))+'   ' +str(np.max(deltaT))
            self.p.T = self.currentSoln()

            i+=1
        if i>=self.maxIter:
            print 'iteration limit reached'

        self.itercount+=[i]
        



    def stoppingCriterion(self):
        """
        stop the sim if the particle flies past xMax or takes longer than tMax
        """
        maxPositionReached = (self.p.axialPosition()>=self.xMax)
        maxStepsReached = (self.currentStep>=(self.maxSteps-2*self.p.timestepRatio))
        if maxStepsReached:
            print 'max steps exceeded'

        return  maxPositionReached or  maxStepsReached or self.p.terminate

    def updateCoeffs(self):
        """
        use the particle's properties, current temperature, and mesh to get the coeffs
        """

        #-----------------get ready---------------------

        #sundries for the rhs and boundary condition
        tm = self.p.nodalThermalMass()
        vmf = self.p.vaporMassFlux()
        htf = self.p.heatTransferCoefficient()*self.p.mesh.interfaceArea[-1]
        #------------------------------------------------

        #-----------------interior nodes are easy---------------------
        #get coeffs for the interior nodes
        self.coeffs = float32(self.p.mesh.weightedAdjacencyMatrix(
                    self.p.conductivity,[0,htf]))#conductivity forms the basis for the coeffs
        self.coeffs[1]+=tm#self coeff also contains a mass term
        #-------------------------------------------------------------


        #-------------------rhs is also simple-----------------------------
        #rhs of inner nodes is the thermal mass from the last step
        self.rhs = tm*self.lastSoln()

        #account for heat in due to convection, loss due to radiation and vaporization
        conv = self.coeffs[-1,-1]*self.boundaryTemperature()
        vap = vmf*self.p.heatOfVaporization

        if self.p.plasmaTemperature>=self.p.surfaceTemperature:
            vap = max([0,min([vap,
                       -conv-diff(self.p.T[-2:])*self.p.conductivity[-1]*( abs(diff(1/self.p.mesh.xNode[-2:]))**-1 )*(4*pi)])])
#            if vmf>2e-7:
#                pdb.set_trace()

		
        sigma = 5.67e-8
        eps = .3
        if self.radiationOn:
            rad = sigma*eps*self.p.mesh.interfaceArea[-1]*(self.p.surfaceTemperature**4 - 300**4)
        else:
            rad = 0.0

#            with open('radTestDump', 'a') as f:
#                f.write('%.2e %.2e %.2e %.2e\n'%(htf,conv, rad, vap))


        self.rhs[-1]-= conv#last node's rhs has a known heat flux from convection/conduction
        self.rhs[-1]-= rad
        self.rhs[-1]-= vap


        #print(conv, rad, self.rhs[-1])
        #------------------------------------------------

		
        #-------------------mop up!-----------------------------
        #coeffs and rhs must be 3-d
        self.coeffs = expand_dims(self.coeffs,0)
        self.rhs = expand_dims(self.rhs,0)
        #------------------------------------------------

        
        #------------------debug statements
        #print str(self.coeffs[0,-1,-1])+'   '+str(self.boundaryTemperature())
        #print self.lastSoln()
        #print self.boundaryTemperature()
        #-----------------------------------


    def updateBCs(self):
        """
        move the particle through the plasma field, creating a new boundary temperature for the solver
        """
        #moving the particle through the field updates the plasma temperature
        self.p.fly()
        self.setCurrentSolution(self.p.T)
        self.bcHistory+=[self.p.plasmaTemperature]

    def plotSoln(self):
        """
        surface plot of the particle's temperature profile
        """
        T = self.getSoln()
        t,r = self.getSolnCoords()

        f = mlm.fig('theAnswer')
        ml.clf(f)
        mlm.mesh3(t,r,T,f)


        mlm.axe([ 't, us','r, um',  'T, K'], mlm.rgs(t,r,T),f)
        mlm.out()
        return f

    def getSoln(self):
        return array(self.solvedSteps)[:,0,:]
    def getTrajectory(self):
        return transpose(self.p.X[:self.nSteps(),:]).copy()
    def getSolnCoords(self):
        dims = array(self.solvedSteps)[:,0,:].shape
        r = []
        for m in self.p.meshHistory:
            r+= [m.xNode]


        r = double(int32(array(r[:-1])*1e7))/10

        t = double(int32(cumsum(self.p.thermalTimestepHistory)*1e7))/10
        t = t[...,newaxis]*ones(dims)

        return (t,r)
    def getEnthalpy(self):
        return sum(self.p.mass*self.p.enthalpy())

    def nSteps(self):
        return self.p.stepCount-1


class explicitParticleSolver(particleSolver):
    def __init__(self, particle = part.particle(thermalTimestep = 1e-6),
                            xMax = .1,
                            maxSteps = 1000,
                            verbose = False,
                            sema = None):
            """
            initialize using the marching solver's init method and the particle's number of nodes
            """
            #whether to print extra messages
            self.verbose = verbose

            #set the particle datamember
            self.p = particle

            #stop the sim if we pass xMax
            self.xMax = xMax

            #stop the sim if we

            self.maxSteps = fix(double(self.p.X.shape[0])/self.p.timestepRatio)
            sol.marchingSolver.__init__(self, self.p.H, self.p.mesh.nNodes,
                                        nVars = 1, verbose = verbose, sema = sema)#last node doesn't count since it's the boundary temperature



            self.convHistory = [0]
            self.radHistory = [0]
            self.vapHistory = [0]


            #debugging messages
            if self.verbose:
                print 'injection at position '+ str(self.p.position())
                print 'with velocity '+ str(self.p.velocity())
                #print 'at temperature '+ str(self.currentSoln())
                print 'ambient temperature is ' + str(self.p.plasmaTemperature)
                print 'with velocity ' + str(self.p.plasmaVelocity)
                
            self.radiationOn = True
    def solveStep(self):
        """
        like marchingSolver.solveStep, but
        update the particle's temperature mesh at every step
        """

        i = 0
#        alpha = double(self.alpha)
        while self.stepError()>self.stepErrorThres and i<self.maxIter or i==0:
            self.updateCoeffs()
            self.prevIterSoln = self.soln.copy()

            #nVars is actually just 1, for this class
            for j in arange(self.nVars):
                try:
#                    if abs(sum(self.coeffs[j][2]))>.001:
#                        pdb.set_trace()
                    deltaH = (self.tdmaSolve(self.coeffs[j], self.rhs[j])-self.soln[j])
                except ValueError:
                    pdb.set_trace()
                self.soln[j] += self.alpha*deltaH
#                deltaH = self.iterEnthalpyChange()
#                hN = self.lastSoln()[-1]+deltaH[-1]
#                Ts = self.p.material.Temperature(hN)


#                if (Tvap-Ts)<100:
#                    hBoil = self.p.material.Enthalpy(Tvap-100)
#                    self.soln[j] = self.lastSoln()+((hBoil-self.lastSoln()[-1])/deltaH[-1])*deltaH
#                else:
#                    self.soln[j] = self.lastSoln()+ alpha*deltaH
            #pdb.set_trace()
#            try:
            self.p.setEnthalpy(self.currentSoln())
#            except ValueError:
#                pdb.set_trace()
                
            i+=1
        #print 'used %d iter'%i
        if i>=self.maxIter:
            #print '%.2e, %.2e, %.2e'%(self.q_conv, self.q_rad, self.q_vap)
            print 'iteration limit reached'

        self.convHistory += [0]
        self.radHistory+= [0]
        self.vapHistory+= [0]

        self.itercount+=[i]

    def iterEnthalpyChange(self):
        q = zeros_like(self.q_cond)
        m = self.p.thermalMassConstant
        q[:-1]+=self.q_cond[1:]
        q-=self.q_cond
        q[-1]+=self.q_conv-self.q_rad-self.q_vap
        deltah = (1./m)*(q)

        

        return deltah

    def updateCoeffs(self):
        """
        use the particle's properties, current temperature,
        and mesh to get the coeffs
        """
        
        #alias for convenience
        m = self.p.mesh
        r = m.xNode
        R = m.xInt[-1]
        V = m.volume()
        T = self.p.T
        H = self.p.H
               
        #-----shell conductivities and explicit conduction term--------
        self.sigma = r_[0,4*pi*self.p.conductivity[1:]/(1./r[:-1]-1./r[1:])]
        self.q_cond = self.sigma*(T-r_[T[0],T[:-1]])
        
        
        cp = self.p.specificHeat()
        cpBar = 1./(V+r_[0,V[:-1]])
        cpBar *= V*cp+r_[0,V[:-1]*cp[:-1]]
        self.sigmah = self.sigma/cpBar
        self.q_condh = self.sigmah*(H-r_[0,H[:-1]])
        #-------------------------------------------------------------
        
        

        #################convective heating#########################
        htf = self.p.heatTransferCoefficient()*self.p.mesh.interfaceArea[-1]
        #surface temperature 
        self.Ts = self.p.T[-1]+(R-r[-1])*diff(self.p.T[-2:])/diff(r[-2:])
        #convective heating in W
        self.q_conv = htf*(self.p.plasmaTemperature-self.Ts)
        ############################################################

        #>>>>>>>>>>>>>>radiative heat loss<<<<<<<<<<<<<<<<<<<<<<<<<<
        sigmaRad = 5.67e-8
        eps = .6
        #radiative heating in W
        if self.radiationOn:        
            self.q_rad = sigmaRad*eps*m.interfaceArea[-1]*(self.Ts**4 - 300**4)
        else:
            self.q_rad = 0.0
        #>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        #::::::::::::::::::: evaporative cooling :::::::::::::::::::::
        vmf = self.p.vaporMassFlux()
        q_boil = self.q_conv-self.q_rad
        Tvap = self.p.vaporizationTemperature
        
        try:        
            q_vap_old = self.vapHistory[-2]
        except IndexError:
            q_vap_old = 0.0
            
        if self.Ts>(Tvap-100) and self.Ts<=(Tvap-50) and self.p.plasmaTemperature>self.Ts:
            q_boil = self.q_conv-self.q_rad
            q_nonBoil =  vmf*self.p.heatOfVaporization
            self.q_vap = (q_boil+q_nonBoil)*.5 
            
        elif self.Ts>(Tvap-100) and self.p.plasmaTemperature>self.Ts:
            self.q_vap = self.q_conv-self.q_rad
            self.p.vmf = self.q_vap/self.p.heatOfVaporization
        else:
            self.q_vap = vmf*self.p.heatOfVaporization
        
        self.q_vap = .4*q_vap_old+.6*min([self.q_vap, max([.993*(self.q_conv-self.q_rad),0])])
        
        self.p.vmf = self.q_vap/self.p.heatOfVaporization
        
        #self.B_h = self.p.B_h_0()*(1-(self.q_conv-self.q_rad-self.q_vap)/self.q_conv)
        #print "BhFactor %.2f"%(1+self.B_h)**-.7
        

        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        
        #================== bookkeeping=====================
        self.convHistory[-1] = double(self.q_conv)
        self.radHistory[-1] = double(self.q_rad)
        self.vapHistory[-1] = double(self.q_vap)
        #====================================================

        
        #$$$$$$$$$$$$$$$$$$$ RHS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        #right hand side has a thermal mass term from the previous timestep
        self.rhs[0] = self.p.density*V*self.lastSoln()/self.p.thermalTimestep
        #treat the boundary using source terms
        self.rhs[0][-1]+=self.q_conv-self.q_rad-self.q_vap
        
        #add the difference between (lagged) true conduction and 
        #enthalpy based conduction if we are using that correction
        self.explicitCorrection = True
        if self.explicitCorrection:
            self.rhs[0][:-1]+= self.q_cond[1:]-self.q_condh[1:]
            self.rhs[0]-= self.q_cond-self.q_condh
        
        #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        
        
        #//////////////////////// matrix coefficients///////////////////
        #diagonal coeff has the mass term and sum of forward and backward coeffs
        self.coeffs[0, 1, :] = (self.p.density*V)/self.p.thermalTimestep+self.sigmah+r_[self.sigmah[1:], 0]
        self.coeffs[0, 0, :] = -self.sigmah
        self.coeffs[0, 2, :] = -r_[self.sigmah[1:], 0]
        
#        pdb.set_trace()
        
        #///////////////////////////////////////////////////////////////
        
    def updateBCs(self):
        """
        move the particle through the plasma field, creating a new boundary temperature for the solver
        """
        #moving the particle through the field updates the plasma temperature
        self.p.fly()
        self.setCurrentSolution(self.p.h)
        self.bcHistory+=[self.p.plasmaTemperature]

    def plotSoln(self):
        """
        surface plot of the particle's temperature profile
        """
        h = self.getSoln()
        t,r = self.getSolnCoords()

        f = mlm.fig('theAnswer')
        ml.clf(f)
        mlm.mesh3(t,r,h,f)


        mlm.axe([ 't, us','r, um',  'h, J/kg'], mlm.rgs(t,r,h),f)
        mlm.out()
        return f
        
    def getSoln(self):
        h= array(self.solvedSteps)[:,0,:]
        T = self.p.material.Temperature(h)
        return T
    def getSolnH(self):
        h= array(self.solvedSteps)[:,0,:]
       
        return h



if __name__== "__main__":
    solverTestOn = False
    close('all')
    if solverTestOn:
        with open('radTestDump', 'w') as f:
            pass
        maxSteps = 5000

        #start = time.time()
        PS = particleSolver(part.particle(diameter = array([.1e-6, 35e-6]),
                                                        T =300,
                                                        thermalTimestep = 1e-5,
                                                        kinematicTimestep = 1e-5,
                                                        x0 = array([0.006,.008,0.]),
                                                        v0 = array([0.,-30.,0.]),
                                                        nNodes = 50,
                                                        nStep = maxSteps) )
        #create the plasmaField
        try:
            PF.density(array([0,0,0]))
        except NameError:
            PF = part.plasmaField(fastOn = True, temperatureScaling = 1.35)
       # pf = part.plasmaField(fastOn = True, temperatureScaling = 1.35)
#        PF = part.lavaField('sharpTemp_600_38.pkl',
#                            fastOn = True, gasPickle = None, temperatureScaling = 1.)
#
        #PS.p.inject(PF,x0 = array([0.001,.004,0.]), v0 = array([0.,-10.,0.]))
        PS.p.inject(PF)

        #//////////////// Timestepping options //////////////////////////
        PS.p.adaptiveTimestepping = True
        PS.p.splitScales = True
        d0 = PS.p.diameter()



        PS.reset()
        start = time.time()
        print 'solving...'
        #cProfile.run('PS.solve()', 'prof')
        PS.solve()
        #pdb.set_trace()
        print 'solved!'
        elapsed = (time.time()-start)
        d1 = PS.p.diameter()

        #f = PS.plotSoln()
        T = PS.getSoln()
        t,r = PS.getSolnCoords()
        x = (PS.getTrajectory()[:-1])
        D = dict(zip(['T','t','r','x'],[T,t,r,x]))

        g = part.gm.defaultMaterial

        H = sum(g.Enthalpy(T)*PS.p.mass, axis = 1)

        fig()
        plot(t[:-1,0], diff(H)/diff(t[:,0]*1e-6))
        xlab('time, microseconds')
        ylab('power, W')

        print 'distance traveled %d'%(PS.p.axialPosition()*100)
        print 'd1/d0: %f'%(d1/d0)
        print 'time elapsed: '+str(elapsed)
        print 'time per step: ' +str(elapsed/len(PS.itercount))
        #print 'mean iterations per step:  '+str(mean(array(PS.itercount)))
        #print 'max iterations per step:   '+str(max(PS.itercount))
        print 'mean final temperature: '+str(mean(T[-1]))


#        pstats.Stats('prof').strip_dirs().sort_stats('time').print_stats(15)
#        pstats.Stats('prof').strip_dirs().sort_stats('cumulative').print_stats(15)

        PS.p.plasmaField.showField()
        subplot(2,1,1)
        plot(x[0], x[1],'.', linewidth = 3, color = 'k')
        subplot(2,1,2)
        plot(x[0], x[1],'.', linewidth = 3, color = 'k')

        # figure(figsize = (20,22))
        # plot(PS.E, '.')


        stepSize =pnorm(diff(PS.p.X,axis=0),dim = 1)[:PS.nSteps()]
        thermalStepSize = diff(array(PS.bcHistory))

        figure(figsize = (15,20))
        subplot(3,1,1)
        plot(array(PS.p.thermalTimestepHistory)*1e6,'.')
        xlabel('step')
        ylabel('timestep, $\mu$s')
        subplot(3,1,2)
        plot(stepSize*1e3,'.')
        xlabel('step')
        ylabel('spatial step size, mm')
        subplot(3,1,3)
        plot(thermalStepSize,'.')
        xlabel('step')
        ylabel('thermal step size, K')

        fig()
        plot(t[:,0], array(PS.bcHistory), linewidth = 2,label = 'plasma temp')
        plot(t[:,0],T[:,-1], linewidth = 2,label = 'surface temp' )
        legend()
        xlabel('Time, $\mu$s', fontsize = 14)
        sf('surfaceTempAndBC_lava')

        #print [mean(Dstar['T'][-1]), mean(Dslow['T'][-1]), mean(Dfast['T'][-1]), mean(Dsplit['T'][-1])]
        fig()
        plot(t[:,0],(T[:,-1]),'-o', label = 'surface')
        plot(t[:,0],(T[:,0]),'-o', label = 'center')
        xlabel('Time, $\mu$s', fontsize = 14)
        ylabel('Temperature', fontsize = 14)

        legend()

#        with open('radTestDump', 'r') as f:
#            radOut = float64([float64(L.split()) for L in f.readlines()])
#        radLabs = ['heat transfer coeff*Area', 'convective heat transfer', 'radiative heat transfer', 'evaporative heat transfer']
#        fig()
#        for i in range(radOut.shape[1]):
#            subplot(2,2,i+1)
#            plot(radOut[:,i], linewidth = 2)
#            tit(radLabs[i])
#
#        with open('vaporModelDump', 'r') as f:
#            vapOut = float64([float64(L.split()) for L in f.readlines()])
#        vapLabs = ['Tsurf', 'Tplas', 'Sherwood Nr', 'vapor Pressure, Pa', 'film vapor fraction', 'B', 'diffusion Coeff', 'vapor mass flux, cubic mics/microsecond', 'fVap']
#        fig()
#        for i in range(vapOut.shape[1]):
#            subplot(3,3,i+1)
#            plot(vapOut[:,i], linewidth = 2)
#            title(vapLabs[i])
        #show()
