import pysces
import pymc
from numpy import *
from numpy.random import randn 
from matplotlib import *
import pdb
from pylab import *
import time
import cPickle as pkl
from scipy import stats
from scipy.interpolate import interp1d
from inverseModel import forwardModel, inverseModel
from multiprocessing import Process
import scipy

class constrainedPysMod(pysces.PyscesModel.PysMod):  
    def Simulate(self, f = None,userinit=0, verbose = False ):
        """
        solve the ode of the system with constrained species values
        f(t) must return a vector the same length as the number of constrained species, 
        """  
        
        self.f = f
        
        specArr = array(self.species, dtype = 'object')
        
        if self.f!=None:
            if type(self.f.constrainedSpeciesName) == type('string'):
                self.constrainedSpeciesIdx = [flatnonzero(specArr==self.f.constrainedSpeciesName)[0]]
            else:            
                self.constrainedSpeciesIdx = [flatnonzero(specArr==sp)[0] for sp in self.f.constrainedSpeciesName]
        R = pysces.PyscesModel.PysMod.Simulate(self, userinit, verbose = verbose)
        return R
        
    def getConstrainedSpecies(self):
        return [self.f(self._TIME_)]
        
    def _EvalODE(self,s,Vtemp):
        """
        overwrite this method so we can enforce the constraint at the end

        _EvalODE(s,Vtemp)

        Core ODE evaluation routine evaluating the reduced set of ODE's. Depending on whether mass conservation is present or not either N*v (_EvalREq) or Nr*v (_EvalREq2) is used to automatically generate the ODE's.
        Returns sdot.

        Arguments:

        s: species vector
        Vtemp: rate vector

        """
        if self.__HAS_RATE_RULES__:
            S, self.__rrule__ = numpy.split(s,[self.Nrmatrix.shape[0]])
        else:
            S = s.copy()
        
        if self.__HAS_MOIETY_CONSERVATION__ == True:
            for x in range(len(S)):
                self.__SALL__[self.nrmatrix_row[x]] = S[x]      
            for x in range(len(self.lzeromatrix_row)):
                self.__SALL__[self.lzeromatrix_row[x]] = self.__tvec_a__[x] + numpy.add.reduce(self.__lzeromatrix__[x,:]*S)
            self.__vvec__ = self._EvalREq(self.__SALL__,Vtemp)
            s_dot = numpy.add.reduce(self.__nrmatrix__*self.__vvec__,1)
        else:
            self.__SALL__ = S
            self.__vvec__ = self._EvalREq(self.__SALL__,Vtemp)
            s_dot = numpy.add.reduce(self.__nmatrix__*self.__vvec__,1)
        if self.__HAS_RATE_RULES__:
            s_dot = numpy.concatenate([s_dot, self.__rrule__]).copy()
         
        
        #compute what s_dot should be for the constrained species
        if self.f!=None:
            f = self.getConstrainedSpecies()
            for i, idx in enumerate(self.constrainedSpeciesIdx):
                cidx = flatnonzero(self.nrmatrix_row==idx)[0]
                s_dot[cidx] = f[i]
        return s_dot
 

class timeSeriesData(object):
    def __init__(self, datafile):
        pass
    
    def __call__(self):
        return self.D.copy()
    
    def weightedDistance(self, D):
        d = ((self.D-D)/self.sigma)**2
        return sum(d)
    
    def logp(self, D):
        return self.logNormalizingConstant-self.weightedDistance(D)*.5
    def time(self):
        return self.t.copy()
        
class erk_timeseries(timeSeriesData):
    def __init__(self, data, N = 100):
        if type(data)==type(float64(arange(2))):
            t = data[0]
            D = data[1]
            sig = data[2]
        else:
            t = []
            D = []
            sig = []
            with open(data, 'r') as f:
                for L in f.readlines():
                    x = float64(L.split())
                    t+=[x[0]]*60          #assuming time given in minutes (only in the datafile)!
                    D+=[x[1]]
                    sig+=[x[2]]
            
            
        self.t = linspace(t[0], t[-1], N)     
        self.D = interp1d(t, D, kind = 'linear')(self.t)
        self.sigma = interp1d(t, sig, kind = 'linear')(self.t)
        
        assert not any(self.sigma<1e-15), "small std dev detected!"
        self.logNormalizingConstant = -.5*(self.sigma.shape[0]*log(2*pi)+log(prod(self.sigma)))
        

        
class ret_erk_experiment(object):
    def __init__(self, datafile, speciesDictFile = None):
        if speciesDictFile == None:
            self.speciesDictFile = 'ret_erk_species.dic'
        else:
            self.speciesDictFile = speciesDictFile 
        try:
            self.spName = self.loadSpeciesDict()[3]['RetStar']
        except KeyError:
            self.spName = self.loadSpeciesDict()[3]['R_star']
        
        
        with open(datafile, 'r') as f:
            self.L = f.readlines()
            self.currentLine = self.L[0]
            self.nextLine = self.L[1]
            
        D = []
        self.cnt = 1
        self.done = False            
        while not self.done:
            D+=[self.getSeries()]
            
        self.ERK = D[slice(1, None, 2)]
        self.RET = D[slice(0, None, 2)]
        self.t = self.ERK[0][0]
            
        self.D = [erk_timeseries(erk) for erk in self.ERK]
        self.F = [singleSpeciesConstraintFunction(self.spName, ret) for ret in self.RET]
    def logp(self, D):
        LP = 0
        for i,d in enumerate(self.D):
            LP+=d.logp(D[i])
        if isnan(LP):
            LP = -1e20
        return LP
            
    def val(self):
        return array([d() for d in self.D])
        
    def log_normConst(self):
        LNC = 0
        for d in self.D:
            LNC+=d.logNormalizingConstant
        return LNC
        
    def getSeries(self):
        S = []
        while not (double(self.nextLine.split()[0])<1e-12):
            S+=[float64(self.currentLine.split())]
            self.currentLine = str(self.nextLine)
            self.cnt+=1
            if self.cnt == (len(self.L)):
                S+= [float64(self.currentLine.split())]
                self.done = True
                return float64(S).T
            else:
                self.nextLine = self.L[self.cnt]
        
        S+=[float64(self.currentLine.split())]        
        self.currentLine = str(self.nextLine)
        self.cnt+=1
        self.nextLine = self.L[self.cnt]
        
        A = float64(S).T
        
        #time in minutes in the datafile!        
        A[0] = A[0]*60
        
        return A
        
    def loadSpeciesDict(self):
        idx = []
        s = []
        nm = []    
        with open(self.speciesDictFile, 'r') as f:
            for L in f.readlines():
                x = L.split()   
                s +=[x[0]]
                nm +=[x[1]]
                idx+= [self.parseIdx(x[0])]
                
        return (dict(zip(nm, idx)), dict(zip(idx, nm)), dict(zip(idx, s)), dict(zip(nm, s)))
    def parseIdx(self, speciesOrRateConstantString):
        idx = int(speciesOrRateConstantString[1:])-1
        assert idx>=0, "index must be semipositive!"
        return idx            
            
        
    def time(self):
        return self.t.copy()
    def erkTime(self):
        return self.D[0].time()
    def constraintFunctionList(self):
        return self.F
        

        
class constraintFunction(object):
    def __init__(self, constrainedSpeciesName, constraint):
        """
        constraint - 3 x nTimepoints array of (time, constraintValue, stdDev)
        """
        self.constrainedSpeciesName = constrainedSpeciesName
            
        self.t = constraint[0].copy()
        self.s = constraint[1].copy()
        self.stdDev = constraint[2].copy()
        
        self.initInterpolator()
        
    def __getstate__(self):
        """
        need to exclude unpickleable members
        """
        black =['f']
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
        gonna have to call the methods that reconstruct the unpickleable members
        """
        self.__dict__= D.copy()
        self.initInterpolator()
        
    def __call__(self, t):
        pass
        
    def initInterpolator(self):
        #expect the interpolator to be named f for compatability with __getstate__        
        self.f = None
 
        
class singleSpeciesConstraintFunction(constraintFunction):
    def __call__(self, t):
        return self.f(t)

    def initInterpolator(self):
        self.T = linspace(self.t[0], self.t[-1], 200)        
        self.S = interp1d(self.t, self.s, kind = 'cubic')(self.T)
        
        dSdT = diff(self.S)/diff(self.T)
        
        self.f = interp1d(r_['0,1', array([-1]), self.T[:-1]+diff(self.T)*.5, array([self.T[-1]*inf])],
                             r_['0,1', ones(1)*dSdT[0], dSdT, ones(1)*dSdT[-1]],
                                kind = 'linear')
                
    

class pyscesForwardModel(forwardModel):
    def __init__(self, 
                 pscFile, 
                 speciesDictFile, 
                 measuredSpeciesFile,
                 referenceParameterFile,
                 measuredTimes, 
                 constraintFunctionList
                 ):
        """
        initialize the pysces model
        
        pscFile - model file for pysces
        speciesDictFile - contains correspondence between the long species names
                            and the abstract species names, e.g.
                            s1 RET
                            s2 RET*
                            .
                            .
                            .
        measuredSpeciesFile - long names of the species that are measured in the data
        

        referenceParameterFile - file with reference rate constants and initial 
                                concentrations for those species with variable 
                                initialValues, e.g.
                                
                                k1 k1_0
                                k2 k2_0
                                .
                                .
                                .
                                speciesNameX speciesX_0
                                
        measuredTimes - vector of times, in s, where we have measured the species concentrations
        
        constraintFunctionList - each F[i] is a forcing function suitable for use with timeForcedPysMod
        """
        self.pscFile = pscFile
        self.referenceParameterFile = referenceParameterFile
        self.speciesDictFile = speciesDictFile
        self.measuredSpeciesFile = measuredSpeciesFile
        
        self.measuredTimes = measuredTimes
        
        self.F = constraintFunctionList
        
        self.initIdx()
        self.initReferenceValues()
        self.initPyscesModel()     
        self.initModelTimes()
        
    def modelDir(self):
        try:
            with open(os.path.join(os.getcwd(), 'model.dir'), 'rb') as f:
                md = f.readLine().strip().replace('\n', '')
        except:
            md = os.path.split(self.pscFile)[0]
        return md
        
    def initPyscesModel(self):
        self.pscMod = constrainedPysMod(os.path.split(self.pscFile)[1], self.modelDir())
        self.pscMod.__settings__['mode_sim_max_iter'] = 1
        self.pscMod.__settings__['lsoda_mxstep'] = 1500
        
                
    def initModelTimes(self):
        """
        use measuredTimes to set simulation start and end times
        """
        self.pscMod.sim_end = self.measuredTimes[-1]
        self.pscMod.sim_points = self.measuredTimes.shape[0]
        
    def initShortModelTimes(self):
        self.pscMod.sim_end = self.measuredTimes[-1]*1e-2
        self.pscMod.sim_points = 100
        
    def initIdx(self):
        self.measuredSpecies = []        
        with open(self.measuredSpeciesFile, 'r') as f:
            for L in f.readlines():
                self.measuredSpecies += L.split()
        #speciesDic:     speciesName->speciesNumber
        #idxDic:         speciesNumber->speciesName
        #sDic:           speciesName->numberedSpeciesName
        self.speciesDic, self.idxDic, bluh, self.sDic = self.loadSpeciesDict()
        
        #measuredList: [numberedSpecies1, numberedSpecies2, ...]          
        self.measuredList = [self.sDic[k] for k in self.measuredSpecies]
        self.midx = array([self.speciesDic[k] for k in self.measuredSpecies])
        

    def initReferenceValues(self, referenceParameterFile = None):
        if referenceParameterFile == None:
            referenceParameterFile = self.referenceParameterFile
        with open(referenceParameterFile, 'r') as f:
            K = []
            vk = []
            s = []
            vs = []
            for L in f.readlines():
                x = L.split()
                nm = x[0]
                val = double(x[1])
                if x[0][0]=='k':
                    K+=[nm]
                    vk+=[val]
                elif x[0][0]=='s':
                    s +=[nm]
                    vs+=[val]
        self.k0Dic = dict(zip(K, vk))
        self.s0Dic = dict(zip(s, vs))
        
        self.k0 = zeros(len(K))
        self.s0 = zeros(len(s))
        
        self.initNumberOfParameters()
        
        for k in K:
            self.k0[self.kidx(k)] = self.k0Dic[k]
        for k in s:
            self.s0[self.parseIdx(k)] = self.s0Dic[k]

        
    def initNumberOfParameters(self):
        self.nK = len(self.k0Dic.keys())
        self.nS = len(self.s0Dic.keys())
        self.nM = self.nK+self.nS
        self.M = zeros(self.nM)

              
    def kidx(self, spcStr):
        return self.parseIdx(spcStr)
    def sidx(self, spcStr):
        return self.parseIdx(spcStr)+self.nK
            
    def parseIdx(self, speciesOrRateConstantString):
        idx = int(speciesOrRateConstantString[1:])-1
        assert idx>=0, "index must be semipositive!"
        return idx
        
    def idx2k(self, idx):
        return 'k'+str(idx+1)
    def idx2s_init(self, idx):
        return 's'+str(idx-self.nK+1)+'_init'
        
    def idx2s(self, idx):
        return 's'+str(idx-self.nK+1)
        
    def loadSpeciesDict(self):
        idx = []
        s = []
        nm = []    
        with open(self.speciesDictFile, 'r') as f:
            for L in f.readlines():
                x = L.split()   
                s +=[x[0]]
                nm +=[x[1]]
                idx+= [self.parseIdx(x[0])]
                
        return (dict(zip(nm, idx)), dict(zip(idx, nm)), dict(zip(idx, s)), dict(zip(nm, s)))
    def __getstate__(self):
        """
        need to exclude unpickleable members
        """
        black =['pscMod']
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
        gonna have to call the methods that reconstruct the unpickleable members
        """
        self.__dict__= D.copy()
        self.initPyscesModel()
        self.initModelTimes()
        
    def __call__(self, k, s):
        """
        solve the forward model for the parameters M
        """
        self.set_m(r_[k,s])
        self.solve()
        self.cook()
        return self.D.copy()
        
    def solve(self):
        """
        for all sets of forcing concentrations, simulate the system forward in time 
        and record the resulting timeseries
        """
        nSim = len(self.F)
        self.S = zeros((nSim, self.nS ,self.pscMod.sim_points))
        self.D = zeros((nSim, self.midx.shape[0], self.pscMod.sim_points))
#        print "solving..."
        for i in range(nSim):
#            print "    solving %d"%i            
            #pre-solve to avoid stiffness issues
            self.initShortModelTimes()
            self.pscMod.__settings__['mode_sim_init'] = 0
            start = time.time()
            self.pscMod.Simulate(self.F[i])
#            print time.time()-start
                        
            #switch back to the real modeling times
            self.initModelTimes()
            self.pscMod.__settings__['mode_sim_init'] = 3
            start = time.time()
            self.pscMod.Simulate(self.F[i])
#            print time.time()-start
            self.D[i], self.S[i] = self.observe()
#        print "done"

        
            
    def observe(self):
        S = self.pscMod.data_sim.getSpecies().T[1:]
        
        return S[self.midx].copy(), S.copy()
        
        
        
    def cook(self):
        """
        after solving, pick out the measured species and do any processing to
        make the result comparable to any measured data
        """
        for i in range(self.D.shape[0]):
#            plot(self.D[i][0].copy())
            tmp = self.D[i].copy()
            mx = amax(self.D[i], axis = -1)
            tmp /= mx[..., newaxis]
            d0 = tmp[:,0].copy()
            self.D[i] = tmp/d0[..., newaxis]
            
#            plot(self.D[i][0].copy())
#            pdb.set_trace()
            
        
    def m(self):
        """
        assemble variable initial concentrations and rate constants into a single vector of 
        model parameters
        """
        for k in self.k0Dic.keys():
            self.M[self.kidx(k)] = log2(self.pscMod.__dict__[k]/self.k0Dic[k])
        for k in self.s0Dic.keys():
            self.M[self.sidx(k)] = log2(self.pscMod.__dict__[k]/self.s0Dic[k])
        return self.M
            
    def set_m(self, M):
        """
        given a vector of modelParameters m, parse the vector and set the proper fields of
        pscMod
        """
        for i in range(self.nK):
            k = self.idx2k(i)            
            self.pscMod.__dict__[k] = self.k0Dic[k]*2**M[i]
        for i in range(self.nK, self.nM):
            ki = self.idx2s_init(i)      
            k = self.idx2s(i)      
            self.pscMod.__dict__[ki] = self.s0Dic[k]*2**M[i]
            
    def deLogM(self, Mstar):
        M = zeros_like(Mstar)
        for i in range(self.nK):
            k = self.idx2k(i)            
            M[i] = self.k0Dic[k]*2**Mstar[i]
        for i in range(self.nK, self.nM):
            k = self.idx2s(i)    
            M[i] = self.s0Dic[k]*2**Mstar[i]   
        return M
        
    def splitM(self, M):
        return (M[:self.nK], M[self.nK:])
    def unsplitM(self, k, s0):
        return r_['0,1', k, s0]
        
    def t(self):
        return self.measuredTimes.copy()
            
class inversePyscesModel(inverseModel):
    def __init__(self, forwardModel, observedData, 
                 logRateConstantStdDev = 2.0, 
                 logInitialConcentrationsStdDev = 2.0, 
                 verbose = False):
        print "inverse model initilizing..."
        inverseModel.__init__(self, forwardModel, observedData)
        
         #suppress printing of LSODA message
        self.verbose = verbose
        if not self.verbose:
            self.FM.pscMod.__settings__["lsoda_mesg"] = 0
        
        self.logRateConstantStdDev = logRateConstantStdDev
        self.logInitialConcentrationsStdDev = logInitialConcentrationsStdDev
        
        print "mcmc model initializing..."
        self.init_mcModel()
        
        
    def __call__(self, N, burn = 0, thin = 1, pklName = 'blah'):
        self.pklName = pklName    
        inverseModel.__call__(self, N, burn = burn, thin = thin)
        
        
    def __getstate__(self):
        """
        need to exclude unpickleable members
        """
        black =['mcmc', 'mcModel']
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
        gonna have to call the methods that reconstruct the unpickleable members
        """
        self.__dict__= D.copy()
        self.init_mcModel()
        
    def trace(self):
        self.tK = self.mcmc.trace('k')[:].T
        self.tS = self.mcmc.trace('s')[:].T
        self.tD = self.mcmc.trace('D_mod')[:]
        self.tLogP = array([self.modelLikelihood(None, self.tD[i]) for i in range(self.tD.shape[0])])
        self.save(self.pklName)
    
    def bounceReferenceParameters(self, kStar, s0Star, fname):
        Mstar = self.FM.unsplitM(kStar, s0Star)
        M = self.FM.deLogM(Mstar)
        K, S0 = self.FM.splitM(M)
        with open(fname, 'wb') as f:
            [f.write('k%d %.2e\n'%(i+1,k)) for i,k in enumerate(K)]
            [f.write('s%d %.2e\n'%(i+1,s)) for i,s in enumerate(S0)]
        
        
    def init_mcModel(self):
        #############init pymc model######################
        print "    stochastics initializing"
        #priors for rate constants
        k = pymc.Normal('k',
                         zeros(self.FM.nK),
                         1./self.logRateConstantStdDev**2)
        print "        k is done"
                         
        #priors for initial concentrations
        s = pymc.Normal('s',
                         zeros(self.FM.nS),
                         1./self.logInitialConcentrationsStdDev**2)
        print "        s is done"
        #modeled species concentrations
        D_mod = pymc.Deterministic( self.FM,
                                  'modeled species concentrations', 
                                  'D_mod' , 
                                  {'k':k, 's':s},
                                  verbose = 0, 
                                  plot = False)
        print "        D_mod is done"
        #observed species concentrations                                  
        D_obs = pymc.Stochastic(self.modelLikelihood,
                               'measured species concentrations',
                               'D_obs',
                               {'D_mod': D_mod},
                               value = self.D_obs.val(),
                               verbose = 0,
                               plot = False,
                               observed = True)
        print "        D_obs is done"
        
        print "    mcModel initializing"                          
        self.mcModel = pymc.Model([k,s, D_mod, D_obs])
        print "    mcmc initializing"                          
        self.mcmc = pymc.MCMC(self.mcModel)
        
    def save(self, name):
        with open(name, 'wb') as f:
            pkl.dump(self, f, -1)
        
class inverseModelPool(object):
    def __init__(self, invModList, pklName):
        self.invModList = invModList
        self.nm = pklName
        self.writeModelHeader()
     
    def writeModelHeader(self):
        with open(os.path.join(os.getcwd(), 'model.dir'), 'wb') as f:
            f.write(os.path.split(self.nm)[0])
            
        print os.getcwd()
            
    def __call__(self, N, burn = 0, thin = 1):
        P = [Process(target = IM, args = (N, burn, thin, self.pklName(i))) 
                            for i,IM in enumerate(self.invModList)]
        [p.start() for p in P]
        [p.join() for p in P]
        self.trace()
        
    def trace(self):
        self.load()
        self.tK = self.invModList[0].tK.copy()
        self.tS = self.invModList[0].tS.copy()
        self.tD = self.invModList[0].tD.copy()
        self.tLogP = self.invModList[0].tLogP.copy()
        K = ['tK', 'tS', 'tD', 'tLogP']
        for i in range(1, len(self.invModList)):
            for k in K:
                if k=='tD':
                    self.__dict__[k] = r_['0', self.__dict__[k], self.invModList[i].__dict__[k]]
                else:
                    self.__dict__[k] = r_['-1', self.__dict__[k], self.invModList[i].__dict__[k]]
        
    def pklName(self, i):
        return self.nm + "."+str(i)+".pkl"
        
    def load(self):
        self.invModList = []        
        for i in range(self.nPkl()):
            with open(self.pklName(i), 'rb') as f:
                self.invModList += [pkl.load(f)]

    def nPkl(self):
        cnt = 0
        while True:
            if os.path.isfile(self.pklName(cnt)):
                cnt+=1
            else:
                break
        if cnt>0:
            return cnt
        else:
            pdb.set_trace()
            raise ValueError, 'no files found!'

class savedInverseModelPool(inverseModelPool):
    def __init__(self, pklName):
        self.nm = pklName
        self.writeModelHeader()

            
        
        
if __name__== "__main__":
    ########################################################    
    ###################MCMC parameters######################
    ########################################################    
    N = 2
    burn = 0
    thin = 1
    
    logRateConstantStdDev = 2.0
    logInitialConcentrationsStdDev = 2.0
    ########################################################
    
    force = False
    nP = 2
    
    
    modelPath = '/media/ubuntuData/CODE/pysces/ret_erk_scratch_model'
    if not os.path.isdir(modelPath):
        modelPath = 'D:\\CODE\\pysces\\ret_erk_scratch_model'
    fileDic = {
        'pklName':'test_pERK',
        'dataFile':'ret_erk_small_timeseries.csv',
        'pscFile':'R-ERK_s.psc', 
        'speciesDictFile':'R-ERK_s_species.dic',
        'measuredSpeciesFile':'R-ERK_s_MeasuredSpecies.txt',
        'referenceParameterFile':'R-ERK_s_referenceParameters.dic'
        }
    for k in fileDic.keys():
        fileDic[k] = os.path.join(modelPath, fileDic[k])
    
    
    
    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    #             initialize data
    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    print "data initializing..."
    D = [ret_erk_experiment(fileDic['dataFile'], fileDic['speciesDictFile'])for i in range(nP)]    
    measuredSpecies = ['RetStar', 'ppERKHalf']
    print "done..."
    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    
    #-------------------------------------------------------
    #             initialize the pysces model
    #-------------------------------------------------------


    
    print "initializing forward models..."
    FM = [pyscesForwardModel(fileDic['pscFile'], 
                 fileDic['speciesDictFile'], 
                 fileDic['measuredSpeciesFile'],
                 fileDic['referenceParameterFile'],
                 D[i].erkTime(),
                 D[i].constraintFunctionList()
                 ) for i in range(nP)]
    print "done. prepare to be patient."
    #-------------------------------------------------------
    




    #////////////////////////////////////////////////////////
    #               get initial parameters
    #////////////////////////////////////////////////////////
        #determine number of k's    
    k0 = FM[0].k0
    s0 = FM[0].s0
    #////////////////////////////////////////////////////////
    
    
    ########################################################    
    ########### read from disk or run simulation############
    ########################################################
    if not force:
        try:
            IMPool = savedInverseModelPool(fileDic['pklName'])
            IMPool.trace()
        except:
            pdb.set_trace()
            print "load failed! forcing computation."
            force = True
        
    if force:            
        print "computation is forced! will now initialize inverse model list..."        
        #initialize inverse model
        IMList = []
        
        for i in range(nP):
            print "initializing %d of %d"%(i+1,nP)
            IMList += [inversePyscesModel(FM[i], 
                                    D[i], 
                                    logRateConstantStdDev = logRateConstantStdDev, 
                                    logInitialConcentrationsStdDev = logInitialConcentrationsStdDev, 
                                    verbose = False)]
        print "list initialized, setting up the pool..."
        IMPool = inverseModelPool(IMList, fileDic['pklName'])
        print "done. now for the tough part..."
        
        #sample from the posterior    
        start = time.time()    
        print "here we go!"
        IMPool(N, burn = burn, thin = thin)
        elapsed = time.time()-start
        print 'sampled %d in %.2f s: %.2e samples/s'%(N*nP,elapsed,nP*N/elapsed  )
        

        
    #break out results for easy plotting
    tK = IMPool.tK
    tS = IMPool.tS
    tD = IMPool.tD
    tLogP = IMPool.tLogP
    logNormConst = IMPool.invModList[0].D_obs.log_normConst()
    
    ########################################################
    
    
    
    
    #////////////////////////////////////////////////////////
    #               posterior histograms
    #////////////////////////////////////////////////////////
    subN = (3,3)
    nFig = ceil(k0.shape[0]/prod(subN)) 
    x = linspace(-3*logRateConstantStdDev, 3*logRateConstantStdDev, 200)
    for i in range(tK.shape[0]):
        idx = mod(i, prod(subN))        
        if idx==0:
            figure()
        subplot(subN[0],subN[1],idx+1)
        
        hist(tK[i],50, label = 'posterior', normed = True)
        plot(x, stats.norm.pdf(x,loc = 0, scale = logRateConstantStdDev ), label = 'prior')
        title('k%d, truth = %.1f'%(i+1, k0[i]))
        xlabel(r'log(k/k$_{true}$)')
     #////////////////////////////////////////////////////////   
        
    
    
    #-------------------------------------------------------
    #             random models
    #------------------------------------------------------
    didx = randint(0,tD.shape[0], 4)
    print "randomly selected trace points are:"    
    print didx
    figure()
    for i, idx in enumerate(didx):
        subplot(2,2,i+1) 
        plot(D[0].erkTime(), D[0].val().T, 'k')
        plot(D[0].erkTime(), squeeze(tD[idx]).T, 'y')
        title('log(P(D|m)) = %.2f'%tLogP[idx])
    #------------------------------------------------------
    
    
    ########################################################    
    ################### correlation coeffs##################
    ########################################################
    figure()
    pcolor(arange(1.,tK.shape[0]+1)-.5, 
           arange(1.,tK.shape[0]+1)-.5,
           corrcoef(tK)-eye(tK.shape[0]))
    colorbar()
    xlabel('constant index')
    ylabel('constant index')
    title('correlation coefficient between rate constants')
    ########################################################
    
    #.......................................................
    #................goodness of fit histograms.............
    #.......................................................
    figure()
    hist(tLogP, 100)
    xlabel('log(P(D|m))')
    title('goodness of fit, max possible log(P(D|m)) = %.2e'%logNormConst)
