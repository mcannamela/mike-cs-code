# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 17:02:30 2011

@author: wichtelwesen
"""
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

rateConstantPctError = .5
sigmaExp = .75
logkSigma = .75

class pyscesInverseModel(object):
    def __init__(self, pyscesMod, experimentalData, verbose = False):
        
        self.pcMod = pyscesMod#forward model        
        
        self.D = experimentalData#data        
        
        self.nK = 0
        for p in self.pcMod.parameters:
            if p[0]=='k':
                self.nK+=1
        self.k_ = zeros(self.nK)     
        
        self.k0 = self.k().copy()#initial parameter values
        
        self.verbose = verbose
        #suppress printing of LSODA message
        if not self.verbose:
            self.pcMod.__settings__["lsoda_mesg"] = 0             
        
#        self.sigma = matrix(identity(self.D.flatten().shape[0])*sigmaExp)
#        self.sigmaInv = linalg.inv(self.sigma)
#        self.detSigma = linalg.det(self.sigma)
#        u, s, v = linalg.svd(self.sigma)
#        self.rankSigma = sum(s > 1e-12)
#        self.logNormConst = -.5*(self.rankSigma*log(2*pi)+ log(self.detSigma))


        self.sigma = sqrt(self.D.size*sigmaExp**2)
        self.logNormConst = log((2*pi*self.sigma**2)**-.5)
        
        
        self.init_mcModel()
        
    def __call__(self, N, burn = 10000, thin = 2):
        self.mcmc.sample(iter = burn+thin*N, burn = burn, thin = thin)

                    
    def init_mcModel(self):
        self.logpTrace  = zeros(10)  
        self.cnt = 0
        #############init pymc model######################
        
        #priors for rate constants
        k = pymc.Normal('k',
                         zeros_like(self.k0),
                         1./logkSigma**2)
        
        #modeled species concentrations
        Dmod = pymc.Deterministic(self.forwardSoln,
                                  'modeled species concentrations', 
                                  'Dmod' , 
                                  {'k':k},
                                  verbose = 0, 
                                  plot = False)
       
        #observed species concentrations                                  
        Dobs = pymc.Stochastic(self.modelLogProbability,
                               'measured species concentrations',
                               'Dobs',
                               {'Dmod': Dmod},
                               value = self.D,
                               verbose = 0,
                               plot = False,
                               observed = True)
        
        
                                  
        self.mcModel = pymc.Model([k, Dmod, Dobs])
        self.mcmc = pymc.MCMC(self.mcModel)
        
    def k(self):
        for p in self.pcMod.parameters:
            if p[0]=='k':
                self.k_[int(p[1])-1] = self.pcMod.__dict__[p]            
        return self.k_
        
    def set_k(self, K):
        for i, k in enumerate(K):
            self.pcMod.__dict__['k%d'%(i+1)] = k
            
    def logK(self):
        return log2(self.k()/self.k0)
        
    def setLogK(self, LK):
        self.set_k(self.k0*(2**LK) )

    def modelLogProbability(self,value, Dmod):
#        x = matrix((Dmod-self.D).flatten()).T
#        logp = float64(self.logNormConst-.5*x.T*self.sigmaInv*x)

        r = sum((Dmod-self.D)**2)
        
        logp =  self.logNormConst+ -.5*r/self.sigma**2
        
#        print logp
#        print self.cnt
        
        return logp
    
    def forwardSoln(self, k):
        self.setLogK(k)     
        self.pcMod.Simulate(verbose = self.verbose)
        
        S = self.pcMod.data_sim.getSpecies()[:,1:]
#        print 'called!'
        return S
        
    
        
if __name__== "__main__":
    #nr samples to generate    
    N = 10000
    burn = 5000
    thin = 5
    
    force = False
    fname = 'testLogInverseModel1.pkl'
    
    #initialize the pysces model
    pcMod = pysces.model('pysces_test_linear2') 
    pcMod.sim_end = 5
    pcMod.sim_points = 60

    #initialize "data"
    pcMod.Simulate()
    D = pcMod.data_sim.getSpecies()[:,1:] 
    D +=sigmaExp*randn(D.shape[0], D.shape[1])
    
    #get initial parameters
    nK = 0
    for p in pcMod.parameters:
        if p[0]=='k':
            nK+=1
    k0 = zeros(nK)     
    for p in pcMod.parameters:
            if p[0]=='k':
                k0[int(p[1])-1] = pcMod.__dict__[p]
    
    
    if not force:
        try:
            with open(fname, 'rb') as f:
                [T,TD, logp, logNormConst] = pkl.load(f)
        except:
            force = True
        
    if force:            
        #initialize inverse model
        invMod = pyscesInverseModel(pcMod, D, verbose = False)
        
        #sample from the posterior    
        start = time.time()    
        invMod(N, burn = burn, thin = thin)
        elapsed = time.time()-start
        print 'sampled %d in %.2f s: %.2e samples/s'%(N,elapsed,N/elapsed  )
            
        T = invMod.mcmc.trace('k')[:].T
        TD = invMod.mcmc.trace('Dmod')[:]
        logp = array([invMod.modelLogProbability(None, TD[i]) for i in range(TD.shape[0])])
        logNormConst = invMod.logNormConst
        with open(fname, 'wb') as f:
            pkl.dump((T, TD, logp, logNormConst), f, -1)
        
    

    
    figure()
    x = linspace(-3*logkSigma, 3*logkSigma, 200)
    for i in range(T.shape[0]):
        subplot(3,3,i+1)
        hist(T[i],50, label = 'posterior', normed = True)
        plot(x, stats.norm.pdf(x,loc = 0, scale = logkSigma ), label = 'prior')
        title('k%d, truth = %.1f'%(i+1, k0[i]))
        xlabel(r'log(k/k$_{true}$)')
        
        
    
    didx = randint(0,N, 4)
    print "randomly selected trace points are:"    
    print didx
    figure()
    for i, idx in enumerate(didx):
        subplot(2,2,i+1) 
        plot(D, 'k')
        plot(TD[idx], 'y')
        title('log(P(D|m)) = %.2f'%logp[idx])
    
    figure()
    pcolor(arange(1.,T.shape[0]+1)-.5, 
           arange(1.,T.shape[0]+1)-.5,
           corrcoef(T))
    colorbar()
    xlabel('constant index')
    ylabel('constant index')
    title('correlation coefficient between rate constants k1-k8')
    
    try:
        figure()
        k2 = T[1]
        k3 = T[2]
        k4 = T[3]
        k8 = T[7]
        subplot(1,2,1)
        hist(k3[(abs(k4)<.05)], 50)
        xlabel('k3')
        title('Conditional distribution \n of k3: P(k3|k4=1)')
        subplot(1,2,2)
        hist(k4[abs(k3)<.05], 50)
        xlabel('k4')
        title('Conditional distribution \n of k4: P(k4|k3=5)')       
    except:
        pass
    
    
    
    
