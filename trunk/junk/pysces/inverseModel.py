# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 16:31:06 2011

@author: wichtelwesen
"""

from numpy import *
import cPickle as pkl

class forwardModel(object):
    def __init__(self, *args, **kwargs):
        """
        initialize the model and any data members needed to drive it
        """
        pass
    def __call__(self, modelParameters):
        """
        evaluate the model specified by modelParameters, 
        return the model output
        """
        pass
    def black(self):
        """
        build the blacklist comprised of unpickleable data members
        """
        B = []        
        for k in self.__dict__.keys():
            try:
                pkl.dumps(self.__dict__[k])
            except:
                B+= [k]
        self.pklBlackList = B
        return B
    
    def m(self):
        """
        return the current model parameters
        """
        pass
    
    def set_m(self, modelParameters):
        """
        set the model parameters according to the argument modelParameters
        """
        pass
    
    
class inverseModel(object):
    def __init__(self, forwardModel, observedData):
        """
        initialize with a forwardModel object and some observed data with which 
        to perform the inversion
        """        
        self.FM = forwardModel
        self.D_obs = observedData
        
    def initMC(self):
        """
        initialize the mcmc model that will compute the posterior probabilities
        """
        pass
    
    
    def __call__(self, N, burn = 0, thin = 1):
        """
        sample the posterior distribution
        N - number of samples to generate
        burn - number of samples to discard before collecting N samples
        thin - keep only every thin samples once we have burned in to reduce 
                autocorrelation
        """
        print "N, burn, thin are %d, %d, %d"%(N, burn, thin)
        self.mcmc.sample(iter = burn+thin*N, burn = burn, thin = thin)
        self.trace()
        
    def trace(self):
        """
        record traces here for easy access
        """
        pass
    
    def modelLikelihood(self, value, D_mod):
        """
        the model likelihood function (see Tarantola)
        D_mod - the modeled data 
        value - initial value for D_mod, can be None
        """
        
        #by default, assume that D_obs is an object representing the observed
        #data and that it has a logp(self, other) method return the log(probability) 
        #that the data represented by D_mod are actually the same as D_obs
        return self.D_obs.logp(D_mod)
    