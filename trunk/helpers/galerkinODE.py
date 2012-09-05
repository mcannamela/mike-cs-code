# -*- coding: utf-8 -*-
"""
Created on Tue Apr 05 10:52:08 2011

@author: wichtelwesen
"""

from numpy import *

class odeBasis(object):
    def __init__(self, t0 = 0, tN = 1, N = 100):        
        self.dt = (tN-t0)/N
        self.t = t0+arange(N)*self.dt
        self.tStar = self.t[:-1] + diff(self.t)
    
    def __call__(self, a):
        pass
    
    def ddt(self, a):
        pass

class firstOrderLagrangeBasis(object):
    def __init__(self, t0 = 0, tN = 1, N = 100):
        odeBasis.__init__(self, t0, tN, N)
        self.psi = zeros_like(self.tStar)
        self.psiDot = zeros_like(self.tStar)
        
    def __call__(self, a):
        self.psi = .5*(a[:-1]+a[1:])
        return self.psi.copy()
    
    def ddt(self, a):
        self.psiDot = (a[1:]-a[:-1])/self.dt
        return self.psiDot.copy()
        
class firstOrderExplicitGalerkinODE(object):
    def __init__(self, diffOp, phi, N):
        """
        represent ODE Galerkin style, which must have the form
        
        d(psi)/dt = diffOp(psi, t)
        
        diffOp - a function of the solution psi[i](t) 
        phi - the basis functions, must be of type 
                firstOrderLagrangeBasis
        N - dimension of 
        """
        self.a = zeros((N, phi.t.shape[0]))
        self.R = zeros(phi.t.shape[0])
        self.diffOp = diffOp
        self.phi = phi
        
    def __call__(self, a, psi0):
        self.a = a        
        self.a[:, 0] = psi0
        self.R = self.phi.ddt(a) - self.diffOp(self.phi(a))
        return self.R
        