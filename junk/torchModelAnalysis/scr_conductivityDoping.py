# -*- coding: utf-8 -*-
"""
Created on Wed Dec 28 12:35:34 2011

@author: wichtelwesen
"""
from numpy import *
from pylab import *
A = 100.
E0 = 10000.
tau = 200.
m = 1./100.

E = linspace(1000, 15000, 3000)
S_sigmoid = A/(1+exp(-(E-E0)/tau))
S_lin = m*(E-E0)
S = S_sigmoid
S[E>E0]+=S_lin[E>E0]

plot(E, S)
show()

#100/(1+exp(-(E.internalField()[celli]-10000)/200))
#if (E.internalField()[celli]>10000)
#   sigmaCells[celli]+=(E.internalField()[celli]-10000)/100