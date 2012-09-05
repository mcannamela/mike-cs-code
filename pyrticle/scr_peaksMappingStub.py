#scr_peaksMappingStub
import cPickle as pkl
import helpers as hp
from interpolators import cubicInterpolator as interpolator
import gasMixture as gm
import numpy as np
import scipy
scipy.pkgload()
import scipy.stats as stat
import mlabMacros as mlm
from plotMacros import *
from pylab import *



with open('crinkle-SharpT.pkl', 'rb') as f:
    PE = pkl.load(f)
    
name = 'sharpTemp_600_38.pkl'
with open(name, 'rb') as f: 
            J=pkl.load(f)
            
            
D = PE.parameterVectors[...,1]
V = abs(PE.parameterVectors[..., -2])
d = D[:,0]
v = V[0]

z = J['axialCoord']
r = J['radialCoord']

S = PE.summary(plotOn = False)

with open('twoPeaksMappingSharpT.pkl','wb') as f:
    pkl.dump(dict(zip(['D','V','d','v','z','r','S'],[D,V,d,v,z,r,S])), f, -1)