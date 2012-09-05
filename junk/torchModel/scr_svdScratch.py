import os 
import scipy
scipy.pkgload('io')
from numpy import *
import numpy as np
import cPickle as pkl
import time
from pylab import *
from plotMacros import *
import mlabMacros as mlm

outputDir = "D:\\00PLASMA_DATA\\trellesTorch"
tName = 'trellesTorch_'
VAR = 'Th'
N = 500
compute = False
if compute:
	with open(os.path.join(outputDir,tName+VAR+'.pkl'),'rb') as f:
		D = pkl.load(f)

	y = D['y'][:N].copy()
	t = D['t'][:N].copy()
	del D

	yBar = np.mean(y,axis = 0)
	y-= yBar

	U,s,Vh = linalg.svd(y, full_matrices = False)
	with open(os.path.join(outputDir,tName+VAR+'_svd_%d.pkl'%N),'wb') as f:
		pkl.dump(dict(zip(['U','s','Vh'],[U,s,Vh])), f, -1)
else:
	D = pkl.load(open(os.path.join(outputDir,tName+VAR+'_svd_%d.pkl'%N),'rb'))
	U,s,Vh = (D['U'].copy(),D['s'].copy(),D['Vh'].copy())
	del D
fig()
subplot(1,2,1)
plot(20*log10(s/s[0]), '.-',linewidth = 2, label = 'singular values')
plot(arange(s.shape[0]), -3*ones_like(s), 'r', label = '-3dB mark')
ylab('sigular value, dB ref s$_0$')
xlab('rank')
tit('SVD of '+ VAR); axisFontSize()


subplot(1,2,2)
plot(cumsum(s**2)/sum(s**2), '.-',linewidth = 2, label = 'fractional energy')
plot(arange(s.shape[0]), .9*ones_like(s), 'r',label = '90% mark')
ylab('fractional energy')
xlab('rank')
tit('Cumulative Energy in /nSingular Values of '+ VAR); axisFontSize()
sf('singularValuePlots')


