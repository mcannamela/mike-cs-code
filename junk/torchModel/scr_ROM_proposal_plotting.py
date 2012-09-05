import os 
import scipy
from numpy import *
import numpy as np
import cPickle as pkl
from plotMacros import *
import mlabMacros as mlm
import trellesHelpers as thp
import qhullCommands as qh
pl = mlm.ml.pipeline
ml = mlm.ml

var = 'Th'
N = 400
x,el,bx,bel = thp.fetchMesh(elements = True, boundaries = True)
t,y = thp.fetchVar(var, N = N)
U,s,Vh,yBar = thp.fetchSVD(var,N = N, compute = False)


bidx = int32([])
blab = int32([])
for i in range(len(bx)):
	bidx = r_[bidx, bx[i]]
	blab = r_[blab, ones_like(bx[i])*i]

xend = x[:-1,bx[-2]]*1e3
mend = Vh[0,bx[-2]]

cidx = bx[-2][argmin(sum(xend**2, axis = 0))]

ROM = lambda i_: yBar+array(matrix(U[:,:i_])*matrix(diag(s[:i_]))*matrix(Vh[:i_,:]))
coeffs = lambda i_: array(matrix(U[:,:i_])*matrix(diag(s[:i_])))
soln = y[:,cidx]
approx1 = ROM(5)[:,cidx]
approx10 = ROM(10)[:,cidx]
approx30 = ROM(20)[:,cidx]
#
fig()
plot(t*1e3,soln, linewidth = 3, color = 'k', label = 'Exact Solution')
plot(t*1e3,approx1, linewidth = 2, color = 'b', alpha = 1, label = '5-term approx.')
plot(t*1e3,approx10, linewidth = 2, color = 'g',alpha = 1, label = '10-term approx')
plot(t*1e3,approx30, linewidth = 2, color = 'r', alpha = 1,label = '20-term approx')
xlab('time, ms')
ylab('Temperature, K')
tit('Solution at Outlet Centerline with SVD Approximations')
axisFontSize()
legend(loc = 'lower right')
sf('OutletCenterSolution')

fig()
n = r_[arange(3), arange(5,50, 5)]
e = zeros_like(n)
for i in range(n.shape[0]):
    e[i] = mean((ROM(n[i])[:,x[2]>.0025]-y[:,x[2]>.0025])**2)**.5

plot(n,e, linewidth = 3, color = 'k', label = 'Error')
xlab('approximation order')
ylab('Mean Error, K')
tit('Decay of Error with Increasing Approximation Order\nmax order = %d'%N)
axisFontSize()
sf('errorVsOrder')



fig()
r = range(5)
C = coeffs(r[-1]+1)
for i in range(5):
    plot(t*1e3, C[:,i], linewidth = 2, alpha = 1, color = cm.copper(int32(i*256./r[-1]))[:3])
    print cm.jet(int32(i*256./r[-1]))[:3]
xlab('time, ms')
ylab('Modal Amplitude')
tit('First %d Modal Coefficients Over Time\nbrighter = higher mode number'%(1+r[-1]))
axisFontSize()

sf('modalCoeffs')




show()