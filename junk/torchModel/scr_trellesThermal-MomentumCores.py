
from numpy import *
from trellesHelpers import fetchMesh, fetchVar, varDict
import mlabMacros as mlm
from pylab import *
from qhullCommands import delaunay
from scipy.interpolate import Rbf


########## location of the data#####################
meshDir = r"/home/wichtelwesen/PHOTOS/00PLASMA_DATA/juansMatlabData"
varDir = r"/home/wichtelwesen/PHOTOS/00PLASMA_DATA/trellesTorch"
########## location of the data#####################

###macro to summon up flow variables###
def fv(v):
    return fetchVar(v, outputDir = varDir, precision = 'double')
#######################################

#retrieve computational mesh
X = fetchMesh(outputDir = meshDir, precision = 'double')[0]

#indices corresponding to the torch outlet
outIdx = flatnonzero(X[2]>.03)
x = X[:2, outIdx]

#retrieve the temperatures and the velocity at the outlet
t, T_ = fv(varDict['heavyTemperature'])
t, Te_ = fv(varDict['electronTemperature'])
U_ = zeros(T_.shape+(3,))
for i in range(3):
    t, U_[:,:,i] = fv(varDict['velocity'][i])
    
#time averaged outlet quantities
Tbar = mean(T_[:,outIdx], axis = 0)
Tebar = mean(Te_[:,outIdx], axis = 0)
thetaBar = mean((Te_/T_)[:,outIdx], axis = 0)
Ubar = mean(U_[:,outIdx, 2], axis = 0)

#radial coordinate for the nodes
r = sum(X[:2, outIdx]**2, axis = 0)**.5

#make a quick scatter plot of T and U
figure()
plot(r, Tbar/mean(Tbar), 'r.', alpha = .5, label = 'temperature/meanTemp')
plot(r, Ubar/mean(Ubar), 'b.', alpha = .25, label = 'velocity/meanVel')
legend()

figure()
plot(r, thetaBar, 'b.', alpha = .3, label = 'electron temperature/ heavy temperature)')
legend()

G = mgrid[amin(x[0]):amax(x[0]):200j,
          amin(x[1]):amax(x[1]):200j]
R = sum(G**2, axis = 0)**.5

T = Rbf(x[0], x[1], Tbar, function = 'cubic')(G[0], G[1]); T[R>amax(r)] = mean(Tbar[r==amax(r)])
U = Rbf(x[0], x[1], Ubar,function = 'cubic')(G[0], G[1]); U[R>amax(r)] = mean(Ubar[r==amax(r)])
Te = Rbf(x[0], x[1], Tbar*thetaBar,function = 'cubic')(G[0], G[1]); Te[R>amax(r)] = mean((Tbar*thetaBar)[r==amax(r)])

eps = 1e-3
dT = (gradient(T)[0]**2+gradient(T)[1]**2)**.5; dT[R>(amax(r)-eps)] = 0
dU = (gradient(U)[0]**2+gradient(U)[1]**2)**.5; dU[R>(amax(r)-eps)] = 0




figure()
contourf(dT, 30, cmap = cm.copper); colorbar()
contour(dU, 10, cmap = cm.cool ); colorbar()


mlm.surf3(G[0], G[1], U, f = mlm.fig('axial velocity'))
mlm.surf3(G[0], G[1], T, f = mlm.fig('temperature'))
mlm.surf3(G[0], G[1], Te, f = mlm.fig('electron temperature'))

show()


 
 



    
    