from numpy import *
seterr(invalid='raise')
seterr(divide='raise')
import time
import particleSolver as ps
import particle as part

import mlabMacros as mlm
import pdb
from plasmaField import *
from pylab import *
from scipy import integrate
import scipy.stats.distributions as dist
from carrierGasParameters import *
import cPickle as pkl
import warnings
warnings.simplefilter('error')

def plotOverSpace(X, F, skip):
    N = double(len(X[0]))
    for i in range(len(X[0])):
        if mod(i,skip)==0:        
            plot(F[i], '.-', color = [(N-i)/N, (N-i)/N, (N-i)/N], 
                 linewidth = 2,alpha = .5, label = str(i))
                 
def meanify(X):
    xBar = mean(X, axis = 0)[newaxis,...]
    return r_['0', xBar, xBar, xBar, xBar]
            
    
xlab = lambda s_: xlabel(s_, fontsize = 24)
ylab = lambda s_: ylabel(s_, fontsize = 24)
tit = lambda s_: title(s_, fontsize = 20)
sf = lambda figName_:  savefig(figName_+'.png', format = 'png', transparent = True)
fig = lambda:  figure(figsize = (20,22))

start = time.time()

#create the plasmaField
with open('/media/raidArray/CODE/jetModel/jet_arPilot_40_500_highDensity.pkl') as f:
    PF_t = pkl.load(f)

PF = clonedOpenFoamPlasma(meanify(PF_t.T), meanify(PF_t.U), meanify(PF_t.f), PF_t.X, linspace(0,1,4), PF_t.gas )
del PF_t

N = 10000
Qc = 7
d = exp(dist.norm(loc = -10.2, scale = .3).rvs(N))
v0 = dist.norm(loc = carrierGasParameters(Qc)[0], scale = carrierGasParameters(Qc)[1]).rvs(N)
while any(v0<0) :
    v0[v0<0]= dist.norm(loc = carrierGasParameters(Qc)[0], scale = carrierGasParameters(Qc)[1]).rvs(sum(v0<0))
    
t0 = dist.uniform(loc = PF.t[0], scale = PF.t[-1]).rvs(N)
theta = dist.uniform(loc = -15.,scale = 15.).rvs(N)*pi/180.
nNodes = 20

p = part.timeDependentJetParticle(      diameter=array([1e-9, d[0]]),
                        x0=float32([0.,0.,0.]),
                        v0=float32([0.,0.,0.]),
                        T=float32(300.), 
                        nNodes = nNodes,
                        thermalTimestep = float32(1.e-6),
                        kinematicTimestep = float32(1.e-6),
                        nStep = uint16(2000),
                        solidus_liquidus = float32([2912.,  3062.]),
                        vaporizationTemperature = float32(4600.),
                        heatOfVaporization = 5227000., 
                        vaporPressureFun = lambda T_: 101325.*exp(-4.5176e4/T_+10.088),
                        materialName = 'zirconiaGauss',
                        material = None,
                        plasmaField = PF,
                        t0=0.0)

p.adaptiveTimestepping=True
p.splitScales = 0

p.thermalStep   = float32(200)
p.inject(PF, array([0,  0.000,.004 ]),array([0, 0, -v0[0]]))
PS = ps.explicitParticleSolver(p, xMax = .11)
PS.alpha = .6
PS.stepErrorThres = 2e-5
PS.maxIter = 200
T = zeros( (v0.shape[0],nNodes))
r = zeros_like(T)
h = zeros(v0.shape[0])
timeHistory = []
v = zeros_like(h)
X = []
VMF = []
Tp = []
Ts = []
Tc = []
TsMinus = []
SS = []
R = []
failFlag = []
xEnd = .1
print 'init time: ' +str(time.time()-start)

start = time.time()
for i in range(v0.shape[0]):
    p = part.timeDependentJetParticle(      diameter=array([1e-7, d[i]]),
                        x0=float32([0.,.004,0.]),
                        v0=v0[i]*array([sin(pi*10/180), -cos(theta[i]), sin(theta[i])  ]),
                        T=float32(300.), 
                        nNodes = nNodes,
                        thermalTimestep = float32(1.e-6),
                        kinematicTimestep = float32(1.e-6),
                        nStep = uint16(2000),
                        solidus_liquidus = float32([2912.,  3062.]),
                        vaporizationTemperature = float32(4600.),
                        heatOfVaporization = 5227000., 
                        vaporPressureFun = lambda T_: 101325.*exp(-4.5176e4/T_+10.088),
                        materialName = 'zirconiaGauss',
                        material = None,
                        plasmaField = PF,
                        t0 = t0[i])

    p.adaptiveTimestepping=True
    p.splitScales = 0

    p.thermalStep   = float32(200)
    p.inject(PF, array([0,  0.004,.000 ]),v0[i]*array([sin(pi*10/180), -cos(theta[i]), sin(theta[i])]) )
    PS.p = p
    PS.reset(array([0,  0.004,.000 ]),v0[i]*array([sin(pi*10/180), -cos(theta[i]), sin(theta[i])]), t0 = t0[i] )
#    try:
    PS.solve()
    failFlag += [False]
    print 'solved!'
    
    r[i] = PS.getSolnCoords()[1][-1]
    R+=[PS.getSolnCoords()[1][:,-1]]
    h[i] = PS.getEnthalpy()
    timeHistory += [array(PS.p.t[:(PS.p.stepCount-1)]).copy()]
    X+=[PS.getTrajectory()]
    VMF+= [array(PS.p.vmfHistory)]
    Tp+= [array(PS.p.plasmaTemperatureHistory)]
    v[i] = sqrt(sum(PS.p.velocity()**2))
    SS+=[PS.getSoln()]
    S = PS.getSoln()
    
    Ts+=[S[:,-1]]
    Tc+=[S[:,0]]
    
    
    j = flatnonzero(X[-1][0]<xEnd)[-1]
    x = X[-1][0][j]
    try:
        xp = X[-1][0][j+1]
    except IndexError:
        T[i] = S[j]
        continue
        

    T[i] = S[j]+(xEnd-x)*(S[j+1]-S[j])/(xp-x)
#    except ValueError:
#        failFlag += [True]
#        print 'failed!'
#    
    elapsed = time.time()-start
    if mod(i,100)==0:
        print "computed %d in %d min, at %d s/particle"%(i+1, elapsed/60.0, elapsed/(i+1) )

    
elapsed = time.time()-start
print "run time is %d s, that's %d s/particle"%(elapsed, elapsed/max([N,1])) 

with open('timeAveragedDistributionStudy_%.1f_%d.pkl'%(Qc, N), 'wb') as f:
    pkl.dump(dict(zip(['diameter', 'initialVelocity', 'initialTime', 'angle',
                       'finalRadius', 'outerRadius', 'finalEnthalpy', 'timeHistory',
                       'trajectory', 'vaporMassFlux', 'plasmaTemperature', 'finalVelocity', 
                       'solution', 'surfaceTemperature', 'centerTemperature', 'finalTemperature'],
                      [d, v0, t0, theta, 
                       r,R, h, timeHistory, 
                       X, VMF, Tp, v, 
                       SS, Ts, Tc, T])),f, -1)
                       
#f = mlm.fig('Temperature Profiles for Different Injection Velocities')
#mlm.mesh3(v0[:,newaxis]*ones_like(r), r, T, f = f, 
           #axeLab = ['injection velocity, m/s', 'radius, um', 'temperature, K'])
#T08 = []
#for i, Ti in enumerate(T):
    
close('all')
PF.showField(velocity = False, planeNormal = 'z', t = t0[1], withColor = True)

N = double(v0.shape[0])
cu = cm.copper
for i in range(v0.shape[0]):
    plot(X[i][0], X[i][1],'.-', color = cu(double(i/N)),#'k',#[(N-i)/N, (N-i)/N, (N-i)/N], 
         linewidth = 2,alpha = .6, label = str(v0[i]))

xlab('x, m')
ylab('y, m')
tit('Temperature field with y trajectories')
xlim((0,.1))

PF.showField(velocity = False, planeNormal = 'y', t = t0[1], withColor = True)
for i in range(v0.shape[0]):
    plot(X[i][0], X[i][2],'.-', color = cu(double(i/N)),#'k',#[(N-i)/N, (N-i)/N, (N-i)/N], 
         linewidth = 2,alpha = .6, label = str(v0[i]))
xlab('x, m')
ylab('z, m')
tit('Temperature field with z trajectories')
xlim((0,.1))

#figure()
#for i in range(v0.shape[0]):
#    plot(X[i][0], VMF[i],'.-', color = [(N-i)/N, (N-i)/N, (N-i)/N], 
#         linewidth = 2,alpha = .5, label = str(v0[i]))
#xlab('x, m')
#ylab('vapor mass flux, kg/s')
#tit('vapor mass flux')

figure()
for i in range(v0.shape[0]):
    plot(X[i][0], Tp[i],':', color = cm.Blues(double(i/N)),#[(N-i)/N, (N-i)/N, (N-i)/N], 
         linewidth = 2,alpha = .8, label = str(v0[i]))
    plot(X[i][0], Ts[i],'--', color = cm.Reds(double(i/N)),#,[(N-i)/N, (N-i)/N, (N-i)/N], 
         linewidth = 2,alpha = .8, label = str(v0[i]))
    plot(X[i][0], Tc[i],'-', color = cm.Greens(double(i/N)),#[1-(N-i)/N, (N-i)/N, (N-i)/N], 
         linewidth = 2,alpha = .8, label = str(v0[i]))
xlab('x, m')
ylab('Temperature, K')
tit('Plasma and particle surface/center temperature histories')

figure()
for i in range(v0.shape[0]):
    plot(r[i], T[i],'.-', color = cu(double(i/N)),#[(N-i)/N, (N-i)/N, (N-i)/N], 
         linewidth = 2,alpha = 1, label = str(v0[i]))
xlab(r'radius, $\mu$m')
ylab('Temperature, K')
tit('Particle temperature profiles at x = 10 cm')



#figure()
#plot(v0, T[:,-1],'ko-', linewidth = 2)
#xlabel('injection velocity, m/s')
#ylabel('temperature, K')
#title('surface temperature vs injection velocity')
#
#figure()
#plot(v0, h,'ko-', linewidth = 2)
#xlabel('injection velocity, m/s')
#ylabel('enthalpy, J/kg')
#title('net enthalpy aquired vs injection velocity')
#
#figure()
#plot(v0,[t[-1] for t in timeHistory],'ko-', linewidth = 2)
#xlabel('injection velocity, m/s')
#ylabel('travel time, s')
#title('travel time vs injection velocity')
#
#figure()
#plot(v0, v, 'ko-', linewidth = 2)
#xlabel('injection velocity, m/s')
#ylabel('final speed, m/s')
#title('final speed vs injection velocity')
